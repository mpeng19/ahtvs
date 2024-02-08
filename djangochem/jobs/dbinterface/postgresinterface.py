import uuid
from django.utils import timezone
from django.db import transaction
from ..models import Job, JobConfig
from pgmols.models import Method, Geom, Calc, Mol
from .dbinterface import DatabaseInterface, WORKER_NAME_KEY, JOB_KEY, PROJECT_KEY, LockError, JobNotFoundError

CLAIMED = 'claimed'
DONE = 'done'
ERROR = 'error'
CLEAR = ''

PROP_MIGRATE = {"orbital_energy_list": "orbitals",
                "normal_mode_list": "normalmodes",
                "excited_states": "excitedstates",
                "velocity_list": "velocities"}


class PostgresInterface(DatabaseInterface):

    def get_new_jobs(self, config_name, project_name, min_priority=0,
                     max_priority=10, limit=None, extraquery={},
                     job_to_jobspec=None):
        """
        returns calculation query objects containing calcs that this worker can work on
        """
        with transaction.atomic():
            query = Job.objects.filter(group__name=project_name,
                                       config__name=config_name,
                                       status=CLEAR,
                                       **extraquery)
            if min_priority > 0:
                query = query.filter(priority__gte=min_priority)
            if max_priority < 10:
                query = query.filter(priority__lte=max_priority)
            if limit is not None:
                query = query.all()[:limit]
        if not job_to_jobspec:
            job_to_jobspec = self._make_jobspec
        return [job_to_jobspec(job) for job in query]

    def _make_jobspec(self, job):
        js = {}
        try:
            js.update(dict(coords=job.parent.get_coords()))
        except AttributeError:
            pass
        js[WORKER_NAME_KEY] = job.config.name
        try:
            mol = job.parent.mol
        except AttributeError:
            mol = job.parent

        js[JOB_KEY] = job.id
        js['uuid'] = str(job.uuid)
        js['priority'] = job.priority
        js['inchikey'] = mol.inchikey
        js['smiles'] = mol.smiles
        js['details'] = job.details
        js[PROJECT_KEY] = job.group.name
        js['charge'] = mol.molecular_charge
        js['parent_class'] = job.config.parent_class_name
        return js

    def claim_job(self, jobspec, extra_details={}):
        """ claiming a job prevents other interfaces from allowing it get claimed"""
        try:
            with transaction.atomic():
                job = Job.objects.select_for_update().filter(uuid=jobspec['uuid'], status__in=[CLEAR, ERROR]).get()
                job.status=CLAIMED
                if extra_details:
                    job.details.update(extra_details)
                job.save()
        except Job.DoesNotExist:
            raise LockError("Unable to claim Job {}".format(jobspec[JOB_KEY]))
        #     raise LockError("Unable to claim Job {}".format(jobspec[JOB_KEY]))
        # if 1 != job_query.update(status=CLAIMED):
        #     raise LockError("Unable to claim Job {}".format(jobspec[JOB_KEY]))
        # elif extra_details:
        #     raise NotImplementedError
        #     job.refresh_from_db()  # update to include the status CLAIMED
        #     job.details.update(extra_details)
        #     job.save()

    def clear_job(self, jobspec):
        with transaction.atomic():
            job = Job.objects.select_for_update().filter(uuid=jobspec['uuid'], status__in=[CLAIMED, ERROR]).get()
            job.status = CLAIMED
            job.save()

    def confirm_exists(self, jobspec, claimed=True):
        """ Just run a job query to see if the job is in the db and claimed"""
        try:
            self._query_job(jobspec, force=False, claimed=True)
            return True
        except:
            return False

    def release_job(self, jobspec, errors=False, force=False):
        job = self._query_job(jobspec, force, claimed=not force)
        if errors:
            job.status = ERROR
        else:
            job.status = DONE
        job.save()

    def add_result(self, result, jobspec, suffix=None, force=False):
        # if coords is there, make a geom
        method, created = Method.objects.get_or_create(**result['theory'])
        job = self._update_job(result, jobspec, force)

        geom = self._make_geom(result, job, method)
        self._make_calc(result, geom, job, method)

    def _make_props(self, result):
        props = {}
        if 'properties' in result:
            for k, v in result['properties'].items():
                if v is not None and not k.startswith('_'):
                    new_key = k.replace('_', '')
                    if new_key.endswith('list'):
                        new_key = new_key[:-4] + 's'
                    props[new_key] = v
        for attr, pg_attr in PROP_MIGRATE.items():
            val = result.get(attr)
            if val:
                props[pg_attr] = val
        return props


    def _query_job(self, jobspec, force, claimed=False):
        try:
            if claimed:
                job = Job.objects.get(uuid=jobspec['uuid'], status=CLAIMED)
            else:
                job = Job.objects.get(uuid=jobspec['uuid'])
        except (KeyError, Job.DoesNotExist) as err:  # handle Mol jobs created without uuid or with wrong uuid
            mol_parent = jobspec.get('parent_class') == 'Mol' or JobConfig.objects.get(name=jobspec[WORKER_NAME_KEY]).parent_class_name == 'Mol'
            if force and mol_parent:
                # go try and find the parent mol and get it's job
                mol = Mol.objects.get(group__name=jobspec[PROJECT_KEY], inchikey=jobspec['inchikey'])
                if len(mol.childjobs.all()) == 1:
                    job = mol.childjobs.all()[0]
                    if 'uuid' not in jobspec:
                        jobspec['uuid'] = str(uuid.uuid4())
                        job.uuid = jobspec['uuid']
                        job.save()
                else:
                    raise Exception("Attempted force find of Job, {} found".format(len(mol.childjobs.all())))
            else:
                if type(err) == Job.DoesNotExist:
                    if claimed:
                        raise JobNotFoundError(msg="No CLAIMED job found with uuid {}.  maybe try to force with -f".format(jobspec['uuid']))
                    else:
                        raise JobNotFoundError(msg="No job found with uuid {}".format(jobspec['uuid']))
                else:
                    raise
        return job

    def _update_job(self, result, jobspec, force=False):
        job = self._query_job(jobspec, force)
        job_details = {'software':result['meta_data'].get('program'),
                           'version':result['meta_data'].get('version'),
                           'nodename':result['meta_data'].get('host'),
                           **result['meta_data']
        }
        if job.details is None:
            job.details = {}
        job.details.update(job_details)
        job.duration = result['meta_data'].get('wall_clock_time')
        job.cores = result['meta_data'].get('num_procs')
        job.completetime = timezone.now()
        job.save()
        return job

    def _make_geom(self, result, job, method):
        geom = None
        if 'coords' in result:
            geom = Geom(method=method, parentjob=job)
            if 'details' in result:
                geom.details = result['details']
            geom.set_coords(result['coords'])
            if type(job.parent) == Mol:
                geom.mol = job.parent
            elif type(job.parent) == Geom:
                geom.mol = job.parent.mol
            else:
                raise Exception("Don't know how to handle job parent type {}".format(type(job.parent)))
            geom.save()
            if type(job.parent) == Geom:
                geom.parents.add(job.parent)
        return geom

    def _make_calc(self, result, geom, job, method):
        props = self._make_props(result)
        if props:
            calc = Calc(method=method, parentjob=job)
            if type(job.parent) not in [Mol, Calc, Geom]:
                raise Exception("Don't know how to handle job parent type {}".format(type(job.parent)))
            if type(job.parent) == Mol:
                calc.mol = job.parent
            else:
                calc.mol = job.parent.mol
            calc.props = props
            calc.save()
            if type(job.parent) == Calc:
                calc.parents.add(job.parent)
            if geom:
                calc.geoms.add(geom)
            elif type(job.parent) == Geom:
                calc.geoms.add(job.parent)

