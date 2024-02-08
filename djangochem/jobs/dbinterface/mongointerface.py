from aag_python.molecular_storage import molecular_data_models as mdm
from .dbinterface import *


class MongoInterface(DatabaseInterface):
    def __init__(self, job_class_name='Calculation'):
        self.job_class_name = job_class_name

    def get_new_jobs(self, worker_name, project_name, min_priority=0, max_priority=10, limit=None):
        """
        returns calculation query objects containing calcs that this worker can work on
        """
        job_class = getattr(mdm, self.job_class_name)

        job_query = job_class.objects(meta_data__status__current_work_order=worker_name,
                                      meta_data__status__current_status__ne=mdm.Status.Value.CALC_RUNNING,
                                      meta_data__status__current_work_priority__gte=min_priority,
                                      meta_data__status__current_work_priority__lte=max_priority,
                                      meta_data__project_name=project_name
                                      ).order_by("meta_data.status.current_work_priority")
        if limit is not None:
           job_query = job_query.limit(limit)
        job_query = job_query.no_cache()
        return [self._make_job(p, worker_name, project_name) for p in job_query]

    def _make_job(self, parent, worker_name, project_name):
        """
        for parent object return a dictionary with key information for creating job
        """
        js = {}
        if hasattr(parent, 'coord_list'):
            js.update(dict(coords=parent.coord_list))
        confnum = None
        for tag in parent.meta_data.tag_list:
            if re.search("conf", tag):
                confnum = tag.split('_')[-1]
                break
        js[WORKER_NAME_KEY] = worker_name
        js[JOB_CLASS_KEY] = parent._class_name
        js[JOB_KEY] = parent.meta_data.key
        js[PARENT_CLASS_KEY] = parent._class_name
        js['priority'] = parent.meta_data.status.current_work_priority
        js['inchikey'] = parent.meta_data.inchi_key
        js[PARENT_KEY_KEY] = parent.meta_data.key
        js['smiles'] = parent.meta_data.smiles
        js[PROJECT_KEY] = project_name
        if confnum:
            js['details']['confnum'] = confnum
        return js

    def _get_job_parent(self, jobinfo):
        parent_class = getattr(mdm, jobinfo[PARENT_CLASS_KEY])
        return parent_class.objects.get(meta_data__key=jobinfo[PARENT_KEY_KEY])

    def confirm_exists(self, jobspec, claimed=True):
        job_class = getattr(mdm, jobspec[JOB_CLASS_KEY])
        query = job_class.objects.filter(meta_data__key=jobspec[JOB_KEY])
        if claimed:
            query = query.filter(meta_data__status__current_status=mdm.Status.Value.CALC_RUNNING)
        return query.count() == 1

    def claim_job(self, jobspec):
        job_class = getattr(mdm, jobspec[JOB_CLASS_KEY])
        query = job_class.objects(meta_data__key=jobspec[JOB_KEY],
                                  meta_data__status__current_status__ne=mdm.Status.Value.CALC_RUNNING)
        update_count = query.update(set__meta_data__status__current_status=mdm.Status.Value.CALC_RUNNING,
                                    push__meta_data__status__status_list=mdm.Status.Value.CALC_RUNNING)
        if update_count != 1:
            raise LockError("Unable to claim Job {}".format(JOB_KEY))

    def release_job(self, jobinfo, errors=False, force=False):
        if errors:
            new_status = mdm.Status.Value.CALC_FAILED
        else:
            new_status = mdm.Status.Value.CALC_DONE

        parent_class = getattr(mdm, jobinfo[PARENT_CLASS_KEY])
        key = jobinfo[PARENT_KEY_KEY]
        try:
            parent = parent_class.objects.get(meta_data__key=key)
        except:
            raise StandardError("Unable to locate parent object {0}.  lock release failed ".format(key))
        last_status = parent.meta_data.status.current_status

        if not force:
            assert(parent.meta_data.status.current_status == mdm.Status.Value.CALC_RUNNING)
        else:
            self.logger.warning("Ignored status {0} for database object with id {1}".format(last_status, key))
        parent.meta_data.status.current_work_order = ""
        parent.meta_data.status.work_history_list.append(jobinfo[WORKER_NAME_KEY])
        parent.meta_data.status.current_status = new_status
        parent.meta_data.status.status_list.append(new_status)
        parent.save()

    def add_result(self, result, jobinfo, suffix=None, force=False):
        parent = self._get_job_parent(jobinfo)
        worker_name = jobinfo[WORKER_NAME_KEY]

        calc = mdm.Calculation()
        for key in list(result.keys()):
            value = result.pop(key)
            if key == 'meta_data':
                temp = mdm.MetaData(**value)
                calc.meta_data = temp
            elif key == 'properties':
                temp = mdm.Properties(**value)
                calc.properties = temp
            elif key == 'theory':
                temp = mdm.Theory(theory_level=value['name'], theory_desction=value['description'])
                calc.theory = temp
            elif key == 'coords':
                calc.coord_list = value
            elif hasattr(calc, key):
                setattr(calc, key, value)
            else:
                msg = "Mongo Interface doesn't know where to put key {}".format(key)
                raise StandardError(msg)

        calc.meta_data.key = parent.meta_data.key + "_" + worker_name
        if suffix:
            calc.meta_data.key += "_" + suffix
        if calc.meta_data.status is None:
            calc.meta_data.status = mdm.Status(current_work_priority=parent.meta_data.status.current_work_priority)
        else:
            calc.meta_data.status.current_work_priority = parent.meta_data.status.current_work_priority
        calc.meta_data.project_name = parent.meta_data.project_name
        calc.validate()

        if type(parent) == mdm.Molecule:
            calc.molecule = parent
            calc.meta_data.generation = 0
            calc.save()
            parent.calculation_list.append(calc)
            parent.save()
        if type(parent) == mdm.Calculation:
            calc.parent_calculation_list.append(parent)
            calc.molecule = parent.molecule
            if len(calc.coord_list) > 0:
                calc.meta_data.generation = parent.meta_data.generation + 1
            else:
                calc.meta_data.generation = parent.meta_data.generation
            calc.save()
            parent.child_calculation_list.append(calc)
            parent.save()
        return calc.meta_data.key