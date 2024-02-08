from django.contrib.auth.models import Group
from django.contrib.contenttypes.models import ContentType

from jobs.dbinterface.postgresinterface import PostgresInterface

from django.test import TestCase

from jobs.models import *
from pgmols.models import *

SCRATCH_DIR = "scratch"
WORKER_NAME = "test_worker"

SMILES_STR = "c1ccccc1"

FAKE_RESULT = {'meta_data': {'wall_clock_time': 60,
                             'num_procs': 8,
                             'program': 'fakesim',
                             'version': 1.0,
                             'host': 'localhost'},
               'coord_list': [{'element': 'C', 'x': 0, 'y': 4.0, 'z': 6},
                              {'element': 'C', 'x': 0, 'y': 5.0, 'z': 6},
                              {'element': 'C', 'x': 1, 'y': 6.0, 'z': 6},
                              {'element': 'C', 'x': 2, 'y': 7.4, 'z': 6},
                              {'element': 'C', 'x': 3, 'y': 8.1, 'z': 6},
                              {'element': 'C', 'x': 4, 'y': 9.2, 'z': 6}],
               'theory': {'name': 'fake_conformer',
                          'description': 'a pretend conformer generator'}
               }


class TestPostgresInterface(TestCase):

    def setUp(self):
        self.pgdbi = PostgresInterface()
        self.config = JobConfig(name="test_config")
        self.config.save()
        self.test_group = Group(name="test_group")
        self.test_group.save()

    def _make_mol_job(self):
        m = Mol(smiles=SMILES_STR)
        m.save()
        job = Job(parent=m, config=self.config, group=self.test_group)
        job.save()
        return m

    def _make_geom_job(self):
        method = Method(name="test_method")
        method.save()
        m = Mol(smiles=SMILES_STR, group=self.test_group)
        m.save()
        geom = Geom(mol=m, method=method, xyz=[[0, 0, 0, 0]])
        geom.save()
        job = Job(parent=geom, config=self.config, group=self.test_group)
        job.save()
        return geom

    def test_create_job_with_mol_parent(self):
        self._make_mol_job()
        self.assertEquals(Job.objects.all()[0].parent.smiles, SMILES_STR)

    def test_create_job_with_geom_parent(self):
        self._make_geom_job()
        self.assertEquals(Job.objects.all()[0].parent.mol.smiles, SMILES_STR)

    def test_find_one_job(self):
        self._make_geom_job()
        jobs = self.pgdbi.get_new_jobs('test_config', 'test_group')
        self.assertEquals(len(jobs), 1)

    def test_find_job_by_mol(self):
        mol = self._make_mol_job()
        jobs = Job.objects.filter(parentid=mol.id, parentct=ContentType.objects.get_for_model(mol))
        self.assertEquals(len(jobs), 1)

    def test_claim_job(self):
        self._make_mol_job()
        jobspec = self.pgdbi.get_new_jobs('test_config', 'test_group')[0]
        self.pgdbi.claim_job(jobspec)
        job = Job.objects.all()[0]
        self.assertEquals(job.status, 'claimed')

    def test_claim_job_twice_fails(self):
        self._make_mol_job()
        jobspec = self.pgdbi.get_new_jobs('test_config', 'test_group')[0]
        self.pgdbi.claim_job(jobspec)
        with self.assertRaises(Exception):
            self.pgdbi.claim_job(jobspec)

    def test_release_job_twice_with_ignore_succeeds(self):
        self._make_mol_job()
        jobspec = self.pgdbi.get_new_jobs('test_config', 'test_group')[0]
        self.pgdbi.claim_job(jobspec)
        self.pgdbi.release_job(jobspec)
        with self.assertRaises(Exception):
            self.pgdbi.release_job(jobspec)
        self.pgdbi.release_job(jobspec, ignore_lock=True)

    def test_release_job_with_errors(self):
        self._make_mol_job()
        jobspec = self.pgdbi.get_new_jobs('test_config', 'test_group')[0]
        self.pgdbi.claim_job(jobspec)
        self.pgdbi.release_job(jobspec, errors=True)
        job = Job.objects.all()[0]
        self.assertEquals(job.status, 'error')

    def test_reclaim_completed_job_fails(self):
        self._make_mol_job()
        jobspec = self.pgdbi.get_new_jobs('test_config', 'test_group')[0]
        self.pgdbi.claim_job(jobspec)
        self.pgdbi.release_job(jobspec)
        with self.assertRaises(Exception):
            self.pgdbi.claim_job(jobspec)

    def test_claim_and_released(self):
        self._make_mol_job()
        jobspec = self.pgdbi.get_new_jobs('test_config', 'test_group')[0]
        self.pgdbi.claim_job(jobspec)
        self.pgdbi.release_job(jobspec)
        job = Job.objects.all()[0]
        self.assertEquals(job.status, 'done')

    def test_unknown_mol_parent(self):
        mol = self._make_mol_job()
        jobspec = self.pgdbi.get_new_jobs('test_config', 'test_group')[0]
        self.pgdbi.claim_job(jobspec)
        job = Job.objects.all()[0]
        random_job = Job(parent=mol, config=self.config, group=self.test_group)
        job.parent = random_job
        job.save()
        with self.assertRaises(Exception):
            self.pgdbi.add_result(FAKE_RESULT, jobspec)

    def test_unknown_geom_parent(self):
        geom = self._make_geom_job()
        jobspec = self.pgdbi.get_new_jobs('test_config', 'test_group')[0]
        self.pgdbi.claim_job(jobspec)
        job = Job.objects.all()[0]
        random_job = Job(parent=geom, config=self.config, group=self.test_group)
        job.parent = random_job
        job.save()
        with self.assertRaises(Exception):
            self.pgdbi.add_result(FAKE_RESULT, jobspec)

    def test_unknown_calc_parent(self):
        result = {'properties': {'total_energy': 5,
                                 'number_list': [1, 2, 3]},
                  'velocity_list': [0.0, 1.0, 5.0]}

        result.update(FAKE_RESULT)
        del(result["coord_list"])
        geom = self._make_geom_job()
        jobspec = self.pgdbi.get_new_jobs('test_config', 'test_group')[0]
        self.pgdbi.claim_job(jobspec)
        job = Job.objects.all()[0]
        random_job = Job(parent=geom, config=self.config, group=self.test_group)
        job.parent = random_job
        job.save()
        with self.assertRaises(Exception):
            self.pgdbi.add_result(result, jobspec)

    def test_add_geom_no_calc_to_mol(self):
        mol = self._make_mol_job()
        jobspec = self.pgdbi.get_new_jobs('test_config', 'test_group')[0]
        self.pgdbi.claim_job(jobspec)
        self.pgdbi.add_result(FAKE_RESULT, jobspec)
        self.pgdbi.release_job(jobspec)
        geom = Geom.objects.get(mol=mol)
        job = Job.objects.all()[0]
        method = Method.objects.get(name='fake_conformer')
        self.assertEquals(geom.xyz[3][2], 7.4)
        self.assertEquals(job.duration, 60)
        self.assertEquals(geom.method, method)
        self.assertEquals(0, Calc.objects.count())

    def test_add_calc_below_geom_from_seperate_job(self):
        mol = self._make_mol_job()
        jobspec = self.pgdbi.get_new_jobs('test_config', 'test_group')[0]
        self.pgdbi.claim_job(jobspec)
        self.pgdbi.add_result(FAKE_RESULT, jobspec)
        self.pgdbi.release_job(jobspec)
        geom = Geom.objects.get(mol=mol)
        result = {'properties': {'total_energy': 5,
                                 'number_list': [1, 2, 3]},
                  'velocity_list': [0.0, 1.0, 5.0]}
        result.update(FAKE_RESULT)
        result.pop('coord_list')
        calc_job = Job(parent=geom, config=self.config, group=self.test_group)
        calc_job.save()
        jobspec2 = self.pgdbi.get_new_jobs('test_config', 'test_group')[0]
        self.pgdbi.claim_job(jobspec2)
        self.pgdbi.add_result(result, jobspec2)
        self.pgdbi.release_job(jobspec2)
        self.assertEquals(Calc.objects.count(), 1)
        calc = Calc.objects.first()
        self.assertEquals(calc.geoms.first(), geom)

    def test_add_geom_and_calc_to_mol(self):
        result = {'properties': {'total_energy': 5,
                                 'number_list': [1, 2, 3]},
                  'velocity_list': [0.0, 1.0, 5.0]}

        result.update(FAKE_RESULT)
        mol = self._make_mol_job()
        jobspec = self.pgdbi.get_new_jobs('test_config', 'test_group')[0]
        self.pgdbi.claim_job(jobspec)
        self.pgdbi.add_result(result, jobspec)
        self.pgdbi.release_job(jobspec)
        geom = Geom.objects.get(mol=mol)
        job = Job.objects.all()[0]
        method = Method.objects.get(name='fake_conformer')
        self.assertEquals(geom.xyz[3][2], 7.4)
        self.assertEquals(job.duration, 60)
        self.assertEquals(geom.method, method)
        calc = Calc.objects.get(geoms=geom)
        self.assertEquals(calc.props['totalenergy'], 5)
        self.assertEquals(calc.props['numbers'], [1, 2, 3])
        self.assertEquals(calc.props['velocities'], [0.0, 1.0, 5.0])

    def test_add_geom_and_calc_to_geom(self):
        result = {'properties': {'total_energy': 5}}
        result.update(FAKE_RESULT)
        root_geom = self._make_geom_job()
        jobspec = self.pgdbi.get_new_jobs('test_config', 'test_group')[0]
        self.pgdbi.claim_job(jobspec)
        self.pgdbi.add_result(result, jobspec)
        self.pgdbi.release_job(jobspec)
        job = Job.objects.all()[0]
        geom = Geom.objects.get(parentjob=job)
        method = Method.objects.get(name='fake_conformer')
        self.assertEquals(geom.xyz[3][2], 7.4)
        self.assertEquals(job.duration, 60)
        self.assertEquals(geom.method, method)
        calc = Calc.objects.get(geoms=geom)
        self.assertEquals(calc.props['totalenergy'], 5)
        self.assertEquals(Geom.objects.count(), 2)
        self.assertEquals(job.parent, root_geom)
        self.assertEquals(job, root_geom.childjob.first())
        self.assertEquals(job.childgeoms.first(), geom)

    def test_add_calc_to_calc(self):
        result = {'properties': {'total_energy': 5}}
        result.update(FAKE_RESULT)
        self._make_geom_job()
        jobspec = self.pgdbi.get_new_jobs('test_config', 'test_group')[0]
        self.pgdbi.claim_job(jobspec)
        self.pgdbi.add_result(result, jobspec)
        self.pgdbi.release_job(jobspec)

        del(result['coord_list'])
        job = Job.objects.all()[0]
        geom = Geom.objects.get(parentjob=job)
        calc = Calc.objects.get(geoms=geom)
        calc_job = Job(parent=calc, config=self.config, group=self.test_group)
        calc_job.save()
        calcjobspec = self.pgdbi.get_new_jobs('test_config', 'test_group')[0]
        self.pgdbi.claim_job(calcjobspec)
        self.pgdbi.add_result(result, calcjobspec)
        self.pgdbi.release_job(calcjobspec)
        calcjob = Job.objects.filter(parentid=calc.id, parentct=ContentType.objects.get_for_model(calc))[0]
        self.assertEquals(calcjob.parent, calc)
        self.assertEquals(Calc.objects.count(), 2)
