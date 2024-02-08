import json
import os
import unittest

try:
    import boto3
    BOTO_AVAIL = True
    from jobs.storage.s3storage import S3WriteDir, S3ReadDir
except ImportError:
    BOTO_AVAIL = False

from django.test import TestCase

from jobs.models import *
from pgmols.models import *

from jobs.dbinterface import postgresinterface

from jobs.jobdirbuilder import JobDirBuilder
from jobs.jobdirparser import JobDirParser

MY_DIR = os.path.dirname(__file__)
CONFIG_PATH = os.path.join(MY_DIR, "config")



WORKER_NAME = "test_worker"

SMILES_A = "c1ccccc1"
SMILES_B = "c1cc([Mg])ncc1"
SMILES_C = "c1nc([Mg])cc([Mg])c1"
SMILES_D = "[Mg][Cl]"


WORKER_NAME_KEY = "worker_name"
PARENT_KEY_KEY = "parent_key"
PARENT_CLASS_KEY = "parent_class"
PROJECT_NAME = "test_project"
BATCH_KEY = "batch"
DEFAULT_BUCKET_NAME = "a2g2.unittests"

JOB_TEMPLATE_FILENAME = "test.inp"
INBOX_PATH = "inbox"

@unittest.skipIf(not BOTO_AVAIL, reason="requires boto3")
class TestPostgresS3JobDirHandlers(TestCase):
    def _db_interface(self):
        return postgresinterface.PostgresInterface()

    def setUp(self):
        self.db_interface = self._db_interface()
        self.dirbuilder = JobDirBuilder(name="test_worker",
                                        project=PROJECT_NAME,
                                        config_path=CONFIG_PATH,
                                        db_interface=self.db_interface,
                                        job_filename=JOB_TEMPLATE_FILENAME,
                                        template_filenames=[],
                                        storage_class=S3WriteDir,
                                        storage_kwargs={'bucket_name': DEFAULT_BUCKET_NAME}
                                        )
        self.dirbuilder._make_batch_key_name = self._mock_time_stamper

        s3 = boto3.resource('s3')
        self.bucket = s3.Bucket(DEFAULT_BUCKET_NAME)
        self.config = JobConfig(name="test_worker")
        self.config.save()
        self.test_group = Group(name=PROJECT_NAME)
        self.test_group.save()

    def _mock_time_stamper(self):
        return "batch"

    def tearDown(self):
        for key in self.bucket.objects.all():
            key.delete()

    def _make_mol_job(self, smiles=SMILES_A, priority=0):
        mol = Mol(smiles=smiles)
        mol.save()
        job = Job(parent=mol,
                  config=self.config,
                  group=self.test_group,
                  priority=priority)
        job.save()
        return mol, job


    def test_find(self):
        self._make_mol_job()
        self.assertEqual(1, len(self.dirbuilder.db_interface.get_new_jobs("test_worker",
                                                                          "test_project",
                                                                          limit=5)))

    def test_create_job_dir(self):
        """
        create mols in database, and then create job dirs
        create_job_from_db called twice to ensure that key1 is made first
        and get's the 0_ index job dir (just for testing's sake)
        """
        self._make_mol_job(SMILES_A)
        self.dirbuilder.create_job_from_db()

        self._make_mol_job(SMILES_B)
        self.dirbuilder.create_job_from_db()

        jobs = [o.key for o in self.bucket.objects.filter(Prefix=INBOX_PATH)]
        job_dirs = list(set([j.split('/')[1] for j in jobs]))
        job_dirs.sort()
        self.assertTrue(job_dirs[0].startswith("KHHUWSREFUHHDL-IUBNYDRBNA-N"))
        self.assertTrue(job_dirs[1].startswith("UHOVQNZJYSORNB-UHFFFAOYNA-N"))
        self.assertTrue(len(job_dirs) == 2)

        for job in Job.objects.all():
            self.assertEquals(job.status, 'claimed')

    def test_create_job_dir_batch(self):
        """
        create mols in database, and then create job dirs
        create_job_from_db called twice to ensure that key1 is made first
        and get's the 0_ index job dir (just for testing's sake)
        """
        self._make_mol_job(SMILES_A, priority=0)
        self._make_mol_job(SMILES_B, priority=1)
        jobs_returned = list(self.dirbuilder.iter_create_jobs_from_db(batch_size=2))

        jobs = [o.key for o in self.bucket.objects.filter(Prefix=INBOX_PATH)]
        sub_job_dirs = list(set([j.split('/')[2] for j in jobs]))

        self.assertTrue("1_UHOVQNZJYSORNB-UHFFFAOYNA-N" in sub_job_dirs)

        self.assertTrue("0_KHHUWSREFUHHDL-IUBNYDRBNA-N" in sub_job_dirs)

        for job in Job.objects.all():
            self.assertEquals(job.status, 'claimed')

    def test_job_script_content(self):
        self._make_mol_job(SMILES_A, 2)
        storage = self.dirbuilder.create_job_from_db()
        job_file_path = os.path.join(storage.job_path, JOB_TEMPLATE_FILENAME)
        script_content = self.bucket.Object(job_file_path).get()['Body'].read().decode('utf-8')

        expect = """Test Input File\n{}""".format(SMILES_A)
        self.assertEquals(script_content, expect)

    def test_job_script_content_with_jobspec_prep(self):
        self.dirbuilder = JobDirBuilder(name="test_worker",
                                        project=PROJECT_NAME,
                                        config_path=CONFIG_PATH,
                                        db_interface=self.db_interface,
                                        job_filename="test_with_extra.inp",
                                        template_filenames=[],
                                        storage_class=S3WriteDir,
                                        storage_kwargs={'bucket_name': DEFAULT_BUCKET_NAME},
                                        jobspec_prep_module_name='jobspec_prep'
                                        )

        self._make_mol_job(SMILES_A, 2)
        storage = self.dirbuilder.create_job_from_db()
        job_file_path = os.path.join(storage.job_path, "test_with_extra.inp")
        script_content = self.bucket.Object(job_file_path).get()['Body'].read().decode('utf-8')

        expect = """Test Input File\n{}\nThis should make it into the template""".format(SMILES_A)
        self.assertEquals(script_content, expect)

    def test_info_file(self):
        mol, job = self._make_mol_job(SMILES_A)
        storage = self.dirbuilder.create_job_from_db()
        info_file_path = os.path.join(storage.job_path, 'job_info.json')
        json_content = self.bucket.Object(info_file_path).get()['Body'].read().decode('utf-8')
        info_dict = json.loads(json_content)
        expect = {'details': None,
                  'inchikey': 'UHOVQNZJYSORNB-UHFFFAOYNA-N',
                  'priority': 0,
                  'smiles': SMILES_A,
                  'project_name': PROJECT_NAME,
                  'worker_name': 'test_worker',
                  'job_key': job.id,
                  'charge': 0}
        self.assertDictEqual(info_dict, expect)


    def test_job_dir_to_db(self):
        self._make_mol_job(SMILES_A)
        storage = self.dirbuilder.create_job_from_db()
        # make complete file

        with open(os.path.join(CONFIG_PATH, "test.xyz")) as fp:
            storage.dump_file("UHOVQNZJYSORNB-UHFFFAOYNA-N_Conf_1.xyz", fp.read())

        self.dirparser = JobDirParser(name="test_worker",
                                      project=PROJECT_NAME,
                                      loader_path=CONFIG_PATH,
                                      loader_module_name="loader",
                                      storage_class=S3ReadDir,
                                      db_interface=self.db_interface,
                                      storage_kwargs={'bucket_name': DEFAULT_BUCKET_NAME})
        self.dirparser.update_db_for_job(storage.job_path)
        geom_num = Geom.objects.count()
        self.assertEquals(geom_num, 1)

        archive_path = os.path.join('archive', *storage.job_path.split('/')[1:])

        self.assertTrue(len(list(self.bucket.objects.filter(Prefix=archive_path))) > 0)



