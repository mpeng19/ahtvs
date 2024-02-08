import json
import os
import unittest
import logging

try:
    import boto3
    from jobs.storage.s3storage import S3ReadDir, S3WriteDir
    BOTO_AVAIL = True
except ImportError:
    BOTO_AVAIL = False

try:
    from aag_python.molecular_storage import molecular_data_models as mdm
    from jobs.dbinterface import mongointerface
    MONGO_AVAIL = True
except ImportError:
    MONGO_AVAIL = False

from jobs.jobdirbuilder import JobDirBuilder
from jobs.jobdirparser import JobDirParser

MY_DIR = os.path.dirname(__file__)
CONFIG_PATH = os.path.join(MY_DIR, "config")

WORKER_NAME = "test_worker"

SMILES_STR = "c1cccc1c"

WORKER_NAME_KEY = "worker_name"
PARENT_KEY_KEY = "parent_key"
PARENT_CLASS_KEY = "parent_class"
PROJECT_NAME = "test_project"
BATCH_KEY = "batch"
DEFAULT_BUCKET_NAME = "calculario.unittest"

JOB_TEMPLATE_FILENAME = "test.inp"
INBOX_PATH = "inbox"

@unittest.skipIf(not (BOTO_AVAIL and MONGO_AVAIL), reason="requires boto3 and aag_python")
class TestMongoS3JobDirHandlers(unittest.TestCase):
    def _db_interface(self):
        return mongointerface.MongoInterface(job_class_name="Molecule")

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

        s3 = boto3.resource('s3')
        self.bucket = s3.Bucket(DEFAULT_BUCKET_NAME)

    def _make_mol(self, key="key1", priority=0):
        status = mdm.Status(current_work_order=WORKER_NAME, current_work_priority=priority)

        meta = mdm.MetaData(key=key,
                            inchi_key=key,
                            data_type=mdm.MetaData.DataTypeValue.MOL_MOLECULE,
                            status=status,
                            smiles=SMILES_STR,
                            project_name=PROJECT_NAME)

        mol = mdm.Molecule(meta_data=meta)
        mol.save()

    def tearDown(self):
        mdm.Molecule.objects.filter(meta_data__project_name=PROJECT_NAME).delete()
        mdm.Calculation.objects.filter(meta_data__project_name=PROJECT_NAME).delete()
        for key in self.bucket.objects.all():
            key.delete()

    def test_find(self):
        self._make_mol()
        self.assertEqual(1, len(self.dirbuilder.db_interface.get_new_jobs("test_worker",
                                                                          "test_project",
                                                                          limit=5)))

    def test_create_job_dir(self):
        """
        create mols in database, and then create job dirs
        create_job_from_db called twice to ensure that key1 is made first
        and get's the 0_ index job dir (just for testing's sake)
        """
        self._make_mol("key1")
        self.dirbuilder.create_job_from_db()

        self._make_mol("key2")
        self.dirbuilder.create_job_from_db()

        jobs = [o.key for o in self.bucket.objects.filter(Prefix=INBOX_PATH)]

        job_dirs = list(set([j.split('/')[1] for j in jobs]))
        job_dirs.sort()
        self.assertTrue(job_dirs[0].startswith("key1"))
        self.assertTrue(job_dirs[1].startswith("key2"))
        self.assertTrue(len(job_dirs) == 2)

        for m in mdm.Molecule.objects.all():
            self.assertEquals(m.meta_data.status.current_status, 'clc_run')

    def test_job_script_content(self):
        self._make_mol("key1")
        storage = self.dirbuilder.create_job_from_db()
        job_file_path = os.path.join(storage.job_path, JOB_TEMPLATE_FILENAME)
        script_content = self.bucket.Object(job_file_path).get()['Body'].read().decode('utf-8')

        expect = """Test Input File\n{}""".format(SMILES_STR)
        self.assertEquals(script_content, expect)

    def test_info_file(self):
        self._make_mol("key1")
        storage = self.dirbuilder.create_job_from_db()
        info_file_path = os.path.join(storage.job_path, 'job_info.json')
        json_content = self.bucket.Object(info_file_path).get()['Body'].read().decode('utf-8')

        info_dict = json.loads(json_content)
        expect = {'inchikey': 'key1',
                  'job_class': 'Molecule',
                  'job_key': 'key1',
                  'parent_class': 'Molecule',
                  'parent_key': 'key1',
                  'priority': 0,
                  'smiles': 'c1cccc1c',
                  'worker_name': 'test_worker',
                  'project_name': PROJECT_NAME}

        self.assertDictEqual(info_dict, expect)

    def test_job_dir_to_db(self):
        self._make_mol("key1")
        storage = self.dirbuilder.create_job_from_db()
        # make complete file

        with open(os.path.join(CONFIG_PATH, "test.xyz")) as fp:
            storage.dump_file("key1_Conf_1.xyz", fp.read())

        self.dirparser = JobDirParser(name="test_worker",
                                      project=PROJECT_NAME,
                                      loader_path=CONFIG_PATH,
                                      loader_module_name="loader",
                                      storage_class=S3ReadDir,
                                      db_interface=self.db_interface,
                                      storage_kwargs={'bucket_name': DEFAULT_BUCKET_NAME})
        self.dirparser.update_db_for_job(storage.job_path)
        calc_num = mdm.Calculation.objects.count()
        self.assertEquals(calc_num, 1)

        archive_path = os.path.join('archive', *storage.job_path.split('/')[1:])

        self.assertTrue(len(list(self.bucket.objects.filter(Prefix=archive_path))) > 0)
