import json
import os
import unittest
import shutil

try:
    from aag_python.molecular_storage import molecular_data_models as mdm
    from jobs.dbinterface import mongointerface
    MONGO_AVAIL = True
except ImportError:
    MONGO_AVAIL = False

from jobs.jobdirbuilder import JobDirBuilder
from jobs.jobdirparser import JobDirParser

MY_DIR = os.path.dirname(__file__)
SCRATCH_PATH = os.path.join(MY_DIR, "scratch")
CONFIG_PATH = os.path.join(MY_DIR, "config")
WORKER_NAME = "test_worker"

SMILES_STR = "c1cccc1c"

WORKER_NAME_KEY = "worker_name"
PARENT_KEY_KEY = "parent_key"
PARENT_CLASS_KEY = "parent_class"
PROJECT_NAME = "test_project"
BATCH_KEY = "batch"

JOB_TEMPLATE_FILENAME = "test.inp"

@unittest.skipIf(not MONGO_AVAIL, reason="requires mongo via aag_python")
class TestMongoJobDirHandlers(unittest.TestCase):
    def setUp(self):
        self.inbox_path = os.path.join(SCRATCH_PATH, "inbox")
        self.completed_path = os.path.join(SCRATCH_PATH, "completed")
        self.archive_path = os.path.join(SCRATCH_PATH, "archive")
        self.error_path = os.path.join(SCRATCH_PATH, "errors")

        if not os.path.isdir(self.inbox_path):
            os.makedirs(self.inbox_path)
        if not os.path.isdir(self.completed_path):
            os.makedirs(self.completed_path)
        if not os.path.isdir(self.archive_path):
            os.makedirs(self.archive_path)
        if not os.path.isdir(self.error_path):
            os.makedirs(self.error_path)

        self.db_interface = mongointerface.MongoInterface(job_class_name="Molecule")
        self.dirbuilder = JobDirBuilder(name="test_worker",
                                        project=PROJECT_NAME,
                                        config_path=CONFIG_PATH,
                                        db_interface=self.db_interface,
                                        job_filename=JOB_TEMPLATE_FILENAME,
                                        template_filenames=[],
                                        storage_kwargs={"inbox_job_dir": self.inbox_path}
                                        )

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
        shutil.rmtree(SCRATCH_PATH)

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

        jobs = os.listdir(self.inbox_path)
        self.assertIn("0000_test_worker_key1", jobs)
        self.assertIn("0001_test_worker_key2", jobs)

        for m in mdm.Molecule.objects.all():
            self.assertEquals(m.meta_data.status.current_status, 'clc_run')

    def test_create_job_dir_priority(self):
        """
        create mols in database, and then create job dirs
        create_job_from_db called twice to ensure that key1 is made first
        and get's the 0_ index job dir (just for testing's sake)
        """
        self._make_mol("key1", 5)
        self.dirbuilder.create_job_from_db()

        self._make_mol("key2", 5)
        self.dirbuilder.create_job_from_db()

        jobs = os.listdir(self.inbox_path)
        self.assertIn("5000_test_worker_key1", jobs)
        self.assertIn("5001_test_worker_key2", jobs)

    def test_job_script_content(self):
        self._make_mol("key1", 2)
        self.dirbuilder.create_job_from_db()
        with open(os.path.join(self.inbox_path, "2000_test_worker_key1", JOB_TEMPLATE_FILENAME)) as fp:
            script_content = fp.read()
        expect = """Test Input File\n{}""".format(SMILES_STR)
        self.assertEquals(script_content, expect)

    def test_info_file(self):
        self._make_mol("key1")
        self.dirbuilder.create_job_from_db()
        with open(os.path.join(self.inbox_path, "0000_test_worker_key1", "job_info.json")) as fp:
            info_dict = json.load(fp)
        expect = {'parent_class': 'Molecule',
                  'inchikey': 'key1',
                  'parent_key': 'key1',
                  'priority': 0,
                  'smiles': 'c1cccc1c',
                  'worker_name': 'test_worker',
                  'job_class': 'Molecule',
                  'job_key': 'key1',
                  'project_name': PROJECT_NAME}
        self.assertDictEqual(info_dict, expect)

    def test_job_dir_to_db(self):
        self._make_mol("key1")
        self.dirbuilder.create_job_from_db()
        job_dir = "0000_test_worker_key1"
        # make complete file
        pretend_completed_job_dir = os.path.join(self.completed_path, job_dir)
        shutil.move(os.path.join(self.inbox_path, job_dir),
                    pretend_completed_job_dir)
        shutil.copy(os.path.join(CONFIG_PATH, "test.xyz"),
                    os.path.join(pretend_completed_job_dir, "key1_Conf_1.xyz"))

        self.dirparser = JobDirParser(name="test_worker",
                                      project=PROJECT_NAME,
                                      loader_path=CONFIG_PATH,
                                      loader_module_name="loader",
                                      db_interface=self.db_interface,
                                      storage_kwargs={'error_path': self.error_path,
                                                      'archive_path': self.archive_path})
        self.dirparser.update_db_for_job(pretend_completed_job_dir)
        calc_num = mdm.Calculation.objects.count()
        self.assertEquals(calc_num, 1)

        def test_hour(self):
            return "5h-4-14-14"
        self.dirparser._hour_dir_name = test_hour
        os.path.isdir(os.path.join(self.archive_path, "5h-4-14-14", job_dir))

