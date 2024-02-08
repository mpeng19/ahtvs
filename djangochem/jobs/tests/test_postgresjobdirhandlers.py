import json
import glob
import shutil


from jobs.jobdirbuilder import JobDirBuilder
from jobs.jobdirparser import JobDirParser
from jobs.dbinterface import postgresinterface
from django.test import TestCase

from jobs.models import *
from pgmols.models import *

MY_DIR = os.path.dirname(__file__)
SCRATCH_PATH = os.path.join(MY_DIR, "scratch")
CONFIG_PATH = os.path.join(MY_DIR, "config")

WORKER_NAME = "test_worker"

SMILES_A = "c1ccccc1"
SMILES_B = "c1cc([Mg])ncc1"
SMILES_C = "c1nc([Mg])cc([Mg])c1"
SMILES_D = "[Mg][Cl]"


WORKER_NAME_KEY = "worker_name"
PARENT_CLASS_KEY = "parent_class"
PROJECT_NAME = "test_project"
BATCH_KEY = "batch"

JOB_TEMPLATE_FILENAME = "test.inp"


class TestPostgresJobDirHandlers(TestCase):

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

        self.db_interface = postgresinterface.PostgresInterface()
        self.dirbuilder = JobDirBuilder(name="test_worker",
                                        project=PROJECT_NAME,
                                        config_path=CONFIG_PATH,
                                        db_interface=self.db_interface,
                                        job_filename=JOB_TEMPLATE_FILENAME,
                                        template_filenames=[],
                                        storage_kwargs={
                                            "inbox_job_dir": self.inbox_path},
                                        batch_filename='batch.sh'
                                        )

        self.config = JobConfig(name="test_worker", parent_class_name="Mol")
        self.config.save()
        self.test_group = Group(name=PROJECT_NAME)
        self.test_group.save()

        self.dirbuilder._make_batch_key_name = self._mock_time_stamper

    def _safe_move(self, src, tgt):
        """ For some reason, shutil.move doesn't work inside a docker environment.  Instead, this method does a copyree followed by rmtree.
            my guess is that it's something to do with preserving permissions.  To be clear, the call to move works, but then you aren't allowed
            to write any new files into it.  It's something edgy."""
        shutil.copytree(src, tgt)
        shutil.rmtree(src)

    def _mock_time_stamper(self):
        return "batch"

    def _make_mol_job(self, smiles=SMILES_A, priority=0):
        mol = Mol(smiles=smiles)
        mol.save()
        job = Job(parent=mol,
                  config=self.config,
                  group=self.test_group,
                  priority=priority)
        job.save()
        return mol, job

    def tearDown(self):
        shutil.rmtree(SCRATCH_PATH)

    def test_find(self):
        self._make_mol_job()
        self.assertEqual(1, len(self.dirbuilder.db_interface.get_new_jobs("test_worker",
                                                                          "test_project",
                                                                          limit=5)))

    # def test_create_job_dir(self):
    #     """
    #     create mols in database, and then create job dirs
    #     create_job_from_db called twice to ensure that key1 is made first
    #     and get's the 0_ index job dir (just for testing's sake)
    #     """
    #     self._make_mol_job(SMILES_A)
    #     self.dirbuilder.create_job_from_db()
    #
    #     self._make_mol_job(SMILES_B)
    #     self.dirbuilder.create_job_from_db()
    #
    #     jobs = os.listdir(self.inbox_path)
    #     self.assertIn("0000_test_worker_UHOVQNZJYSORNB-UHFFFAOYNA-N", jobs)
    #     self.assertIn("0001_test_worker_KHHUWSREFUHHDL-IUBNYDRBNA-N", jobs)
    #
    #     for job in Job.objects.all():
    #         self.assertEquals(job.status, 'claimed')

    # def test_create_job_dir_batch(self):
    #     """
    #     create mols in database, and then create job dirs
    #     create_job_from_db called twice to ensure that key1 is made first
    #     and get's the 0_ index job dir (just for testing's sake)
    #     """
    #     self._make_mol_job(SMILES_A, priority=0)
    #     self._make_mol_job(SMILES_B, priority=1)
    #     jobs = list(self.dirbuilder.build_job_dirs(batch_size=2))
    #
    #     expect = os.path.join(self.inbox_path, "0000_test_worker_batch")
    #     self.assertTrue(os.path.isdir(expect + "/1_UHOVQNZJYSORNB-UHFFFAOYNA-N"),
    #                     msg="Found: " + "\n".join(os.listdir(expect)))
    #     jobs = os.listdir(self.inbox_path)
    #     self.assertTrue(os.path.isdir(
    #         expect + "/0_KHHUWSREFUHHDL-IUBNYDRBNA-N"))
    #
    #     for job in Job.objects.all():
    #         self.assertEquals(job.status, 'claimed')

    # def test_create_job_dir_priority(self):
    #     """
    #     create mols in database, and then create job dirs
    #     create_job_from_db called twice to ensure that key1 is made first
    #     and get's the 0_ index job dir (just for testing's sake)
    #     """
    #     self._make_mol_job(SMILES_A, priority=5)
    #     self.dirbuilder.create_job_from_db()
    #
    #     self._make_mol_job(SMILES_B, priority=5)
    #     self.dirbuilder.create_job_from_db()
    #     jobs = os.listdir(self.inbox_path)
    #     self.assertIn("5000_test_worker_UHOVQNZJYSORNB-UHFFFAOYNA-N", jobs)
    #     self.assertIn("5001_test_worker_KHHUWSREFUHHDL-IUBNYDRBNA-N", jobs)

    def test_job_script_content(self):
        self._make_mol_job(SMILES_A, 2)
        self.dirbuilder.create_job_from_db()
        file_path = os.path.join(self.inbox_path, "2000_test_worker_UHOVQNZJYSORNB-UHFFFAOYNA-N_*", JOB_TEMPLATE_FILENAME)
        file_path = glob.glob(file_path)[0]
        with open(file_path) as fp:
            script_content = fp.read()
        expect = """Test Input File\n{}""".format(SMILES_A)
        self.assertEquals(script_content, expect)

    def test_job_script_content_with_jobspec_prep(self):
        self.dirbuilder = JobDirBuilder(name="test_worker",
                                        project=PROJECT_NAME,
                                        config_path=CONFIG_PATH,
                                        db_interface=self.db_interface,
                                        job_filename="test_with_extra.inp",
                                        template_filenames=[],
                                        storage_kwargs={
                                            "inbox_job_dir": self.inbox_path},
                                        jobspec_prep_module_name='jobspec_prep'
                                        )

        self._make_mol_job(SMILES_A, 2)
        self.dirbuilder.create_job_from_db()
        file_path = os.path.join(self.inbox_path, "2000_test_worker_UHOVQNZJYSORNB-UHFFFAOYNA-N_*", "test_with_extra.inp")
        file_path = glob.glob(file_path)[0]
        with open(file_path) as fp:
            script_content = fp.read()
        expect = """Test Input File\n{}\nThis should make it into the template""".format(SMILES_A)
        self.assertEquals(script_content, expect)

    def test_info_file(self):
        mol, job = self._make_mol_job(SMILES_A)
        self.dirbuilder.create_job_from_db()
        file_path = os.path.join(self.inbox_path, "0000_test_worker_UHOVQNZJYSORNB-UHFFFAOYNA-N_*", "job_info.json")
        file_path = glob.glob(file_path)[0]
        with open(file_path) as fp:
            info_dict = json.load(fp)
        expect = {'details': None,
                  'inchikey': 'UHOVQNZJYSORNB-UHFFFAOYNA-N',
                  'priority': 0,
                  'smiles': SMILES_A,
                  'project_name': PROJECT_NAME,
                  'parent_class': 'Mol',
                  'worker_name': 'test_worker',
                  'job_key': job.id,
                  'uuid': str(job.uuid),
                  'charge': 0}
        self.assertDictEqual(info_dict, expect)

    def test_parse(self):
        mol, job = self._make_mol_job(SMILES_A)
        self.dirbuilder.create_job_from_db()
        job_dir = "0000_test_worker_UHOVQNZJYSORNB-UHFFFAOYNA-N_*/"
        # make complete file
        inbox_jobdir = os.path.join(self.inbox_path, job_dir)
        inbox_jobdir = glob.glob(inbox_jobdir)[0]
        job_dir = os.path.normpath(inbox_jobdir).split(os.path.sep)[-1]
        pretend_completed_job_dir = os.path.join(self.completed_path, job_dir)
        self._safe_move(inbox_jobdir,
                    pretend_completed_job_dir)
        shutil.copy(os.path.join(CONFIG_PATH, "test.xyz"),
                    os.path.join(pretend_completed_job_dir, "UHOVQNZJYSORNB-UHFFFAOYNA-N_Conf_1.xyz"))

        self.dirparser = JobDirParser(name="test_worker",
                                      project=PROJECT_NAME,
                                      loader_path=CONFIG_PATH,
                                      loader_module_name="loader",
                                      db_interface=self.db_interface,
                                      storage_kwargs={'error_path': self.error_path,
                                                      'archive_path': self.archive_path})
        self.dirparser.update_db_for_job(pretend_completed_job_dir)
        geom_num = Geom.objects.count()
        self.assertEquals(geom_num, 1)

        def test_hour(self):
            return "5h-4-14-14"
        self.dirparser._hour_dir_name = test_hour
        os.path.isdir(os.path.join(self.archive_path, "5h-4-14-14", job_dir))

    # def test_parse_batch(self):
    #     """
    #     create mols in database, and then create job dirs
    #     create_job_from_db called twice to ensure that key1 is made first
    #     and get's the 0_ index job dir (just for testing's sake)
    #     """
    #     self._make_mol_job(SMILES_A, priority=0)
    #     self._make_mol_job(SMILES_B, priority=1)
    #     list(self.dirbuilder.build_job_dirs(batch_size=2))
    #
    #     job_dir = "0000_test_worker_batch"
    #     pretend_completed_job_dir = os.path.join(self.completed_path, job_dir)
    #     self._safe_move(os.path.join(self.inbox_path, job_dir), pretend_completed_job_dir)
    #     targetpath0 = os.path.join(pretend_completed_job_dir, "0_KHHUWSREFUHHDL-IUBNYDRBNA-N", "KHHUWSREFUHHDL-IUBNYDRBNA-N_Conf_1.xyz")
    #     shutil.copy(os.path.join(CONFIG_PATH, "test.xyz"), targetpath0)
    #     targetpath1 = os.path.join(pretend_completed_job_dir, "1_UHOVQNZJYSORNB-UHFFFAOYNA-N", "UHOVQNZJYSORNB-UHFFFAOYNA-N_Conf_1.xyz")
    #     shutil.copy(os.path.join(CONFIG_PATH, "test.xyz"), targetpath1)
    #
    #     self.dirparser = JobDirParser(name="test_worker",
    #                                   project=PROJECT_NAME,
    #                                   loader_path=CONFIG_PATH,
    #                                   loader_module_name="loader",
    #                                   db_interface=self.db_interface,
    #                                   storage_kwargs={'error_path': self.error_path,
    #                                                   'archive_path': self.archive_path})
    #
    #     self.dirparser.update_db_for_completed_jobs(self.completed_path)
    #     geom_num = Geom.objects.count()
    #     self.assertEquals(geom_num, 2)
