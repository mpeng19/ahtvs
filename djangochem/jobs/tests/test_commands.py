import json
import os
import shutil
import glob

from django.core.management import call_command
from django.test import TestCase
from django.utils.six import StringIO

from jobs.models import *
from pgmols.models import *

from rdkit import RDLogger
RDLogger.logger().setLevel(RDLogger.ERROR)


MY_DIR = os.path.dirname(__file__)
FIXTURES_DIR = os.path.join(MY_DIR, 'fixtures')
SCRATCH_PATH = os.path.join(MY_DIR, 'scratch')
SMILES_A = "c1ccccc1"
SMILES_B = "c1cc([Mg])ncc1"
SMILES_C = "c1nc([Mg])cc([Mg])c1"


PROJECT_NAME = "test_project"
CONFGEN_NAME = "pretend_confgen"
CONFIG_NAME = "test_worker"

CONFIG_DIR = os.path.join(MY_DIR, 'config')
CONFIG_PATH = os.path.join(CONFIG_DIR, 'test_config.json')

# JOBS_QUERY_OUTPUT_EXPECT = [
#     "  job_key    charge    details  inchikey                       priority  project_name    smiles          worker_name",
#     "---------  --------  ---------  ---------------------------  ----------  --------------  --------------  -------------",
#     "        1         0             UHOVQNZJYSORNB-UHFFFAOYNA-N           0  test_project    c1ccccc1        test_worker",
#     "        2         0             KHHUWSREFUHHDL-IUBNYDRBNA-N           1  test_project    c1cc([Mg])ncc1  test_worker",
#     ""
# ]


class CommandTest(TestCase):
    # reset_sequences = True

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

        self.config = JobConfig(name=CONFIG_NAME,
                                configpath=CONFIG_PATH,
                                parent_class_name="Mol")
        self.config.save()
        self.test_group = Group(name=PROJECT_NAME)
        self.test_group.save()

    def _safe_move(self, src, tgt):
        """ For some reason, shutil.move doesn't work inside a docker environment.  Instead, this method does a copyree followed by rmtree.
            my guess is that it's something to do with preserving permissions.  To be clear, the call to move works, but then you aren't allowed
            to write any new files into it.  It's something edgy."""
        shutil.copytree(src, tgt)
        shutil.rmtree(src)
    
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

    # def test_build_dir(self):
    #     out = StringIO()
    #
    #     self._make_mol_job(SMILES_A, priority=0)
    #     self._make_mol_job(SMILES_B, priority=1)
    #
    #     call_command('buildjobs',
    #                  PROJECT_NAME,
    #                  self.inbox_path,
    #                  '-c'+CONFIG_NAME,
    #                  stdout=out)
    #
    #     expected = ["0000_test_worker_UHOVQNZJYSORNB-UHFFFAOYNA-N_*",
    #                 "1000_test_worker_KHHUWSREFUHHDL-IUBNYDRBNA-N_*"]
    #
    #     jobs = os.listdir(self.inbox_path)
    #     for exp in expected:
    #         self.assertIn(exp, jobs)
    #
    #     for job in Job.objects.all():
    #         self.assertEquals(job.status, 'claimed')
    #
    #     job_list = out.getvalue().split("\n")
    #     for exp in expected:
    #         self.assertIn("created: {}".format(
    #             os.path.join(self.inbox_path, exp)), job_list)

    def test_build_dir_with_batch(self):
        out = StringIO()

        self._make_mol_job(SMILES_A, priority=0)
        self._make_mol_job(SMILES_B, priority=1)

        call_command('buildjobs',
                     PROJECT_NAME,
                     self.inbox_path,
                     '-c'+CONFIG_NAME,
                     '-b2',
                     stdout=out)

        expected = "0000_test_worker_batch"

        job = os.listdir(self.inbox_path)[0]
        self.assertIn(expected, job)

        reported_job = out.getvalue().split("\n")[0]
        expected_report = "created: {}".format(os.path.join(self.inbox_path, expected))
        self.assertTrue(reported_job.startswith(expected_report))

    def test_query_only(self):
        out = StringIO()

        self._make_mol_job(SMILES_A, priority=0)
        self._make_mol_job(SMILES_B, priority=1)

        call_command('buildjobs',
                     PROJECT_NAME,
                     self.inbox_path,
                     '-c'+CONFIG_NAME,
                     '-q',
                     stdout=out)

        self.assertEquals(os.listdir(self.inbox_path), [])
        for job in Job.objects.all():
            self.assertEquals(job.status, '')
        # query_report = out.getvalue().split("\n")
        # self.assertEquals(JOBS_QUERY_OUTPUT_EXPECT, query_report)

    def test_parse(self):
        mol, job = self._make_mol_job(SMILES_A)
        out = StringIO()
        call_command('buildjobs',
                     PROJECT_NAME,
                     self.inbox_path,
                     '-c'+CONFIG_NAME,
                     stdout=out)
        inchikey = "UHOVQNZJYSORNB-UHFFFAOYNA-N"
        #match time stampt
        job_dir = "0000_test_worker_{}_*/".format(inchikey)
        # make complete file
        inbox_jobdir = os.path.join(self.inbox_path, job_dir)
        inbox_jobdir = glob.glob(inbox_jobdir)[0]
        job_dir = "0000_test_worker_{}_{}".format(inchikey,inbox_jobdir.split("_")[-1])
        pretend_completed_job_dir = os.path.join(self.completed_path, job_dir)
        self._safe_move(inbox_jobdir, pretend_completed_job_dir)
        target_path = os.path.join(pretend_completed_job_dir, "{}_Conf_1.xyz".format(inchikey))
        src = os.path.join(CONFIG_DIR, "test.xyz")
        #import pdb; pdb.set_trace()
        with open(src) as infile:
            with open(target_path, 'w') as outfile:
                outfile.write(infile.read())

        call_command('parsejobs',
                     PROJECT_NAME,
                     self.completed_path,
                     '-c'+CONFIG_NAME,
                     '-r' + SCRATCH_PATH,
                     stdout=out)
        geom_num = Geom.objects.count()
        self.assertEquals(geom_num, 1)

    def test_parse_jobdir(self):
        mol, job = self._make_mol_job(SMILES_A)
        out = StringIO()
        call_command('buildjobs',
                     PROJECT_NAME,
                     self.inbox_path,
                     '-c'+CONFIG_NAME,
                     stdout=out)
        inchikey = "UHOVQNZJYSORNB-UHFFFAOYNA-N"
        #match time stampt
        job_dir = "0000_test_worker_{}_*/".format(inchikey)
        # make complete file
        inbox_jobdir = os.path.join(self.inbox_path, job_dir)
        inbox_jobdir = glob.glob(inbox_jobdir)[0]
        job_dir = "0000_test_worker_{}_{}".format(inchikey,inbox_jobdir.split("_")[-1])
        pretend_completed_job_dir = os.path.join(self.completed_path, job_dir)
        self._safe_move(os.path.join(self.inbox_path, job_dir),
                    pretend_completed_job_dir)
        shutil.copyfile(os.path.join(CONFIG_DIR, "test.xyz"),
                        os.path.join(pretend_completed_job_dir, "{}_Conf_1.xyz".format(inchikey)))

        call_command('parsejobs',
                     PROJECT_NAME,
                     os.path.join(self.completed_path, pretend_completed_job_dir),
                     '-c'+CONFIG_NAME,
                     '-r' + SCRATCH_PATH,
                     stdout=out)
        geom_num = Geom.objects.count()
        self.assertEquals(geom_num, 1)

    def test_parse_no_config_specified(self):
        mol, job = self._make_mol_job(SMILES_A)
        out = StringIO()
        call_command('buildjobs',
                     PROJECT_NAME,
                     self.inbox_path,
                     stdout=out)
        inchikey = "UHOVQNZJYSORNB-UHFFFAOYNA-N"
        job_dir = "0000_test_worker_{}_*/".format(inchikey)
        # make complete file
        inbox_jobdir = os.path.join(self.inbox_path, job_dir)
        inbox_jobdir = glob.glob(inbox_jobdir)[0]
        job_dir = "0000_test_worker_{}_{}".format(inchikey,inbox_jobdir.split("_")[-1])
        # make complete file
        pretend_completed_job_dir = os.path.join(self.completed_path, job_dir)
        self._safe_move(os.path.join(self.inbox_path, job_dir),
                    pretend_completed_job_dir)
        shutil.copyfile(os.path.join(CONFIG_DIR, "test.xyz"),
                        os.path.join(pretend_completed_job_dir, "{}_Conf_1.xyz".format(inchikey)))

        call_command('parsejobs',
                     PROJECT_NAME,
                     self.completed_path,
                     '-r' + SCRATCH_PATH,
                     stdout=out)
        geom_num = Geom.objects.count()
        self.assertEquals(geom_num, 1)

    def test_parse_with_limit(self):
        mol, job = self._make_mol_job(SMILES_A, priority=0)
        mol, job2 = self._make_mol_job(SMILES_B, priority=1)
        out = StringIO()
        call_command('buildjobs',
                     PROJECT_NAME,
                     self.inbox_path,
                     '-c'+CONFIG_NAME,
                     stdout=out)
        inchikey = "UHOVQNZJYSORNB-UHFFFAOYNA-N"
        job_dir = "0000_test_worker_{}_*/".format(inchikey)
        # make complete file
        inbox_jobdir = os.path.join(self.inbox_path, job_dir)
        inbox_jobdir = glob.glob(inbox_jobdir)[0]
        job_dir = "0000_test_worker_{}_{}".format(inchikey,inbox_jobdir.split("_")[-1])
        # make complete file
        pretend_completed_job_dir = os.path.join(self.completed_path, job_dir)
        self._safe_move(os.path.join(self.inbox_path, job_dir),
                    pretend_completed_job_dir)
        shutil.copyfile(os.path.join(CONFIG_DIR, "test.xyz"),
                        os.path.join(pretend_completed_job_dir, "{}_Conf_1.xyz".format(inchikey)))
        out = StringIO()
        call_command('parsejobs',
                     PROJECT_NAME,
                     self.completed_path,
                     '-c'+CONFIG_NAME,
                     '-l1',
                     '-r' + SCRATCH_PATH,
                     stdout=out)
        geom_num = Geom.objects.count()
        self.assertEquals(geom_num, 1)


class RequestTest(TestCase):

    def setUp(self):
        JobConfig.objects.create(name=CONFGEN_NAME,
                                 parent_class_name="Mol")
        # JobConfig.objects.create(name=DFT_NAME)

        test_group = Group.objects.create(name=PROJECT_NAME)

        Mol.objects.create(smiles=SMILES_A, group=test_group)
        Mol.objects.create(smiles=SMILES_B, group=test_group)
        Mol.objects.create(smiles=SMILES_C, group=test_group)

        self.mmff_method = Method.objects.create(name="mmff")

    def test_request_job(self):
        out = StringIO()
        call_command('requestjobs',
                     PROJECT_NAME,
                     CONFGEN_NAME,
                     stdout=out)
        self.assertEquals(Job.objects.count(), 3)

    def test_request_job_with_limit(self):
        out = StringIO()
        call_command('requestjobs',
                     PROJECT_NAME,
                     CONFGEN_NAME,
                     '--limit=1',
                     stdout=out)
        self.assertEquals(Job.objects.count(), 1)

    def test_request_job_with_details(self):
        out = StringIO()
        call_command('requestjobs',
                     PROJECT_NAME,
                     CONFGEN_NAME,
                     '--details={"this": "super"}',
                     stdout=out)
        self.assertEquals(Job.objects.count(), 3)
        self.assertEquals(Job.objects.all()[0].details["this"], "super")


class ScanConfigTest(TestCase):

    def _assert_lines_match(self, out_lines, expect_lines):
        """
        ignoring whitespace, check for equality of every line
        """
        expect_lines = set([s.strip() for s in expect_lines if s])
        out_lines = set([s.strip() for s in out_lines if s])

        self.assertEqual(out_lines, expect_lines)

    def test_scan(self):
        out = StringIO()
        call_command('scanconfigs',
                     CONFIG_DIR,
                     '--update',
                     stdout=out)
        self.assertEquals(JobConfig.objects.count(), 1)

        expect_lines = ["ERROR -- Config: " + CONFIG_DIR + "/bad_config.json:",
                        "Missing Keys",
                        "- loader_module",
                        "Adding: test_worker from path " + CONFIG_DIR + "/test_config.json",
                        "Config Count: 2. Valid Count: 1",
                        "Added 1 new JobConfigs"]
        out_lines = out.getvalue().split('\n')
        self._assert_lines_match(out_lines, expect_lines)

    def test_scan_dryrun(self):
        out = StringIO()
        call_command('scanconfigs',
                     CONFIG_DIR,
                     stdout=out)
        self.assertEquals(JobConfig.objects.count(), 0)
        expect_lines = ["ERROR -- Config: " + CONFIG_DIR + "/bad_config.json:",
                        "Missing Keys",
                        "- loader_module",
                        "Found (dryrun): test_worker at path " + CONFIG_DIR + "/test_config.json",
                        "Config Count: 2. Valid Count: 1"]
        out_lines = out.getvalue().split('\n')
        self._assert_lines_match(out_lines, expect_lines)

    def test_scan_extra(self):
        out = StringIO()
        EXTRA_CONFIG = os.path.join(MY_DIR, 'config_w_extra', 'test_config.json')
        call_command('scanconfigs',
                     EXTRA_CONFIG,
                     '--update',
                     stdout=out)
        self.assertEquals(JobConfig.objects.count(), 1)
        expect_lines = ["Adding: test_worker from path " + EXTRA_CONFIG,
                        "Config Count: 1. Valid Count: 1",
                        "Added 1 new JobConfigs"]
        out_lines = out.getvalue().split('\n')
        self._assert_lines_match(expect_lines, out_lines)

    def test_missing_extra(self):
        out = StringIO()
        EXTRA_CONFIG = os.path.join(MY_DIR, 'config_w_extra', 'test_noprep_config.json')
        call_command('scanconfigs',
                     EXTRA_CONFIG,
                     '--update',
                     stdout=out)
        self.assertEquals(JobConfig.objects.count(), 0)
        expect_lines = ["ERROR -- Config: " + EXTRA_CONFIG + ":",
                        "Template File: " + os.path.dirname(EXTRA_CONFIG) + "/test_with_extra.inp",
                        "Template variable: dict object has no attribute extra_arg (jobspec.extra_arg on line 4)",
                        "Config Count: 1. Valid Count: 0",
                        "Added 0 new JobConfigs"]
        out_lines = out.getvalue().split('\n')
        self._assert_lines_match(expect_lines, out_lines)



