import os

from django.contrib.auth.models import Group
from django.test import TestCase

from jobs import jobrequester
from jobs.models import *
from pgmols.models import *


SMILES_A = "c1ccccc1"
SMILES_B = "c1cc([Mg])ncc1"
SMILES_C = "c1nc([Mg])cc([Mg])c1"
SMILES_D = "[Mg][Cl]"


PROJECT_NAME = "test_project"
CONFGEN_NAME = "pretend_confgen"
CONFGEN_NAME_WITH_DEFAULTS = "pretend_confgen_w_default"
DFT_NAME = "pretend_dft"

MY_DIR = os.path.dirname(__file__)
CONFIG_PATH = os.path.join(MY_DIR, "config_w_extra", 'test_config.json')

# TODO: Make concurrency tests

class TestPostgresJobDirHandlers(TestCase):

    def setUp(self):
        JobConfig.objects.create(name=CONFGEN_NAME_WITH_DEFAULTS,
                                 configpath=CONFIG_PATH,
                                 parent_class_name="Mol")
        JobConfig.objects.create(name=CONFGEN_NAME,
                                 parent_class_name="Mol")
        JobConfig.objects.create(name=DFT_NAME,
                                 parent_class_name="Geom")

        test_group = Group.objects.create(name=PROJECT_NAME)

        Mol.objects.create(smiles=SMILES_A, group=test_group)
        Mol.objects.create(smiles=SMILES_B, group=test_group)
        Mol.objects.create(smiles=SMILES_C, group=test_group)

        self.mmff_method = Method.objects.create(name="mmff")

    def _make_geom(self, parentjob):
        return Geom.objects.create(mol=parentjob.parent,
                                   method=self.mmff_method,
                                   parentjob=parentjob,
                                   xyz=[[0, 0, 0, 0]])

    def test_mol_request(self):
        jobrequester.request_jobs(project=PROJECT_NAME,
                            requested_config_name=CONFGEN_NAME,
                            parent_config_name=None,
                            limit=None,
                            workbatch=None,
                            details={},
                            avoid_dupes=True)
        self.assertEquals(Job.objects.count(), 3)

    def test_mol_request_with_defaults(self):
        jobrequester.request_jobs(project=PROJECT_NAME,
                            requested_config_name=CONFGEN_NAME_WITH_DEFAULTS,
                            parent_config_name=None,
                            limit=None,
                            workbatch=None,
                            details={},
                            avoid_dupes=True)
        job = Job.objects.all()[0]
        self.assertEquals(job.details['pants'], "fancy")

    def test_mol_request_with_limit(self):
        jobrequester.request_jobs(project=PROJECT_NAME,
                            requested_config_name=CONFGEN_NAME,
                            parent_config_name=None,
                            limit=1,
                            workbatch=None,
                            details={},
                            avoid_dupes=True)
        self.assertEquals(Job.objects.count(), 1)

    def test_geom_request(self):
        jobrequester.request_jobs(project=PROJECT_NAME,
                            requested_config_name=CONFGEN_NAME,
                            parent_config_name=None,
                            limit=None,
                            workbatch=None,
                            details={},
                            avoid_dupes=True)
        for job in Job.objects.all():
            self._make_geom(job)
        jobrequester.request_jobs(project=PROJECT_NAME,
                            requested_config_name=DFT_NAME,
                            parent_config_name=CONFGEN_NAME,
                            workbatch=None,
                            details={},
                            avoid_dupes=True)
        self.assertEquals(Job.objects.count(), 6)
