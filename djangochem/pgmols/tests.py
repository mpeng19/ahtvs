from django.test import TestCase
from rdkit.Chem import AllChem as Chem

import logging
from .models import Mol

from rdkit import RDLogger
RDLogger.logger().setLevel(RDLogger.CRITICAL)

BENZENE_SMILES = 'c1ccccc1'
CARBAZOLE_SMILES = 'N1C2=CC=CC=C2C2=CC=CC=C12'
CROWM_ETHER_SMILES = 'C1COCCOCCOCCO1'
VITAMIN_A_SMILES = 'C\C(=C/CO)\C=C\C=C(/C)\C=C\C1=C(C)CCCC1(C)C'

INCHI_OPTION_STRING = " -RecMet  -FixedH "


def smiles_to_inchikey(smiles):
    return Chem.InchiToInchiKey(Chem.MolToInchi(Chem.MolFromSmiles(smiles), options=str(INCHI_OPTION_STRING)))


class TestMol(TestCase):
    # def test_init_with_smiles(self):
    #     mol = Mol(smiles=BENZENE_SMILES)
    #     mol.save()
    #     self.assertEquals(BENZENE_SMILES, Chem.MolToSmiles(mol.rdmol))

    # def test_create_with_smiles(self):
    #     mol = Mol.objects.create(smiles=BENZENE_SMILES)
    #     self.assertEquals(BENZENE_SMILES, Chem.MolToSmiles(mol.rdmol))

    # def test_init_with_smiles2(self):
    #     mol = Mol(rdmol=BENZENE_SMILES)
    #     mol.save()
    #     self.assertEquals(BENZENE_SMILES, Chem.MolToSmiles(mol.rdmol))

    # def test_init_with_mol(self):
    #     mol = Mol(rdmol=Chem.MolFromSmiles(BENZENE_SMILES))
    #     mol.save()
    #     self.assertEquals(BENZENE_SMILES, Chem.MolToSmiles(mol.rdmol))

    def test_inchi(self):
        mol = Mol(smiles=BENZENE_SMILES)
        mol.save()
        self.assertEquals(smiles_to_inchikey(BENZENE_SMILES), mol.inchikey)

    # def test_substructure_search(self):
    #     for sm in [CARBAZOLE_SMILES, VITAMIN_A_SMILES]:
    #         Mol.objects.create(smiles=sm)
    #     query = Mol.objects.filter(rdmol__hassubstruct=BENZENE_SMILES)
    #     self.assertEquals(query.count(), 1)
    #     self.assertEquals(query[0].smiles, CARBAZOLE_SMILES)

    def test_heritage_add_parent(self):
        b = Mol.objects.create(smiles=BENZENE_SMILES)
        c = Mol.objects.create(smiles=CARBAZOLE_SMILES)
        c.parents.add(b)
        self.assertIn(b, c.parents.all())
        self.assertIn(c, b.children.all())

    def test_heritage_add_parentandchild(self):
        a = Mol.objects.create(smiles=VITAMIN_A_SMILES)
        b = Mol.objects.create(smiles=BENZENE_SMILES)
        c = Mol.objects.create(smiles=CARBAZOLE_SMILES)
        b.parents.add(a)
        c.parents.add(b)
        self.assertIn(b, c.parents.all())
        self.assertIn(c, b.children.all())

    def test_heritage_add_child(self):
        b = Mol.objects.create(smiles=BENZENE_SMILES)
        c = Mol.objects.create(smiles=CARBAZOLE_SMILES)
        b.children.add(c)
        self.assertIn(b, c.parents.all())
        self.assertIn(c, b.children.all())
