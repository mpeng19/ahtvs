from pgmols.models import Mol, Group
from rdkit.Chem import AllChem as Chem

def put_block(block, project, tags=[], details={}, return_mol=False):
    '''
    A method for loading a block into a database
    '''
    return put_mol(
        rdkit_mol=block.mol,
        project=project,
        tags=tags,
        details={**{'generation': block.generation}, **details},
        parent_inchikeys=[parent.inchikey for parent in block.parents],
        calc_charge=True,
        return_mol=return_mol,
        inchikey=block.inchikey,
    )

def put_mol(rdkit_mol=None, project=None, tags=[], details={},
            parent_inchikeys=[], calc_charge=True, return_mol=False,
            inchikey=None):
    group, group_created = Group.objects.get_or_create(name=project)
    if not inchikey:
        Chem.SanitizeMol(rdkit_mol)
        inchi_option_string = " -RecMet  -FixedH "
        inchi_string = Chem.MolToInchi(rdkit_mol, options=inchi_option_string)
        inchikey = Chem.InchiToInchiKey(inchi_string)
    mol, mol_created = Mol.objects.get_or_create(
        defaults={'smiles': Chem.MolToSmiles(rdkit_mol),
                  'inchikey': inchikey},
        inchikey=inchikey,
        group=group)
    if mol_created:
        calculated_details = {
            'stoichiometry': Chem.CalcMolFormula(rdkit_mol),
            'mass': Chem.CalcExactMolWt(rdkit_mol),
        }
        if calc_charge:
            calculated_details.update(get_molecular_charge_props(rdkit_mol))
        mol.details = {**calculated_details, **details}
        if mol.tags is None:
            mol.tags = []
    mol.tags = list(set(mol.tags) | set(tags))
    for parent_inchikey in parent_inchikeys:
        try:
            parent_mol = Mol.objects.get(group=group, inchikey=parent_inchikey)
            mol.parents.add(parent_mol)
        except Mol.DoesNotExist:
            pass
    mol.save()
    if return_mol:
        return {'mol': mol, 'created': mol_created}
    else:
        return mol_created

def get_molecular_charge_props(rdkit_mol):
    molecular_charge = Chem.GetFormalCharge(rdkit_mol)
    num_protons = sum([a.GetAtomicNum()
                       for a in Chem.AddHs(rdkit_mol).GetAtoms()])
    electron_count = num_protons + molecular_charge
    return {'molecular_charge': molecular_charge,
            'electron_count': electron_count}


def put_smiles(smiles=None, **kwargs):
    return put_mol(rdkit_mol=Chem.MolFromSmiles(smiles), **kwargs)
