from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import AllChem
from operator import itemgetter
import time

def write_xyz(coords, filename, comment):
    '''
    Write an xyz file from coords
    '''
    with open(filename, "w") as fp:
        fp.write(str(len(coords))+"\n")
        fp.write(str(comment)+"\n")
        for atom in coords:
            fp.write("%s %.4f %.4f %.4f\n" % (atom[0], atom[1][0], atom[1][1], atom[1][2]))


def _extract_atomic_type(confomer):
    '''
    Extracts the elements associated with a conformer, in order that prune_threshy
    are read in
    '''
    elements = []
    mol = confomer.GetOwningMol()
    for atom in mol.GetAtoms():
        elements.append(atom.GetSymbol())
    return elements


def _atomic_pos_from_conformer(conformer):
    '''
    Extracts the atomic positions for an RDKit conformer object, to allow writing
    of input files, uploading to databases, etc.
    Returns a list of lists
    '''
    atom_positions = []
    natoms = conformer.GetNumAtoms()
    for atom_num in range(0, natoms):
        pos = conformer.GetAtomPosition(atom_num)
        atom_positions.append([pos.x, pos.y, pos.z])
    return atom_positions

# The "main" in python starts here 
if __name__ == "__main__":
    start = time.time()
    input_smiles_file = open('input_smiles.dat','r') 
    input_smiles_file_lines = input_smiles_file.readlines()
    smiles = input_smiles_file_lines[0]
    print('Input smiles: ',smiles)
    m = Chem.MolFromSmiles(smiles)
    m_with_h=Chem.AddHs(m)
    
    # Could make an input file that adjusts all of these parameters
    prune_thresh = 0.01 # RMSD where 2 conformers are considered identical...
    max_energy_window = 10.0 # Energy window to take conformers from (kcal/mol)
    initial_conformers = 200
    max_conformers_to_return = 25

    conformer_ids = AllChem.EmbedMultipleConfs(m_with_h, initial_conformers,AllChem.ETKDG())
    
    props = AllChem.MMFFGetMoleculeProperties(m_with_h)

    pre_opt_energy_dict = {}
    for cid in conformer_ids:
        force_field = AllChem.UFFGetMoleculeForceField(m_with_h, confId=cid)
#        force_field = AllChem.UFFGetMoleculeForceField(m_with_h, props, confId=cid)
        energy = force_field.CalcEnergy()
        pre_opt_energy_dict[str(cid)]  = energy
    
    props = AllChem.MMFFGetMoleculeProperties(m_with_h)

    counter = 1
    for cid in conformer_ids:
        if counter % 100 == 0:
            print('Minimizing conformer',str(counter))
#        AllChem.MMFFOptimizeMolecule(m_with_h, confId=cid)
        AllChem.UFFOptimizeMolecule(m_with_h, confId=cid)
        counter +=1

    energy_dict = {} 
    energy_as_key_dicts = list()
    opt_energies = []
    for cid in conformer_ids:
#        force_field = AllChem.MMFFGetMoleculeForceField(m_with_h, props, confId=cid)
        force_field = AllChem.UFFGetMoleculeForceField(m_with_h, confId=cid)
        energy = force_field.CalcEnergy()
        energy_dict[str(cid)]  = energy
        energy_as_key_dict = {'energy': energy, 'cid': cid}
        energy_as_key_dicts.append(energy_as_key_dict)
        opt_energies.append(energy)


    min_energy = min(opt_energies)
    print('Min Energy is',str(min_energy))

    print('Removing conformers above the energy threshold')
    low_energy_conformers = list()
    for cid in conformer_ids:
        if energy_dict[str(cid)] < min_energy + max_energy_window:
            low_energy_conformers.append(cid)


    print('Number of remaining conformers ',str(len(low_energy_conformers)))
    print('Writing out lowest ', str(max_conformers_to_return), 'conformers with clustering')
    print('Hydrogens are removed during this process for speed')    

    final_conf_id_list = list()
    excluded_confs = list()

    energy_as_key_dicts = sorted(energy_as_key_dicts,key=itemgetter('energy'))
    num_confs_added = 0 
    mol_no_h = Chem.RemoveHs(m_with_h)
    n_atom_no_h = mol_no_h.GetNumAtoms()
    atom_map = [(a, a) for a in range(0, n_atom_no_h)]

    for energy_as_key_dict in energy_as_key_dicts:
        has_low_rms_match = False
        for conf_id in final_conf_id_list:
            # This can be made more efficient, on the to-do list
            rms =  AllChem.GetBestRMS(prbMol=mol_no_h,refMol=mol_no_h,refId=energy_as_key_dict['cid'],prbId=conf_id)
            if rms < prune_thresh:
                print('Conformers',str(conf_id),'and ',str(energy_as_key_dict['cid']),'have an rms of',rms)
                has_low_rms_match = True

        if not has_low_rms_match: 
            num_confs_added += 1 
            final_conf_id_list.append(energy_as_key_dict['cid'])
        
        if num_confs_added >= max_conformers_to_return:
            print('Found ', str(max_conformers_to_return), 'breaking')
            break
    
    print('The total conformer list after clustering ',str(len(final_conf_id_list)))

    counter = 0
    for conf_id in final_conf_id_list:
        conformer = m_with_h.GetConformer(id=conf_id)
        pos = _atomic_pos_from_conformer(conformer)
        elements = _extract_atomic_type(conformer)
        pos = [[-float(coor[k]) for k in range(3)] for coor in pos]
        coords = list(zip(elements, pos))
        xyz_file = "Conf_{index}.xyz".format(index=(counter + 1))
        write_xyz(coords=coords, filename=xyz_file,
                              comment=smiles)
        counter += 1 




    print('Program completed')
    end = time.time()
    print('Total time (s):',end-start)
