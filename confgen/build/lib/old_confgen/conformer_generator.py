# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 13:06:49 2014

@author: epyzerknapp
@author: hawk
"""

import os
import random
import subprocess

from rdkit.Chem import (AddHs, MolFromSmiles)
from rdkit.Chem.AllChem import (EmbedMultipleConfs, GetBestRMS,
                                UFFGetMoleculeForceField, MMFFGetMoleculeForceField,
                                MMFFGetMoleculeProperties)
from rdkit.Chem.rdmolops import RemoveHs


def write_xyz(coords, filename, comment):
    '''
    Write an xyz file from coords
    '''
    with open(filename, "w") as fp:
        fp.write(str(len(coords)) + "\n")
        fp.write(str(comment) + "\n")
        for atom in coords:
            fp.write("%s %.4f %.4f %.4f\n" % (atom[0], atom[1][0], atom[1][1], atom[1][2]))


def obfit_rmsd(file1, file2, smarts, path=''):
    cmd = ["obfit", str(smarts),  # Does confgen need this quotes too???
           os.path.join(path, file1 + '.xyz'),
           os.path.join(path, file2 + '.xyz')]
    ret = subprocess.check_output(cmd,
                                  stdin=None,
                                  stderr=subprocess.STDOUT,
                                  shell=False,
                                  universal_newlines=False)
    rmsd = float(ret.decode('utf-8')[5:13])
    return rmsd


def align_rmsd(file1, file2, path='', smarts=None):
    cmd = ["obabel",
           os.path.join(path, file1 + '.xyz'),
           os.path.join(path, file2 + '.xyz'),
           '-o', 'smi',
           '--align',
           '--append',
           'rmsd']
    if smarts:
        cmd += ['-s', str(smarts)]
    ret = subprocess.check_output(cmd,
                                  stdin=None,
                                  stderr=subprocess.STDOUT,
                                  shell=False,
                                  universal_newlines=False)
    rmsd = ret.decode('utf-8').split()[-4]
    return float(rmsd)


class ConformerGenerator(object):
    '''
    Generates conformations of molecules from 2D representation.
    '''

    def __init__(self, smiles, forcefield="mmff"):
        '''
        Initialises the class
        '''
        self.mol = MolFromSmiles(smiles)
        self.full_clusters = []
        self.forcefield = forcefield
        self.conf_energies = []
        self.initial_confs = None

    def generate(self, max_generated_conformers=50, prune_thresh=0.01, output=None):
        '''
        Generates conformers

        Note  the number max_generated _conformers required is related to the
        number of rotatable bonds
        '''
        self.mol = AddHs(self.mol, addCoords=True)
        self.initial_confs = EmbedMultipleConfs(self.mol,
                                                numConfs=max_generated_conformers,
                                                pruneRmsThresh=prune_thresh,
                                                useRandomCoords=False,
                                                # Despite what the documentation says -1 is a seed!!
                                                # It doesn't mean random generatrion
                                                randomSeed=random.randint(1, 10000000)
                                                )

        output.write("Generated " + str(len(self.initial_confs)) + " initial confs\n")
        return self.initial_confs

    def minimise(self, output=None):
        '''
        Minimises conformers using a force field
        '''

        if self.forcefield != "mmff" and self.forcefield != "uff":
            raise ValueError("Unrecognised force field")
        if self.forcefield == "mmff":
            props = MMFFGetMoleculeProperties(self.mol)
            for i in range(0, len(self.initial_confs)):
                output.write("Minimising conformer number {}\n".format(i))
                potential = MMFFGetMoleculeForceField(self.mol, props, confId=i)
                potential.Minimize()
                mmff_energy = potential.CalcEnergy()
                self.conf_energies.append((i, mmff_energy))

        elif self.forcefield == "uff":
            for i in range(0, len(self.initial_confs)):
                potential = UFFGetMoleculeForceField(self.mol, confId=i)
                potential.Minimize()
                uff_energy = potential.CalcEnergy()
                self.conf_energies.append((i, uff_energy))
        self.conf_energies = sorted(self.conf_energies, key=lambda tup: tup[1])
        return self.mol

    def cluster(self, rms_tolerance=0.1, max_ranked_conformers=10,
                energy_window=5, Report_e_tol=10, output=None):
        '''
        Removes duplicates after minimisation
        '''
        self.counter = 0
        self.factormax = 3
        self.mol_no_h = RemoveHs(self.mol)
        n_atom_no_h = self.mol_no_h.GetNumAtoms()
        atom_map = [(a, a) for a in range(0, n_atom_no_h)]
        calcs_performed = 0
        self.full_clusters = []
        confs = self.conf_energies[:]
        ignore = []
        ignored = 0

        for i, pair_1 in enumerate(confs):
            if i == 0:
                index_0, energy_0 = pair_1
            output.write("clustering cluster {} of {}\n".format(i, len(self.conf_energies)))
            index_1, energy_1 = pair_1
            if abs(energy_1 - energy_0) > Report_e_tol:
                msg = 'Breaking because hit Report Energy Window, E was {} Kcal/mol and minimum was {} \n'
                output.write(msg.format(energy_1, energy_0))
                break
            if i in ignore:
                ignored += i
                continue
            self.counter += 1
            if self.counter == self.factormax * max_ranked_conformers:
                output.write('Breaking because hit MaxNConfs \n')
                break
            clustered = [[self.mol.GetConformer(id=index_1), energy_1, 0.00]]
            ignore.append(i)
            for j, pair_2 in enumerate(confs):
                if j > 1:
                    index_2, energy_2 = pair_2
                    if j in ignore:
                        ignored += 1
                        continue
                    if abs(energy_1 - energy_2) > energy_window:
                        break
                    if abs(energy_1 - energy_2) <= 1e-3:
                        clustered.append([self.mol.GetConformer(id=index_2), energy_2, 0.00])
                        ignore.append(j)
                    else:
                        rms = GetBestRMS(self.mol_no_h,
                                         self.mol_no_h,
                                         refId=index_1,
                                         prbId=index_2,
                                         map=[atom_map])
                        calcs_performed += 1
                        if rms <= rms_tolerance:
                            clustered.append([self.mol.GetConformer(id=index_2), energy_2, rms])
                            ignore.append(j)
            self.full_clusters.append(clustered)
        output.write("{} ignore passes made\n".format(ignored))
        output.write("{} overlays needed out of a possible {}\n".format(calcs_performed, len(self.conf_energies) ** 2))
        # print self.full_clusters
        ranked_clusters = []
        for i, cluster in enumerate(self.full_clusters):
            if i < self.factormax * max_ranked_conformers:
                ranked_clusters.append(cluster[0])
        # print(ranked_clusters)
        return ranked_clusters

    def recluster(self,
                  rms_tolerance=0.1,
                  max_ranked_conformers=10,
                  energy_window=5,
                  output=None,
                  clustered_confs=None,
                  molecule=None,
                  key=None,
                  fallback_to_align=False):

        self.removed = []
        self.counter = 0
        i = -1
        # print(clustered_confs)
        for conformerA in clustered_confs:
            i += 1
            j = i
            if self.counter == max_ranked_conformers:
                for k in range(i, len(clustered_confs)):
                    if os.path.isfile(key + "_Conf_" + str(k + 1) + ".xyz"):
                        os.remove(key + "_Conf_" + str(k + 1) + ".xyz")
                        output.write("removed" + key + "_Conf_" + str(k + 1) + ".xyz\n")
                break
            if i in self.removed:
                continue
            self.counter += 1
            for conformerB in clustered_confs[i + 1:]:
                j += 1
                if conformerB[1] - conformerA[1] > energy_window:
                    break
                if j in self.removed:
                    continue
                try:
                    rms = obfit_rmsd(key + "_Conf_" + str(i + 1), key + "_Conf_" + str(j + 1), str(molecule))
                except subprocess.CalledProcessError:
                    if fallback_to_align:
                        output.write('obfit failed, falling back to obabel --align')
                        rms = align_rmsd(key + "_Conf_" + str(i + 1), key + "_Conf_" + str(j + 1))
                    else:
                        raise

                output.write("Comparing " + str(i + 1) + " " + str(j + 1) + ' RMSD ' + str(rms) + "\n")
                if rms > rms_tolerance:
                    pos = _atomic_pos_from_conformer(conformerB[0])
                    elements = _extract_atomic_type(conformerB[0])
                    pos = [[-float(coor[k]) for k in range(3)] for coor in pos]
                    coords = list(zip(elements, pos))
                    write_xyz(coords=coords, filename=key + "_Conf_" + str(j + 1) + "_inv.xyz", comment=conformerB[1])
                    try:
                        file1 = key + "_Conf_" + str(i + 1)
                        file2 = key + "_Conf_" + str(j + 1) + "_inv"
                        rmsinv = obfit_rmsd(file1, file2, str(molecule))
                    except subprocess.CalledProcessError:
                        if fallback_to_align:
                            output.write('obfit failed, falling back to obabel --align')
                            rmsinv = align_rmsd(key + "_Conf_" + str(i + 1), key + "_Conf_" + str(j + 1) + "_inv")
                        else:
                            raise

                    rms = min([rms, rmsinv])
                    os.remove(key + "_Conf_" + str(j + 1) + "_inv.xyz")
                    output.write("Comparing {} {} RMSD after checking inversion {}\n".format(i + 1, j + 1, rms))
                if rms <= rms_tolerance:
                    self.removed.append(j)
                    output.write("Removed Conf_" + str(j + 1) + "\n")
                    # print 'would delete', key + "_Conf_"+str(j+1)+"_inv.xyz", key + "_Conf_"+str(j+1)+".xyz"
                    os.remove(key + "_Conf_" + str(j + 1) + ".xyz")

    @staticmethod
    def cluster_clean_dftb(rms_tolerance=0.1,
                           max_ranked_conformers=10,
                           energy_window=5,
                           output=None,
                           molecule=None,
                           key=None,
                           fallback_to_align=False,
                           ):
        """
        Additional clean up pass using DFTB+. Reads xyz files present on the disk, minimizes, sorts by energy
        and makes same procedure for cleaning as 'recluster' above
        """
        import ase.io
        from ase.calculators.dftb import Dftb
        from ase.optimize import QuasiNewton
        import time
        import glob

        start_time = time.time()
        output.write("Cleanup conformers using DFTB:\n")

        flist = glob.glob("./*.xyz")
        clustered_confs = [None] * len(flist)
        for filename in flist:
            num = int(filename.split('_')[-1][:-4])
            output.write("  Minimising conformer number {}\n".format(num))
            test = ase.io.read(filename)
            # Generate MaxAngular dictionary for dftb input file based on elements
            max_ang = {}
            for elem in test.get_chemical_symbols():
                if elem == "H":
                    max_ang['Hamiltonian_MaxAngularMomentum_' + str(elem)] = '"s"'
                else:
                    max_ang['Hamiltonian_MaxAngularMomentum_' + str(elem)] = '"p"'

            test.set_calculator(Dftb(label='conf',
                                     atoms=test,
                                     run_manyDftb_steps=True,
                                     Driver_='ConjugateGradient',
                                     Driver_MaxForceComponent='1E-4',
                                     Driver_MaxSteps=1000,
                                     Hamiltonian_SCC='Yes',
                                     Hamiltonian_MaxAngularMomentum_='',
                                     **max_ang
                                     ))

            dyn = QuasiNewton(test, trajectory='test.traj')
            try:
                dyn.run(fmax=10, steps=1)
            except Exception:
                output.write(
                    "  DFTB error, probably you have an atom not supported by parameters (see conf.out for details)\n")
                break

            totalenergy = test.get_total_energy()
            output.write("  Energy of conformer number {} is {}\n".format(num, totalenergy))
            # the optimized structure is in dftb outbup file gen and not in ase.Atom (bug?)
            test = ase.io.read('geo_end.gen')
            # Emulate clustered_confs struct for further reclustering, instead conformer
            clustered_confs[num - 1] = [(test.get_chemical_symbols(), test.get_positions()), totalenergy, 0.0]
        else:

            # Clean up old MM xyz files
            for filename in glob.glob("./*.xyz"):
                os.remove(filename)

            # sort conformeres by DFTB energy
            clustered_confs.sort(key=lambda conformer: conformer[1])
            # write conformers xyz files to disk (sorted by lowest energy)
            for i, conformer in enumerate(clustered_confs):
                (elements, pos) = conformer[0]
                coords = list(zip(elements, pos))
                write_xyz(coords=coords, filename=key + "_Conf_" + str(i + 1) + ".xyz", comment=conformer[1])

            # almost the same as recluster method, only modified to use (elements, pos) tuple instead of Conformer
            removed = []
            counter = 0
            i = -1
            # print(clustered_confs)
            for conformerA in clustered_confs:
                i += 1
                j = i
                if counter == max_ranked_conformers:
                    for k in range(i, len(clustered_confs)):
                        if os.path.isfile(key + "_Conf_" + str(k + 1) + ".xyz"):
                            os.remove(key + "_Conf_" + str(k + 1) + ".xyz")
                            output.write("  reached max confs removed" + key + "_Conf_" + str(k + 1) + ".xyz\n")
                    break
                if i in removed:
                    continue
                counter += 1
                for conformerB in clustered_confs[i + 1:]:
                    j += 1
                    if conformerB[1] - conformerA[1] > energy_window: #*0.043: convert to eV fom kcal/mol
                        break
                    if j in removed:
                        continue
                    try:
                        rms = obfit_rmsd(key + "_Conf_" + str(i + 1), key + "_Conf_" + str(j + 1), str(molecule))
                    except subprocess.CalledProcessError:
                        if fallback_to_align:
                            output.write('  obfit failed, falling back to obabel --align')
                            rms = align_rmsd(key + "_Conf_" + str(i + 1), key + "_Conf_" + str(j + 1))
                        else:
                            raise

                    output.write("  Comparing " + str(i + 1) + " " + str(j + 1) + ' RMSD ' + str(rms) + "\n")
                    if rms > rms_tolerance:
                        (elements, pos) = conformerB[0]
                        pos = [[-float(coor[k]) for k in range(3)] for coor in pos]
                        coords = list(zip(elements, pos))
                        write_xyz(coords=coords, filename=key + "_Conf_" + str(j + 1) + "_inv.xyz",
                                  comment=conformerB[1])
                        try:
                            file1 = key + "_Conf_" + str(i + 1)
                            file2 = key + "_Conf_" + str(j + 1) + "_inv"
                            rmsinv = obfit_rmsd(file1, file2, str(molecule))
                        except subprocess.CalledProcessError:
                            if fallback_to_align:
                                output.write('  obfit failed, falling back to obabel --align')
                                rmsinv = align_rmsd(key + "_Conf_" + str(i + 1), key + "_Conf_" + str(j + 1) + "_inv")
                            else:
                                raise

                        rms = min([rms, rmsinv])
                        os.remove(key + "_Conf_" + str(j + 1) + "_inv.xyz")
                        output.write("  Comparing {} {} RMSD after checking inversion {}\n".format(i + 1, j + 1, rms))

                    if rms <= rms_tolerance:
                        removed.append(j)
                        output.write("  Removed Conf_" + str(j + 1) + "\n")
                        # print 'would delete', key + "_Conf_"+str(j+1)+"_inv.xyz", key + "_Conf_"+str(j+1)+".xyz"
                        os.remove(key + "_Conf_" + str(j + 1) + ".xyz")

        output.write('Total time DFTB clean {0:1f}  sec\n'.format(time.time() - start_time))


def _extract_atomic_type(confomer):
    """
    Extracts the elements associated with a conformer, in order that prune_threshy
    are read in
    """
    elements = []
    mol = confomer.GetOwningMol()
    for atom in mol.GetAtoms():
        elements.append(atom.GetSymbol())
    return elements


def _atomic_pos_from_conformer(conformer):
    """
    Extracts the atomic positions for an RDKit conformer object, to allow writing
    of input files, uploading to databases, etc.
    Returns a list of lists
    """
    atom_positions = []
    natoms = conformer.GetNumAtoms()
    for atom_num in range(0, natoms):
        pos = conformer.GetAtomPosition(atom_num)
        atom_positions.append([pos.x, pos.y, pos.z])
    return atom_positions


def rename_xyz_files(path="."):
    namedict = {}
    flist = os.listdir(path)
    for filename in flist:
        if filename.endswith(".xyz") and ("Conf" in filename):
            num = int(filename.split('_')[-1][:-4])
            namedict[num] = filename
    keys = namedict.keys()
    for i, num in enumerate(sorted(keys)):
        oldfilename = namedict[num]
        newfilename = '_'.join(oldfilename.split('_')[:-1]) + '_' + str(i + 1) + ".xyz"
        oldfilepath = os.path.join(path, oldfilename)
        newfilepath = os.path.join(path, newfilename)
        os.rename(oldfilepath, newfilepath)
