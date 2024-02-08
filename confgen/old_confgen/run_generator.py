import argparse
import os
import re
import socket
import subprocess
import time

from rdkit.Chem import MolFromSmiles

from confgen.conformer_generator import (ConformerGenerator,
                                 _extract_atomic_type,
                                 _atomic_pos_from_conformer,
                                 rename_xyz_files,
                                 write_xyz
                                )

def main(smiles_file=None,
         forcefield=None,
         nconf=None,
         nconf_gen=None,
         e_window=None,
         rms_tol=None,
         prun_tol=None,
         interactive=None,
         log=None,
         rep_e_window=None,
         fallback_to_align=None,
         clean_DFTB=None):
    with open(log, "w") as output:
        smiles = []
        if not interactive:
            smiles_f = open(smiles_file)
            for line in smiles_f:
                if len(line.split()) != 0:
                    smiles.append(line.split()[0])
        else:
            for arg in args:
                smiles.append(arg)

        output.write("The smiles strings that will be run are:\n")
        output.write("\n".join(smiles)+"\n")

        for j, molecule in enumerate(smiles):
            if 'B' in molecule:
                output.write("Switching to UFF, since MMFF94 does not have boron\n")
                forcefield = 'uff'

            output.write("Analysing smiles string {}\n".format(molecule))
            MolFromSmiles(molecule)
            # print "There are", NumRotatableBonds(mol)
            output.write("Generating initial conformations\n")
            confgen = ConformerGenerator(smiles=molecule, forcefield=forcefield)
            output.write("Minimising conformations using the {} force field\n".format(forcefield))
            confgen.generate(max_generated_conformers=int(nconf_gen),
                             prune_thresh=float(prun_tol),
                             output=output)
            gen_time = time.time()
            confgen.minimise(output=output)
            min_time = time.time()
            output.write("Minimisation complete, generated conformations with the following energies:\n")
            output.write("\n".join([str(energy[1]) for energy in confgen.conf_energies])+"\n")
            msg = ("Clustering structures using an energy window of {} and a rms "
                   "tolerance of {} and a Report Energy Window of {}\n")
            output.write(msg.format(float(e_window), float(rms_tol), float(rep_e_window)))
            clustered_confs = confgen.cluster(rms_tolerance=float(rms_tol),
                                              max_ranked_conformers=int(nconf),
                                              energy_window=float(e_window),
                                              Report_e_tol=float(rep_e_window),
                                              output=output)

            cluster_time = time.time()
            key = os.path.basename(str(smiles_file)[:-4])
            for i, conformer in enumerate(clustered_confs):
                output.write("Cluster {} has energy {}\n".format(i, conformer[1]))
                pos = _atomic_pos_from_conformer(conformer[0])
                elements = _extract_atomic_type(conformer[0])
                coords = list(zip(elements, pos))
                xyz_file = "{key}_Conf_{index}.xyz".format(
                    key=key, index=(i + 1))
                write_xyz(coords=coords, filename=xyz_file,
                          comment=conformer[1])

            cmd = ["obabel", key + "_Conf_"+str(i+1)+".xyz", "-osmi"]
            molecule = subprocess.check_output(cmd,
                                               stdin=None,
                                               stderr=subprocess.STDOUT,
                                               shell=False,
                                               universal_newlines=False).decode('utf-8')
            molecule = str(molecule.split()[0])
            molecule = re.sub('C', '[#6]', molecule)
            molecule = re.sub('c', '[#6]', molecule)
            molecule = re.sub('N', '[#7]', molecule)
            molecule = re.sub('n', '[#7]', molecule)
            molecule = re.sub('\[\[', '[', molecule)
            molecule = re.sub('\]\]', ']', molecule)
            molecule = re.sub('\]H\]', 'H]', molecule)
            molecule = re.sub('=', '~', molecule)

            confgen.recluster(rms_tolerance=float(rms_tol),
                              max_ranked_conformers=int(nconf),
                              energy_window=float(e_window),
                              output=output,
                              clustered_confs=clustered_confs,
                              molecule=molecule,
                              key=key,
                              fallback_to_align=fallback_to_align)

            rename_xyz_files()
            recluster_time = time.time()

            if clean_DFTB:
                confgen.cluster_clean_dftb(max_ranked_conformers=int(nconf),
                                           rms_tolerance=float(rms_tol),
                                           energy_window=float(e_window),
                                           output=output,
                                           molecule=molecule,
                                           key=key,
                                           fallback_to_align=fallback_to_align)
                rename_xyz_files()

            output.write(socket.gethostname()+"\n")
            output.write('gen time  {0:1f}  sec\n'.format(gen_time - start_time))
            output.write('min time  {0:1f}  sec\n'.format(min_time - gen_time))
            output.write('cluster time  {0:1f}  sec\n'.format(cluster_time - min_time))
            output.write('recluster time  {0:1f}  sec\n'.format(recluster_time - cluster_time))
            output.write('total time  {0:1f}  sec\n'.format(time.time() - start_time))
            output.write('Terminated successfully\n')


if __name__ == "__main__":
    start_time = time.time()

    parser = argparse.ArgumentParser(description="Runs conformer generator")
    parser.add_argument('smiles_file', help="file containing smiles string")
    parser.add_argument('-f', '--forcefield',
                        help="Forcefield used for minimisations",
                        default="mmff")
    parser.add_argument('-n',
                        '--minimised_conformations',
                        dest='nconf',
                        help="Number of low energy conformations to return",
                        default=20)
    parser.add_argument('-g',
                        '--generated_confs',
                        dest="nconf_gen",
                        help="Number of conformations to build in the generation stage",
                        default=200)
    parser.add_argument('-e', '--energy_window',
                        dest="e_window",
                        help="Energy window to use when clustering, Kcal/mol",
                        default=5.0)
    parser.add_argument('-t', '--rms_tolerance',
                        dest="rms_tol",
                        help="RMS tolerance to use when clustering",
                        default=0.1)
    parser.add_argument('-p',
                        '--pruning_tolerance',
                        dest="prun_tol",
                        help="RMS tolerance used in pruning in the generation phase",
                        default=0.01)
    parser.add_argument('-i', '--interactive',
                        dest="interactive",
                        help="Interactive session, expect smiles string as argument instead of smiles file",
                        default=False,
                        action="store_true")
    parser.add_argument("-l",
                        "--log",
                        dest="log",
                        help="name of log file",
                        default="confgen.log")
    parser.add_argument("-E", '--report_energy_window',
                        dest="rep_e_window",
                        help="Energy window to use when reporting, Kcal/mol",
                        default=5.0)
    parser.add_argument('--fallback_to_align',
                        action='store_true',
                        help="whether to use rmsd_align in obfit fails")
    parser.add_argument('--clean_DFTB',
                        action='store_true',
                        default=False,
                        help='run DFTB confomer minimization at the end to clean conformers')

    parser.set_defaults(fallback_to_align=False)

    args = parser.parse_args()
    main(**vars(args))
