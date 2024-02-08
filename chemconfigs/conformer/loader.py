import os
from munch import Munch

XYZ_NAME = "{0}_Conf_{1}.xyz"
LOAD_MAX_CONFORMERS = 50

WORKER = "conformer"
PROGRAM = "confgen (1-06-19)"
VERSION = "rdkit-2019_03_1"


def theory():
    theory = Munch()
    theory.name = 'molecular_mechanics_mmff94_or_uff'
    theory.description = "MMFF or UFF conformer."
    return theory


def load_calc_list(job_dir, kwargs):
    geom_list = []
    old_path = os.path.join(job_dir, "log.log")
    if os.path.exists(old_path):
        logpath = old_path
    else:  # setup new name
        logpath = os.path.join(job_dir, "confgen.log")
    with open(logpath) as fp:
        loglines = fp.readlines()
    if loglines[-1].strip() != "Terminated successfully":
        raise Exception("'Terminated successfully' not found at end of conformer output")
    duration = float(loglines[-2].split()[2])
    host = loglines[-7]
    for i in range(LOAD_MAX_CONFORMERS):
        path = os.path.join(job_dir, XYZ_NAME.format(kwargs["inchikey"], i+1))
        if os.path.isfile(path):
            meta_data = Munch()

            meta_data.generation = 0
            meta_data.program = PROGRAM
            meta_data.version = VERSION
            meta_data.wall_clock_time = duration
            meta_data.host = host

            meta_data.worker_name = WORKER
            meta_data.tag_list = ["mmff_conf_{0}".format(i+1)]
            geom = Munch(meta_data=meta_data,
                         properties=Munch())
            geom.details = {'confnum': i+1}
            geom.theory = theory()
            lines = open(path, 'r').readlines()
            geom.properties.total_energy = float(lines[1])
            geom.coords = []
            for line in lines[2:]:
                atom_char, x, y, z = line.split()
                geom.coords.append(dict(element=atom_char, x=x, y=y, z=z))
            geom_list.append(geom)
    return geom_list

