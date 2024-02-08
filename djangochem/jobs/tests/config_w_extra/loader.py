import os

from munch import Munch

XYZ_NAME = "{}_Conf_{}.xyz"
LOAD_MAX_CONFORMERS = 1


def theory():
    theory = Munch()
    theory.name = "mmff94"
    theory.description = "rdkit MMFF conformer: "
    return theory


def load_calc_list(job_dir, context):
    for i in range(LOAD_MAX_CONFORMERS):
        path = os.path.join(job_dir, XYZ_NAME.format(context['inchikey'], i+1))
        calc_list = []
        if os.path.isfile(path):
            meta_data = Munch()
            meta_data.program = "rdkit"
            meta_data.version = "2015_03_1"
            meta_data.data_type = "clc_cnf"
            meta_data.host = "localhost"
            meta_data.generation = 0
            calc = Munch(meta_data=meta_data,
                         properties=Munch())

            calc.theory = theory()
            lines = open(path, 'r').readlines()
            calc.properties.total_energy = float(lines[1])
            calc.coords = []
            for line in lines[2:]:
                atom_char, x, y, z = line.split()
                calc.coords.append(dict(element=atom_char, x=x, y=y, z=z))
            calc_list.append(calc)
    return calc_list
