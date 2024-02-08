import os
from munch import Munch

from chemconfigs.parsers import gaussian

def theory():
    theory = Munch()
    theory.name = 'dft_hybrid_b3lyp_def2svpp_opt_gaussian'
    theory.description = "Gaussian b3lyp/def2svpp Opt"
    theory.details = Munch()
    return theory

def load_calc_list(job_dir, context=None):
    meta_data = Munch()

    out_filename = 'gaussian_test.log'
    out_path = os.path.join(job_dir, out_filename)
    with open(out_path, 'r') as f:
        lines = f.readlines()

    gaussian.assert_completed(lines)

#    program, version = orca.parse_program_version(lines)
#    program = orca.parse_program(lines)
#    meta_data.program = program
#    meta_data.version = version

#    meta_data.wall_clock_time = orca.duration(lines)
#    meta_data.host = orca.host(lines)
    

    meta_data.worker_name = "b3lyp_devsvpp_opt_gaussian"
    properties = Munch()
    properties.total_energy = gaussian.energy(lines)
#    print('Energy is '+str(properties.total_energy))
    calc = Munch(meta_data=meta_data, properties=properties)
    calc.theory = theory()


    return [calc]
