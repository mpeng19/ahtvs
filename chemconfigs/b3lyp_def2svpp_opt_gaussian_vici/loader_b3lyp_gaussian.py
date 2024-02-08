import os
from munch import Munch
import pprint

from chemconfigs.parsers import gaussian

def theory():
    theory = Munch()
    theory.name = 'dft_hybrid_b3lyp_Def2SVPP_gaussian'
    theory.description = "Gaussian b3lyp/Def2SVPP opt DFT"
    theory.details = Munch()
    
    return theory

def load_calc_list(job_dir, context=None):
    meta_data = Munch()

    out_filename = 'gaussian_test.log'
    out_path = os.path.join(job_dir, out_filename)
    with open(out_path, 'r') as f:
        lines = f.readlines()

    gaussian.assert_completed(lines)

    if '    -- Stationary point found.\n' not in lines:
        raise AssertionError("'-- Stationary point found.' line not found in Gaussian calc output")


#    program, version = orca.parse_program_version(lines)
#    program = orca.parse_program(lines)
#    meta_data.program = program
#    meta_data.version = version

#    meta_data.wall_clock_time = orca.duration(lines)
#    meta_data.host = orca.host(lines)
    

    meta_data.worker_name = "b3lyp_def2svpp_opt_gaussian_vici"
    properties = Munch()
    properties.total_energy = gaussian.energy(lines)
#    print('Energy is '+str(properties.total_energy))


    properties.electric_dipole_moment_norm = gaussian.dipole(lines)
    calc = Munch(meta_data=meta_data, properties=properties)
    calc.theory = theory()
    calc.properties.update(gaussian.parse_orbitals(lines))

    calc.coords = gaussian.last_coordinate_list(lines)

    return [calc]

#if __name__ == "__main__":
#    calc_info = load_calc_list('.','')
#    pprint.pprint(calc_info)

