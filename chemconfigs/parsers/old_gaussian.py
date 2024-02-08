import itertools
import re
from elements import elements

COORD_MARKER = "Standard orientation:"
COORD_ENDING_MARKER = "---------------------------------------------------------------------"

#ORB_MARKER = "The electronic state is"
ORB_MARKER = "Alpha  occ. eigenvalues"
ORB_INT_MARKER = "Alpha virt. eigenvalues"
#ORB_ENDING_MARKER = "-" * 62
ORB_ENDING_MARKER = "Condensed to atoms"

def parse_program(lines):
    for line in lines:
            if "Program Version" in line:
                text = line.split()
                program = "ORCA (" + " ".join(text[:3]) + ")"
#                version = text[5]
                return program #, version


def assert_completed(lines):
    for line in reversed(lines):
        if 'Normal termination of Gaussian' in line:
            return
    raise AssertionError("Normal termination of Gaussian line not found in Gaussian calc output")

def energy(lines):
    for line in reversed(lines):
        if 'SCF Done:' in line:
#            print(re.findall("[-+]?\d+\.\d+", line))
            return float(re.findall("[-+]?\d+\.\d+", line)[0])

    raise AssertionError("'Energy is ' not found in Gaussian calc output")



def parse_coords_starting_at_line(lines, start_line_num):
    data=elements.Elements #get all elements
    data=sorted(data,key=lambda i:i.AtomicNumber) #sort elements by atomic number in ascending order
#    print(data[1].Symbol)
    cur_line = start_line_num
    coords = []
    line = lines[cur_line]
    while COORD_ENDING_MARKER not in line:
        atom_num = line.split()[1]
#        print(str(atom_num))
        atom_char = data[int(atom_num)-1].Symbol
#        print(re.findall("\-?\d+(?:\.\d+)?", line))
        num, empty_var,empty_var2, x, y, z = re.findall("\-?\d+(?:\.\d+)?", line)
        coords.append(dict(element=atom_char, x=x, y=y, z=z))
        cur_line += 1
        line = lines[cur_line]
    return coords


def last_coordinate_list(lines):
    j = 0
    coord_marker_set = set()
    for line in reversed(lines):
        j = j + 1
        if COORD_MARKER in line:
#            print('Found COORD_MARKER')
            normal_start = len(lines) - j
            coord_marker_set.add(normal_start)
    
#    print('Normal Start is '+str(max(coord_marker_set)))
    return parse_coords_starting_at_line(lines, max(coord_marker_set)+5)


def parse_orbitals_from_line(lines, cur_line, method, spin, orb_type, stop_marker):
    orb_list = []
    orb_counter = 0
    start_line = cur_line
    line = lines[cur_line]
    while stop_marker not in line:
        for energy in re.findall("\-?\d+(?:\.\d+)?", line):
            orb_counter += 1
            d = dict(method=method, spin=spin, type=orb_type, number=orb_counter, energy=float(energy))
            orb_list.append(d)
        cur_line += 1
        try:
            line = lines[cur_line]
        except IndexError:
            raise IndexError("Run out of lines looking for {} after line {}".format(stop_marker, start_line))
    return orb_list, cur_line

def parse_orbitals(lines, multiplicity=1):
    cur_line = 0
    orb_list = []
    homo, lumo, somo = None, None, None

    line = lines[cur_line]
    while ORB_MARKER not in line:
        cur_line += 1
        line = lines[cur_line]

    if multiplicity == 1:
        occ_orbs, cur_line = parse_orbitals_from_line(lines,
                                                      cur_line,
                                                      method="rhf",
                                                      spin='alpha',
                                                      orb_type="occ",
                                                      stop_marker=ORB_INT_MARKER)
        orb_list.extend(occ_orbs)
        unocc_orbs, cur_line = parse_orbitals_from_line(lines,
                                                        cur_line,
                                                        method="rhf",
                                                        spin='alpha',
                                                        orb_type="unocc",
                                                        stop_marker=ORB_ENDING_MARKER)
        orb_list.extend(unocc_orbs)
        homo = float(occ_orbs[-1]["energy"])
        lumo = float(unocc_orbs[0]["energy"])
        return {'orb_list': orb_list, "homo": homo, "lumo": lumo}

    else:
        occ_orbs, cur_line = parse_orbitals_from_line(lines,
                                                      cur_line,
                                                      method="uhf",
                                                      spin='alpha',
                                                      orb_type="occ",
                                                      stop_marker=ORB_INT_MARKER)
        orb_list.extend(occ_orbs)
        unocc_orbs, cur_line = parse_orbitals_from_line(lines,
                                                        cur_line,
                                                        method="uhf",
                                                        spin='alpha',
                                                        orb_type="unocc",
                                                        stop_marker=ORB_MARKER)
        orb_list.extend(unocc_orbs)

        somo = float(occ_orbs[-1]["energy"])

        occ_orbs, cur_line = parse_orbitals_from_line(lines,
                                                      cur_line,
                                                      method="uhf",
                                                      spin='beta',
                                                      orb_type="occ",
                                                      stop_marker=ORB_INT_MARKER)
        orb_list.extend(occ_orbs)
        unocc_orbs, cur_line = parse_orbitals_from_line(lines,
                                                        cur_line,
                                                        method="uhf",
                                                        spin='beta',
                                                        orb_type="unocc",
                                                        stop_marker=ORB_ENDING_MARKER)
        orb_list.extend(unocc_orbs)
        return {'orb_list': orb_list, "somo": somo}





def dipole(lines):
    prev_line = None
    for line in reversed(lines):
        if "Dipole moment (field-independent basis, Debye):" in line:
#            print('prev line is')
#            print(prev_line)
            return float(re.findall("(-?\d+(?:\.\d+)?)", prev_line)[3])
        prev_line = line
    raise AssertionError("' Dipole Moment ' not found in Gaussian calc output")

