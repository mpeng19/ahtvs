import itertools
import re
from elements import elements

COORD_MARKER = "Standard orientation:"
COORD_ENDING_MARKER = "---------------------------------------------------------------------"

#ORB_MARKER = "The electronic state is"
ORB_START_MARKER = "The electronic state is"
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
    # count orb markers
    cur_line = 0

    ORB_MARKER_TOTAL_COUNTER = 0

    for line in lines:
        if ORB_START_MARKER in line:
            ORB_MARKER_TOTAL_COUNTER += 1


    cur_line = 0
    orb_list = []
    homo, lumo, somo = None, None, None

    line = lines[cur_line]

    ORB_MARKER_COUNTER = 0

    while ORB_MARKER_COUNTER < ORB_MARKER_TOTAL_COUNTER:
        cur_line += 1
        line = lines[cur_line]
        if ORB_START_MARKER in line:
            ORB_MARKER_COUNTER += 1

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



def electronic_spatial_extent(lines):
    for line in reversed(lines):
        if "Electronic spatial extent" in line:
            return float(re.findall("(-?\d+(?:\.\d+)?)", line)[1])
    raise AssertionError("' Dipole Moment ' not found in Gaussian calc output")



def dipole(lines):
    dipole = dict()
    prev_line = None
    for line in reversed(lines):
        if "Dipole moment (field-independent basis, Debye):" in line:
#            print('prev line is')
#            print(prev_line)
            list_of_dipole_contribs = re.findall("(-?\d+(?:\.\d+)?)", prev_line)
            dipole['x'] = list_of_dipole_contribs[0]
            dipole['y'] = list_of_dipole_contribs[1]
            dipole['z'] = list_of_dipole_contribs[2]
            dipole['electric_dipole_moment_norm'] = list_of_dipole_contribs[3]

            return dipole
#            return float(re.findall("(-?\d+(?:\.\d+)?)", prev_line)[3])
        prev_line = line
    raise AssertionError("' Dipole Moment ' not found in Gaussian calc output")

def quadrupole_traceless(lines):
    quadrupole_traceless = dict()
    prev_line_1 = None
    prev_line_2 = None
    for line in reversed(lines):
        if "Traceless Quadrupole moment (field-independent basis, Debye-Ang):" in line:
#            print('prev line is')
#            print(prev_line)
            list_of_traceless_quadrupole_contribs_1 = re.findall("(-?\d+(?:\.\d+)?)", prev_line_1) # XX YY ZZ
            list_of_traceless_quadrupole_contribs_2 = re.findall("(-?\d+(?:\.\d+)?)", prev_line_2) # XY XZ YX

            quadrupole_traceless['XX'] = list_of_traceless_quadrupole_contribs_1[0]
            quadrupole_traceless['YY'] = list_of_traceless_quadrupole_contribs_1[1]
            quadrupole_traceless['ZZ'] = list_of_traceless_quadrupole_contribs_1[2]

            quadrupole_traceless['XY'] = list_of_traceless_quadrupole_contribs_2[0]
            quadrupole_traceless['XZ'] = list_of_traceless_quadrupole_contribs_2[1]
            quadrupole_traceless['YZ'] = list_of_traceless_quadrupole_contribs_2[2]

            quadrupole_traceless['units'] = "Debye-Ang"
            return quadrupole_traceless
#            return float(re.findall("(-?\d+(?:\.\d+)?)", prev_line)[3])
        prev_line_2 = prev_line_1
        prev_line_1 = line

    raise AssertionError("' Traceless Quadrupole moment ' not found in Gaussian calc output")


def octapole(lines):
    octapole = dict()
    prev_line_1 = None
    prev_line_2 = None
    prev_line_3 = None
    for line in reversed(lines):
        if "Octapole moment (field-independent basis, Debye-Ang**2)" in line:
#            print('prev line is')
#            print(prev_line)
            list_of_octapole_contribs_1 = re.findall("(-?\d+(?:\.\d+)?)", prev_line_1) #  XXX  YYY ZZZ XYY
            list_of_octapole_contribs_2 = re.findall("(-?\d+(?:\.\d+)?)", prev_line_2) # XXY XXZ  XZZ  YZZ
            list_of_octapole_contribs_3 = re.findall("(-?\d+(?:\.\d+)?)", prev_line_3) # YYZ  XYZ

#            print(prev_line_1)
#            print(prev_line_2)
#            print(prev_line_3)

            octapole['XXX'] = list_of_octapole_contribs_1[0]
            octapole['YYY'] = list_of_octapole_contribs_1[1]
            octapole['ZZZ'] = list_of_octapole_contribs_1[2]
            octapole['XYY'] = list_of_octapole_contribs_1[3]

            octapole['XXY'] = list_of_octapole_contribs_2[0]
            octapole['XXZ'] = list_of_octapole_contribs_2[1]
            octapole['XZZ'] = list_of_octapole_contribs_2[2]
            octapole['YZZ'] = list_of_octapole_contribs_2[3]

            octapole['YYZ'] = list_of_octapole_contribs_3[0]
            octapole['XYZ'] = list_of_octapole_contribs_3[1]

            octapole['units'] = "Debye-Ang**2"
            return octapole
#            return float(re.findall("(-?\d+(?:\.\d+)?)", prev_line)[3])
        prev_line_3 = prev_line_2
        prev_line_2 = prev_line_1
        prev_line_1 = line

    raise AssertionError("' Octapole moment ' not found in Gaussian calc output")

def hexadecapole(lines):
    hexadecapole = dict()
    prev_line_1 = None
    prev_line_2 = None
    prev_line_3 = None
    prev_line_4 = None
    for line in reversed(lines):
        if "Hexadecapole moment (field-independent basis, Debye-Ang**3)" in line:
#            print('prev line is')
#            print(prev_line)
            list_of_hexadecapole_contribs_1 = re.findall("(-?\d+(?:\.\d+)?)", prev_line_1) #  XXXX YYYY ZZZZ XXXY
            list_of_hexadecapole_contribs_2 = re.findall("(-?\d+(?:\.\d+)?)", prev_line_2) # XXXZ YYYX YYYZ ZZZX
            list_of_hexadecapole_contribs_3 = re.findall("(-?\d+(?:\.\d+)?)", prev_line_3) #  ZZZY XXYY XXZZ YYZZ
            list_of_hexadecapole_contribs_4 = re.findall("(-?\d+(?:\.\d+)?)", prev_line_4) # XXYZ YYXZ ZZXY

#            print(prev_line_1)
#            print(prev_line_2)
#            print(prev_line_3)

            hexadecapole['XXXX'] = list_of_hexadecapole_contribs_1[0]
            hexadecapole['YYYY'] = list_of_hexadecapole_contribs_1[1]
            hexadecapole['ZZZZ'] = list_of_hexadecapole_contribs_1[2]
            hexadecapole['XXXY'] = list_of_hexadecapole_contribs_1[3]

            hexadecapole['XXXZ'] = list_of_hexadecapole_contribs_2[0]
            hexadecapole['YYYX'] = list_of_hexadecapole_contribs_2[1]
            hexadecapole['YYYZ'] = list_of_hexadecapole_contribs_2[2]
            hexadecapole['ZZZX'] = list_of_hexadecapole_contribs_2[3]
            
            hexadecapole['ZZZY'] = list_of_hexadecapole_contribs_3[0]
            hexadecapole['XXYY'] = list_of_hexadecapole_contribs_3[1]
            hexadecapole['XXZZ'] = list_of_hexadecapole_contribs_3[2]
            hexadecapole['YYZZ'] = list_of_hexadecapole_contribs_3[3]

            hexadecapole['XXYZ'] = list_of_hexadecapole_contribs_4[0]
            hexadecapole['YYXZ'] = list_of_hexadecapole_contribs_4[1]
            hexadecapole['ZZXY'] = list_of_hexadecapole_contribs_4[2]

            hexadecapole['units'] = "Debye-Ang**3"
            return hexadecapole
#            return float(re.findall("(-?\d+(?:\.\d+)?)", prev_line)[3])
        prev_line_4 = prev_line_3
        prev_line_3 = prev_line_2
        prev_line_2 = prev_line_1
        prev_line_1 = line

    raise AssertionError("' Hexadecapole moment ' not found in Gaussian calc output")



