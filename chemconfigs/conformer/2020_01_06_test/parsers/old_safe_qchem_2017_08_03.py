import itertools
import re


TDA_MARKER = "TDDFT/TDA Excitation Energies"
RPA_MARKER = "TDDFT Excitation Energies"
TDA_POLES_MARKER = "TDA Excited-State Multipoles, State"
RPA_POLES_MARKER = "RPA Excited-State Multipoles, State"
GS_POLES_MARKER = "Cartesian Multipole Moments"

ENDING_MARKER = "------------------------"

ORB_MARKER = "Occupied"
ORB_INT_MARKER = "Virtual"
ORB_ENDING_MARKER = "-" * 62
COORD_ENDING_MARKER = "Point Group"
COORD_MARKER = "Coordinates (Angstroms)"


def parse_program_version(lines):
    for line in lines:
            if "Intel X86" in line:
                text = line.split()
                program = "Q-Chem (" + " ".join(text[:4]) + ")"
                version = text[5]
                return program, version


def assert_completed(lines):
    for line in reversed(lines):
        if 'Thank you very much for using Q-Chem.  Have a nice day.' in line:
            return
    raise AssertionError("Have a nice day.' line not found in Q-Chem calc output")


def parse_multipoles(lines):
    cur_line = 2
    line = lines[cur_line]
    while ENDING_MARKER not in line:
        if 'Dipole Moment (Debye)' in line:
            dipole_dict = dict()
            dipole_line = lines[cur_line+1].split()
            for i, j in enumerate(dipole_line[::2]):
                dipole_dict.update({str(dipole_line[2*i]): float(dipole_line[2*i+1])})
        elif 'Quadrupole Moments (Debye-Ang)' in line:
            quad_dict = dict()
            multipole_lines = [lines[cur_line+i+1].split() for i in range(2)]
            quad_lines = list(itertools.chain(*multipole_lines))
            for i, j in enumerate(quad_lines[::2]):
                quad_dict.update({str(quad_lines[2*i]): float(quad_lines[2*i+1])})
        elif 'Octopole Moments (Debye-Ang^2)' in line:
            octa_dict = dict()
            multipole_lines = [lines[cur_line+i+1].split() for i in range(4)]
            octa_lines = list(itertools.chain(*multipole_lines))
            for i, j in enumerate(octa_lines[::2]):
                octa_dict.update({str(octa_lines[2*i]): float(octa_lines[2*i+1])})
        elif 'Hexadecapole Moments (Debye-Ang^3)' in line:
            hexadeca_dict = dict()
            multipole_lines = [lines[cur_line+i+1].split() for i in range(5)]
            hexadeca_lines = list(itertools.chain(*multipole_lines))
            for i, j in enumerate(hexadeca_lines[::2]):
                hexadeca_dict.update({str(hexadeca_lines[2*i]): float(hexadeca_lines[2*i+1])})
        cur_line += 1
        line = lines[cur_line]

    return [dipole_dict, quad_dict, octa_dict, hexadeca_dict]


def duration(lines):
    for line in reversed(lines):
        if "Total job time" in line:
            return float(re.findall("(\d+(?:\.\d+)?)", line)[0])


def host(lines):
    for line in lines:
        if "Host:" in line:
            text = line.split()
            if len(text) > 1:
                return text[1]


def nprocs(lines):
    for line in lines:
        if "Parallel job on" in line:
            return int(re.findall("(\d+(?:\.\d+)?)", line)[0])


def energy(lines):
    for line in reversed(lines):
        if 'Final energy is' in line:
            return float(re.findall("(-?\d+(?:\.\d+)?)", line)[0])

    for line in reversed(lines):
        if "Total energy in the final basis set" in line:
            return float(re.findall("(-?\d+(?:\.\d+)?)", line)[0])

    energy_terms = ['Total Coulomb',
                    'Alpha Exchange',
                    'Beta  Exchange',
                    'DFT   Exchange',
                    'DFT Correlation',
                    'Nuclear Repu',
                    'Nuclear Attr',
                    'Kinetic']
    vals = []
    for term in energy_terms:
        for line in reversed(lines):
            if term in line:
                vals.append(float(re.findall("(-?\d+(?:\.\d+)?)", line)[0]))
                break
    return sum(vals)
    raise AssertionError("'Energy is ' not found in QChem calc output")


def dipole(lines):
    for line in reversed(lines):
        if " Tot " in line:
            return float(re.findall("(-?\d+(?:\.\d+)?)", line)[0])
    raise AssertionError("' Tot ' not found in QChem calc output")


def parse_coords_starting_at_line(lines, start_line_num):
    cur_line = start_line_num
    coords = []
    line = lines[cur_line]
    while COORD_ENDING_MARKER not in line:
        atom_char = line.split()[1]
        num, x, y, z = re.findall("\-?\d+(?:\.\d+)?", line)
        coords.append(dict(element=atom_char, x=x, y=y, z=z))
        cur_line += 1
        line = lines[cur_line]
    return coords


def last_coordinate_list(lines):
    for j, line in enumerate(lines):
        if COORD_MARKER in line:
            normal_start = j
    return parse_coords_starting_at_line(lines, normal_start+2)


def parse_states_from_line(lines, start_line_num, expect_state_count):
    cur_line = start_line_num
    state_list = []
    current_state = None
    line = lines[cur_line]
    while ENDING_MARKER not in line:
        if "Excited state" in line:
            if current_state is not None:
                state_list.append(current_state)
            current_state = {}
            current_state["orbital_list"] = []
            num, energy = re.findall("(\d+(?:\.\d+)?)", line)
            current_state["state_number"] = int(num)
            current_state["energy"] = float(energy)

        elif "Multiplicity" in line:
            current_state["multiplicity"] = line.split()[1].lower()
        elif "Strength" in line:
            strength = re.findall("(\d+(?:\.\d+)?)", line)[0]
            current_state["oscillator_strength"] = float(strength)
        elif 'PBHT overlap' in line:
            pbht = re.findall("(\d+(?:\.\d+)?)", line)[0]
            current_state["pbht"] = float(pbht)
        elif 'Trans. Mom.' in line:
            trans_dipole = re.findall("(\d+(?:\.\d+)?)", line)[:3]
            current_state["trans_dipole"] = [float(i) for i in trans_dipole]
        elif "-->" in line:
            i, f, a = re.findall("\-?\d+(?:\.\d+)?", line)
            initial = int(i)
            final = int(f)
            amp = float(a)
            current_state["orbital_list"].append(dict(initial=initial, final=final, amplitude=amp))
        cur_line += 1
        line = lines[cur_line]
    state_list.append(current_state)  # put on last excited state
    assert(len(state_list) == expect_state_count)
    return state_list


# def tddft_excitations(lines):
#     for j, line in enumerate(lines):
#         if 'cis_n_roots' in line:
#             total_states = int(line.split()[1])
#         if TDA_MARKER in line:
#             tda_start = j + 2
#     tda_excited_states = parse_states_from_line(lines, tda_start, total_states)
#     return tda_excited_states


def tddft_excitations(lines, include_rpa=True, triplets=True):
    for j, line in enumerate(lines):
        if 'cis_n_roots' in line:
            total_states = int(line.split()[1])
        if TDA_MARKER in line:
            tda_start = j + 2
        if include_rpa and RPA_MARKER in line:
            rpa_start = j + 2
    if triplets:
        total_states *= 2

    tda_excited_states = parse_states_from_line(lines, tda_start, total_states)
    if include_rpa:
        rpa_excited_states = parse_states_from_line(lines, rpa_start, total_states)
        return tda_excited_states, rpa_excited_states
    else:
        return tda_excited_states

def update_tddft_excitations_w_multipole(lines, tda_excited_states, rpa_excited_states=None):
    for j, line in enumerate(lines):
        if TDA_POLES_MARKER in line:
            tda_start = j
            break

    for j, line in enumerate(lines[tda_start:]):
        if 'TDA Excited-State Multipoles, State' in line:
            n_state = int(re.findall("(\d+(?:\.\d+)?)", line)[0]) - 1
            multipoles = parse_multipoles(lines[tda_start+j:tda_start+j+25])
            tda_excited_states[n_state].update({'multipoles': multipoles})

    if rpa_excited_states:
        for j, line in enumerate(lines):
            if RPA_POLES_MARKER in line:
                rpa_start = j
                break

        for j, line in enumerate(lines[rpa_start:]):
            if 'RPA Excited-State Multipoles, State' in line:
                n_state = int(re.findall("(\d+(?:\.\d+)?)", line)[0]) - 1
                multipoles = parse_multipoles(lines[rpa_start+j:rpa_start+j+25])
                rpa_excited_states[n_state].update({'multipoles': multipoles})

        return tda_excited_states, rpa_excited_states
    else:
        return tda_excited_states


def ground_state_multipoles(lines):
    for j, line in enumerate(lines):
        if GS_POLES_MARKER in line:
            multipoles = parse_multipoles(lines[j:j+25])
            return multipoles


def get_lines_for_interstate_dipoles(lines):
    state_to_state_string = '                    STATE-TO-STATE TRANSITION MOMENTS\n'
    start_line = lines.index(state_to_state_string)

    start_lines_dict = {'rpa': [], 'tda': []}

    for j, line in enumerate(lines[start_line:]):
        if 'Within CIS/TDA Excited States:' in line:
            method = 'tda'
        elif 'Within RPA/TDDFT Excited States:' in line:
            method = 'rpa'
        elif 'Transition Moments Between Triplet Excited States' in line:
            start_lines_dict[method].append(j+4+start_line)
        elif 'Transition Moments Between Ground and Singlet Excited States' in line:
            start_lines_dict[method].append(j+4+start_line)
        elif 'Transition Moments Between Singlet Excited States' in line:
            start_lines_dict[method].append(j+4+start_line)

    return start_lines_dict


def parse_interstate_strengths(lines, dictionary):
    cur_line = 0
    line = lines[cur_line]
    while ENDING_MARKER not in line:
        initial, final, x, y, z, strength = re.findall("[+-]?(?=\d*)(?=\.?\d)\d*\.?\d*(?:[eE][+-]?\d+)?", line)
        initial, final = int(initial), int(final)
        if initial == 0:
            initial = final
            final = 0
        if initial in dictionary:
            dictionary[int(initial)][int(final)] = float(strength)
        else:
            dictionary[int(initial)] = {int(final): float(strength)}
        cur_line += 1
        line = lines[cur_line]
    return dictionary


def update_tddft_dict_w_intestate_strengths(lines, startlines, excited_states_dict):
    additional_states_dipoles = {0: dict()}
    for start in startlines:
        additional_states_dipoles = parse_interstate_strengths(lines[start:], additional_states_dipoles)
    for i in excited_states_dict:
        if i['state_number'] in additional_states_dipoles:
            i['interstate_strengths'] = additional_states_dipoles[i['state_number']]
    return excited_states_dict


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


