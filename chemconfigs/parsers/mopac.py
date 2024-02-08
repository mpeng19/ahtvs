import re


def duration(lines):
    for line in lines:
        if "COMPUTATION TIME" in line:
            text = line.strip().split()
            return float(text[3])


def version(lines):
    for line in lines:
        if "MOPAC2012" in line:
            text = line.strip().split()
            return text[2]


def dipole(lines):
    for line in lines:
        if "DIPOLE" in line:
            text = line.strip().split()
            return float(text[2])


def molwt(lines):
    for line in lines:
        if "MOLECULAR WEIGHT" in line:
            text = line.strip().split()
    return float(text[3])


def energies(lines):
    for line in lines:
        if "TOTAL ENERGY" in line:
            text = line.strip().split()
            totale = float(text[3])
        if "ELECTRONIC ENERGY" in line:
            text = line.strip().split()
            elece = float(text[3])
        if "CORE-CORE" in line:
            text = line.strip().split()
            nuce = float(text[3])
        if "HEAT OF FORMATION " in line:
            text = line.strip().split()
            heatof = float(text[4])
    return (totale, elece, nuce, heatof)


def charge(lines):
    for line in lines:
        if "CHARGE=" in line:
            charge = re.findall("\-?\d+(?:\.\d+)?", line.split('CHARGE=')[1])[0]
            return charge


def geometry(arcfile):
    lines = arcfile.split('\n')
    lines = [x.strip() for x in lines] 
    for i,line in enumerate(lines):
        if 'FINAL GEOMETRY OBTAINED' in line:
            coord_start = i

    cur_line = coord_start+4
    coords = []
    line = lines[cur_line]
    done = False
    while not done:
        parsed_coords = line.split()
        atom_char = parsed_coords[0]
        x = parsed_coords[1]
        y = parsed_coords[3]
        z = parsed_coords[5]
        coords.append(dict(element=atom_char, x=x, y=y, z=z))
        cur_line += 1
        line = lines[cur_line]
        if not line.strip():
            done = True
#    coords = re.findall("(\w*)\s*(\-?\d+(?:\.\d+)?)\s*\+1\s*(\-?\d+(?:\.\d+)?)\s*\+1\s*(\-?\d+(?:\.\d+)?)", arcfile)
#    print(len(coords))
#    for i in coords:
#        coords.append(dict(element=i[0], x=float(i[1]), y=float(i[2]), z=float(i[3])))
    return coords


def natoms(lines):
    for line in lines:
        if "Empirical Formula" in line:
            number = int(re.findall("=\s*(\-?\d+(?:\.\d+)?)", line)[0])
    return number


def fmos(lines):
    for line in lines:
        if "HOMO LUMO ENERGIES" in line:
            homo, lumo = re.findall("\-?\d+(?:\.\d+)?", line)
            return (float(homo), float(lumo))


def rhf_eigenvalues(outfile):
    eigens = []
    with open(outfile, 'r') as infile:
        copy = False
        for line in infile:
            if "FILLED" in line.strip():
                n_occ_orbs = int(line.strip().split()[5])
            if "EIGENVALUES" in line.strip():
                copy = True
            elif "NET ATOMIC CHARGES AND DIPOLE CONTRIBUTIONS" in line.strip():
                copy = False
            elif copy:
                eigens.append(line.strip().split())
    flattedeigens = [float(eigens[iter1][iter2]) for iter1 in range(len(eigens)) for iter2 in range(len(eigens[iter1]))]
    # print flattedeigens[n_occ_orbs]
    occ = [dict(method="rhf", spin="alpha", type="occ", number=num+1, energy=i) for num, i in enumerate(flattedeigens[0:n_occ_orbs])]
    unocc = [dict(method="rhf", spin="alpha", type="unocc", number=num+1, energy=i) for num, i in enumerate(flattedeigens[n_occ_orbs:])]
    return (n_occ_orbs, occ + unocc)


def vibrations(loglines, atom_number):
    """"
    build this field normal_mode_list = ListField(DictField())
    # each dict contains keys: 'number', 'freq' in cm-1, 'dipole',
    'travel', 'red_mass', 'eff_mass', 'coords'
     which is a list of floats in mass-weighted coordinates
    """
    normal_modes = []
    # number_modes = 3 * atom_number - 6

    start_line = loglines.index('          DESCRIPTION OF VIBRATIONS\n')
    end_line = loglines.index('           FORCE CONSTANT IN CARTESIAN COORDINATES (Millidynes/A)\n')

    m = old_m = 1
    for line in loglines[start_line:end_line+1]:
        line = line.strip()

        if line.startswith('VIBRATION'):
            m = int(re.findall("VIBRATION\s*(\-?\d+(?:\.\d+)?)", line)[0])
        if line.startswith('FREQ.'):
            f = float(re.findall("FREQ.\s*(\-?\d+(?:\.\d+)?)", line)[0])
        if line.startswith('T-DIPOLE'):
            d = float(re.findall("T-DIPOLE\s*(\-?\d+(?:\.\d+)?)", line)[0])
        if line.startswith('TRAVEL'):
            t = float(re.findall("TRAVEL\s*(\-?\d+(?:\.\d+)?)", line)[0])
        if line.startswith('RED. MASS'):
            r = float(re.findall("RED. MASS\s*(\-?\d+(?:\.\d+)?)", line)[0])
        if line.startswith('EFF. MASS'):
            e = float(re.findall("EFF. MASS\s*(\-?\d+(?:\.\d+)?)", line)[0])

        if m != old_m or 'FORCE' in line:
            item = dict(number=m, freq=f, dipole=d, travel=t, red_mass=r, eff_mass=e)
            normal_modes.append(item)
            old_m = m

    start_line = loglines.index('           MASS-WEIGHTED COORDINATE ANALYSIS (NORMAL COORDINATES)\n')
    end_line = loglines.index('          DESCRIPTION OF VIBRATIONS\n')

    coord_motions = [[] for i in range(3*atom_number)]

    for line in loglines[start_line:end_line]:
        for coord in range(3*atom_number):
            if line.strip().startswith(str(coord+1)+'  '):
                coord_motions[coord] = coord_motions[coord] + line.split()[1:]

    motions_per_mode = zip(*coord_motions)

    for i, coords in enumerate(motions_per_mode):
        normal_modes[i]['coords'] = coords
    return normal_modes
