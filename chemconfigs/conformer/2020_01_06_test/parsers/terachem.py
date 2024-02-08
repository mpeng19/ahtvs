import re


def duration(lines):
    for line in lines:
        if "Total processing time" in line:
            text = line.split()
            return float(text[3])


def host(lines):
    for line in lines:
        if "On" in line:
            text = line.split()
            return text[1]


def energy(lines):
    for line in reversed(lines):
        if "FINAL ENERGY:" in line:
            return float(re.findall("(-?\d+(?:\.\d+)?)", line)[0])
    raise AssertionError("'FINAL ENERGY' not found in TeraChem calc output")


def dipole(lines):
    for line in reversed(lines):
        if "DIPOLE MOMENT: " in line:
            return float(re.findall("(-?\d+(?:\.\d+)?)", line)[-1])
    raise AssertionError("'DIPOLE MOMENT: ' not found in TeraChem calc output")


def orbitals(molden_content):
    """
    build this field: orbital_energy_list = ListField(DictField())
    # each dict contains keys: 'method', 'spin', 'type', 'number', 'energy'
    """
    orbitals = re.findall("Ene=\s*(\-?\d+(?:\.\d+)?)\s*Spin=\s*(\w*)\s*Occup=\s*(\-?\d+(?:\.\d+)?)", molden_content)
    occ_num = 0
    unocc_num = 0
    orb_list = []
    homo = None
    lumo = None
    for orb in orbitals:
        spin = orb[1].lower()
        energy = float(orb[0])
        if float(orb[2]) == 0:
            unocc_num += 1
            number = unocc_num
            occtype = 'unocc'
            if lumo is None:
                lumo = energy
        else:
            occ_num += 1
            number = occ_num
            occtype = 'occ'
            homo = energy  # overwritten every time, left at final
        d = dict(method="rhf", spin=spin, type=occtype, number=number, energy=energy)
        orb_list.append(d)
    return orb_list, homo, lumo
