import json
import re
import pandas as pd

def _parse_xyz_population(s):
    """
    Parse XYZ format and return element numbers and Cartesian coordinates.

    Arguments
    ---------
    s: Input coordinate string

    Returns
    -------
    coords :

    """

    # Create regexes
    re_atom = re.compile(r'^\s*[a-zA-z]+')
    re_num = re.compile(r'[-]?\d+\.\d+')

    # Skip the first two lines
    _, _, atoms = s.split('\n', 2)

    coords = []


    # Iterate over atoms
    for line in atoms.splitlines():
        if not line.strip():
            continue
        m = re_atom.match(line)
        if m:
            atom = m.group(0).strip()
            atx, aty, atz, population = [float(number) for number in re_num.findall(line)]
            coords.append(dict(element=atom, x=atx, y=aty, z=atz))
    return coords


def get_wall_time(lines):
    for line in reversed(lines):
        if "real" in line:
            return float(re.findall("(\d+(?:\.\d+)?)", line)[0])

def get_cpu_time(lines):
    for line in reversed(lines):
        if "user" in line:
            user = float(re.findall("(\d+(?:\.\d+)?)", line)[0])
        if "sys" in line:
            sys = float(re.findall("(\d+(?:\.\d+)?)", line)[0])
    return sys+user

def get_release(lines):
    dftb_release=""
    for line in lines:
        if "Release:" in line:
            dftb_release = line[line.find("Release:") + len("Release:") + 1:].strip()
    return dftb_release

def get_excited_states_dict(file):
    df1 = pd.read_fwf(file, widths=[14, 16, 20, 13, 10, 7], header=1, comment='=', skip_blank_lines=True)
    return json.loads(df1.dropna(axis=0).reset_index(drop=True).to_json())
