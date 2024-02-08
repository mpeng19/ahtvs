import re


def duration(lines):
    for line in reversed(lines):
        if "PWSCF        :" in line:
            # PW wall clock timing is printed in a variable format:
            #   8m21.31s WALL
            #   1h 4m WALL
            line = line.split('WALL')[0]
            time = line.split('CPU')[1]
            # Convert WALL time from QE's h:m:s form into seconds of run time
            time = re.sub('m', '*60+', time)
            time = re.sub('s', '+', time)
            time = re.sub('h', '*60*60+', time)
            time = time+'0'
            time_min = eval(time)
            return time_min


def nprocs(lines):
    for line in lines:
        if "Number of MPI processes" in line:
            return int(re.findall("(\d+)", line)[0])


def mbd_energy(lines):
    for line in reversed(lines):
        if "Dispersion MBD Correction" in line:
            return float(re.findall("(-?\d+(?:\.\d+)?)", line)[0])
    raise AssertionError("'Dispersion MBD Correction' not found in QE calc output")


def total_energy(lines):
    for line in reversed(lines):
        if "!    total energy" in line:
            return float(re.findall("(-?\d+(?:\.\d+)?)", line)[0])
    raise AssertionError("'Total energy is ' not found in QE calc output")
