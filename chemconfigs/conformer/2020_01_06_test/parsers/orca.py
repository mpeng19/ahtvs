import itertools
import re


def parse_program(lines):
    for line in lines:
            if "Program Version" in line:
                text = line.split()
                program = "ORCA (" + " ".join(text[:3]) + ")"
#                version = text[5]
                return program #, version


def assert_completed(lines):
    for line in reversed(lines):
        if '****ORCA TERMINATED NORMALLY****' in line:
            return
    raise AssertionError("****ORCA TERMINATED NORMALLY****' line not found in Orca calc output")

def energy(lines):
    for line in reversed(lines):
        if 'FINAL SINGLE POINT ENERGY' in line:
            return float(re.findall("(-?\d+(?:\.\d+)?)", line)[0])

    raise AssertionError("'Energy is ' not found in Orca calc output")


