#!/usr/bin/env python3.5
# -*- coding: utf-8 -*-
"""
Calculate HR channel excitation energy ranges based on separation energies.
Outputs bash variable assignments that can be sourced.
"""

import configparser
import math
import os
import sys

# Read config file
config = configparser.ConfigParser(inline_comment_prefixes=(';',))
config.read('reac_info.txt')

# Read separation energies
Sn_CN = float(config['separation_energies']['Sn_CN'])
Sn_1nDght = float(config['separation_energies'].get('Sn_1nDght', '0'))
Sn_2nDght = float(config['separation_energies'].get('Sn_2nDght', '0'))
Sn_3nDght = float(config['separation_energies'].get('Sn_3nDght', '0'))
Sn_4nDght = float(config['separation_energies'].get('Sn_4nDght', '0'))

recoil_Z = int(float(config["recoil_info"]["recoil_Z"]))
recoil_A = int(float(config["recoil_info"]["recoil_A"]))

# Read excitation energy parameters (required)
if 'excEn_start' not in config['recoil_info'] or 'excEn_stop' not in config['recoil_info'] or 'excEn_bin' not in config['recoil_info']:
    print("Error: Missing required parameters excEn_start, excEn_stop, or excEn_bin in [recoil_info] section", file=sys.stderr)
    sys.exit(1)

excEn_start = float(config['recoil_info']['excEn_start'])
excEn_stop = float(config['recoil_info']['excEn_stop'])
excEn_bin = float(config['recoil_info']['excEn_bin'])

# Read run_HR_modes to determine which channels to calculate
if 'run_HR_modes' not in config['recoil_info']:
    print("Error: Missing required parameter run_HR_modes in [recoil_info] section", file=sys.stderr)
    sys.exit(1)

run_HR_modes_str = config['recoil_info']['run_HR_modes']
run_HR_modes = [mode.strip() for mode in run_HR_modes_str.split(',')]

# Generate excitation energy list: start, start+bin, start+2*bin, ... up to and including stop
# Use a small epsilon to handle floating point precision when including stop
exc_ens = []
current = excEn_start
while current <= excEn_stop + 1e-6:  # Add small epsilon to handle floating point precision
    exc_ens.append(current)
    current += excEn_bin

# Read overlap range parameter (3.0 MeV)
overlap_range = 3.0

# Find indices in exc_ens array that correspond to these energies
# Find the first index where exc_ens[i] >= threshold
def find_start_index(threshold):
    for i, E in enumerate(exc_ens):
        if E >= threshold:
            return i
    return len(exc_ens) - 1  # If threshold is beyond all energies, use last index

# Find the first index where exc_ens[i] >= threshold (inclusive of next value)
def find_stop_index(threshold):
    for i, E in enumerate(exc_ens):
        if E >= threshold:
            return i
    return len(exc_ens) - 1  # If threshold is beyond all energies, use last index

def load_bjornholm_uA(path):
    """
    Get Bjornholm Barrier map:
    columns = Z  A  Ua(MeV)  Ub(MeV)
    Returns dict[(Z,A)] -> Ua (float or None)
    """
    table = {}
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue

            parts = s.split()
            if len(parts) < 4:
                continue

            Z = int(parts[0])
            A = int(parts[1])

            try:
                ua = float(parts[2])
                if math.isnan(ua):
                    ua = None
            except ValueError:
                ua = None

            table[(Z, A)] = ua

    return table


# Calculate energy thresholds and indices only for modes in run_HR_modes
HR_ranges = {}

if 'HRg' in run_HR_modes:
    HRg_start_E = min(exc_ens)
    HRg_stop_E = Sn_CN + overlap_range
    HR_ranges['HRg'] = (find_start_index(HRg_start_E), find_stop_index(HRg_stop_E))

if 'HR1n' in run_HR_modes:
    HR1n_start_E = Sn_CN
    HR1n_stop_E = Sn_CN + Sn_1nDght + overlap_range
    HR_ranges['HR1n'] = (find_start_index(HR1n_start_E), find_stop_index(HR1n_stop_E))

if 'HR2n' in run_HR_modes:
    HR2n_start_E = Sn_CN + Sn_1nDght
    HR2n_stop_E = Sn_CN + Sn_1nDght + Sn_2nDght + overlap_range
    HR_ranges['HR2n'] = (find_start_index(HR2n_start_E), find_stop_index(HR2n_stop_E))

if 'HR3n' in run_HR_modes:
    HR3n_start_E = Sn_CN + Sn_1nDght + Sn_2nDght
    HR3n_stop_E = Sn_CN + Sn_1nDght + Sn_2nDght + Sn_3nDght + overlap_range
    HR_ranges['HR3n'] = (find_start_index(HR3n_start_E), find_stop_index(HR3n_stop_E))

if 'HR4n' in run_HR_modes:
    HR4n_start_E = Sn_CN + Sn_1nDght + Sn_2nDght + Sn_3nDght
    HR4n_stop_E = Sn_CN + Sn_1nDght + Sn_2nDght + Sn_3nDght + Sn_4nDght + overlap_range
    HR_ranges['HR4n'] = (find_start_index(HR4n_start_E), find_stop_index(HR4n_stop_E))

# HRf (fission window)
barrier_file = os.environ.get(
    "BJORNHOLM_BARRIER_FILE",
    "Bjornholm_Fission_Barrier.txt"
)

if not os.path.exists(barrier_file):
    sys.stderr.write("ERROR: barrier file not found: {}\n".format(barrier_file))
    sys.exit(1)

barriers = load_bjornholm_uA(barrier_file)

Ua = barriers.get((recoil_Z, recoil_A), None)

if Ua is None:
    sys.stderr.write("ERROR: Ua not found for Z={} A={}\n".format(recoil_Z, recoil_A))
    sys.exit(2)

HRf_start_E = Ua - 2.0
HRf_stop_E  = max_excEn

if HRf_start_E < min(exc_ens):
    HRf_start_E = min(exc_ens)

HRf_start = find_start_index(HRf_start_E)
HRf_stop  = find_stop_index(HRf_stop_E)

HRf_stop = min(HRf_stop, max_index)

if HRf_start > HRf_stop:
    HRf_start = max_index + 1
    HRf_stop  = max_index

print("excitation_Ens=(%s)" % exc_ens_str)

# Output run_HR_modes array
run_HR_modes_str = ' '.join(run_HR_modes)
print("run_HR_modes=(%s)" % run_HR_modes_str)

# Output HR mode ranges only for modes in run_HR_modes
for mode in run_HR_modes:
    if mode in HR_ranges:
        start_idx, stop_idx = HR_ranges[mode]
        print("%s_start=%d" % (mode, start_idx))
        print("%s_stop=%d" % (mode, stop_idx))
        
print("HRf_start={}".format(HRf_start))
print("HRf_stop={}".format(HRf_stop))
print("Ua={}".format(Ua))