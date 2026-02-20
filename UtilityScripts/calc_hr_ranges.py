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

# -------------------------------------------------
# Read reac_info.txt config file 
# -------------------------------------------------

config = configparser.ConfigParser(inline_comment_prefixes=(';',))
config.read('reac_info.txt')

# Read separation energies
Sn_CN = float(config['separation_energies']['Sn_CN'])
Sn_1nDght = float(config['separation_energies'].get('Sn_1nDght', '0'))
Sn_2nDght = float(config['separation_energies'].get('Sn_2nDght', '0'))
Sn_3nDght = float(config['separation_energies'].get('Sn_3nDght', '0'))
Sn_4nDght = float(config['separation_energies'].get('Sn_4nDght', '0'))

# Read excitation energy list
exc_ens_str = config['recoil_info']['recoil_excEns']
exc_ens = [float(x.strip()) for x in exc_ens_str.split(',')]

# Read maximum excitation energy limit
max_excEn = float(config['recoil_info']['max_excEn'])

recoil_Z = int(float(config["recoil_info"]["recoil_Z"]))
recoil_A = int(float(config["recoil_info"]["recoil_A"]))

# -------------------------------------------------
# Sub-functions / Helpers
# -------------------------------------------------

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

# -------------------------------------------------
# HR neutron/gamma channels
# -------------------------------------------------

# Calculate energy thresholds for each channel (3 MeV beyond next channel opening)
# HRg: from min to (Sn_CN + 3)
# HR1n: from Sn_CN to (Sn_CN + Sn_1nDght + 3)
# HR2n: from (Sn_CN + Sn_1nDght) to (Sn_CN + Sn_1nDght + Sn_2nDght + 3)
# HR3n: from (Sn_CN + Sn_1nDght + Sn_2nDght) to (Sn_CN + Sn_1nDght + Sn_2nDght + Sn_3nDght + 3)
# HR4n: from (Sn_CN + Sn_1nDght + Sn_2nDght + Sn_3nDght) to (Sn_CN + Sn_1nDght + Sn_2nDght + Sn_3nDght + Sn_4nDght + 3)

HRg_start_E = min(exc_ens)
HRg_stop_E = Sn_CN + 3.0

HR1n_start_E = Sn_CN
HR1n_stop_E = Sn_CN + Sn_1nDght + 3.0

HR2n_start_E = Sn_CN + Sn_1nDght
HR2n_stop_E = Sn_CN + Sn_1nDght + Sn_2nDght + 3.0

HR3n_start_E = Sn_CN + Sn_1nDght + Sn_2nDght
HR3n_stop_E = Sn_CN + Sn_1nDght + Sn_2nDght + Sn_3nDght + 3.0

HR4n_start_E = Sn_CN + Sn_1nDght + Sn_2nDght + Sn_3nDght
HR4n_stop_E = Sn_CN + Sn_1nDght + Sn_2nDght + Sn_3nDght + Sn_4nDght + 3.0

HRg_start = find_start_index(HRg_start_E)
HRg_stop = find_stop_index(HRg_stop_E)

HR1n_start = find_start_index(HR1n_start_E)
HR1n_stop = find_stop_index(HR1n_stop_E)

HR2n_start = find_start_index(HR2n_start_E)
HR2n_stop = find_stop_index(HR2n_stop_E)

HR3n_start = find_start_index(HR3n_start_E)
HR3n_stop = find_stop_index(HR3n_stop_E)

HR4n_start = find_start_index(HR4n_start_E)
HR4n_stop = find_stop_index(HR4n_stop_E)

# Find the maximum allowed index based on max_excEn
max_index = find_stop_index(max_excEn)

# Cap all stop indices to not exceed max_index
HRg_stop = min(HRg_stop, max_index)
HR1n_stop = min(HR1n_stop, max_index)
HR2n_stop = min(HR2n_stop, max_index)
HR3n_stop = min(HR3n_stop, max_index)
HR4n_stop = min(HR4n_stop, max_index)

# -------------------------------------------------
# HRf (fission window)
# -------------------------------------------------
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


# -------------------------------------------------
# Output for Bash (eval)
# -------------------------------------------------

print("HRg_start={}".format(HRg_start))
print("HRg_stop={}".format(HRg_stop))
print("HR1n_start={}".format(HR1n_start))
print("HR1n_stop={}".format(HR1n_stop))
print("HR2n_start={}".format(HR2n_start))
print("HR2n_stop={}".format(HR2n_stop))
print("HR3n_start={}".format(HR3n_start))
print("HR3n_stop={}".format(HR3n_stop))
print("HR4n_start={}".format(HR4n_start))
print("HR4n_stop={}".format(HR4n_stop)) 
print("HRf_start={}".format(HRf_start))
print("HRf_stop={}".format(HRf_stop))
print("Ua={}".format(Ua))

