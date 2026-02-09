#!/usr/bin/env python3.5
# -*- coding: utf-8 -*-
"""
Calculate HR channel excitation energy ranges based on separation energies.
Outputs bash variable assignments that can be sourced.
"""

import configparser
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

# Read excitation energy list
exc_ens_str = config['recoil_info']['recoil_excEns']
exc_ens = [float(x.strip()) for x in exc_ens_str.split(',')]

# Read maximum excitation energy limit
max_excEn = float(config['recoil_info']['max_excEn'])

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

# Output bash variable assignments
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
