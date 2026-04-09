#!/usr/bin/env bash
set -euo pipefail

for Estar in $(seq 5 26); do
    EstarFmt=$(LC_NUMERIC=C printf "%.1f" "$Estar");
    echo "Running GEF_Reader for E* = ${EstarFmt} MeV";
    root -q -b "./UtilityScripts/GEF_Reader.cpp+(92,239,${EstarFmt},1)";
done