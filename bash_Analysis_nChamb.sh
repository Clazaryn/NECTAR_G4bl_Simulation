#!/bin/bash

reaction=$(grep '^reaction' reac_info_nChamb.txt | awk -F'=' '{gsub(/^ +| +$/,"",$2); print $2}' | awk '{print $1}')

echo "   ####################################################   "
echo "                                                          "
echo "     SIMULATIONS FOR NEW NECTAR REACTION CHAMBER          "
echo "     Analysis of Detector output for  HR modes            "
echo "     Reaction: $reaction - check scripts!                 "
echo "                                                          "
echo "   ####################################################   "

# Read excitation energies from reac_info_nChamb.txt
excitation_Ens=($(grep '^recoil_excEns' reac_info_nChamb.txt | awk -F'=' '{print $2}'))
echo "Excitation energies: ${excitation_Ens[@]}"

# 206Pb(d,p): g = 0-9 MeV, 1n = 6.7-17 MeV, 2n = 15-23 MeV, 3n = 21.5-30 MeV
# 206Pb(p,p')/(d,d'): g = 0-10 MeV, 1n = 8-17 MeV, 2n = 15-20 MeV
# \!/\!/\!/ ------- DISABLE HR3n for (p,p') or (d,d'), otherwise the job hangs!!! ------- \!/\!/\!/
if [ "$reaction" == "206Pbdp" ]; then
  HRg_start=0;    HRg_stop=9;
  HR1n_start=7;   HR1n_stop=17;
  HR2n_start=15;  HR2n_stop=23;
  HR3n_start=22;   HR3n_stop=26;
elif [[ "$reaction" == "206Pbdd" || "$reaction" == "206Pbpp" ]]; then
  HRg_start=0;    HRg_stop=10;
  HR1n_start=8;   HR1n_stop=17;
  HR2n_start=15;  HR2n_stop=20;
fi

# Loop over excitation energies and run all applicable HR modes sequentially
for i in $(seq 0 $((${#excitation_Ens[@]}-1))); do

    en=${excitation_Ens[$i]}
    lbl=$(printf "%02d" $(echo "$en" | awk '{print int($1 + 0.5)}'))  # Round to nearest integer

    if [ "$i" -ge $HRg_start ] && [ "$i" -le $HRg_stop ]; then
        echo "Running analysis_nChamb.C($i, \"HRg\")"
        root -l -q "analysis_nChamb.C($i, \"HRg\")"
    fi

    if [ "$i" -ge $HR1n_start ] && [ "$i" -le $HR1n_stop ]; then
        echo "Running analysis_nChamb.C($i, \"HR1n\")"
        root -l -q "analysis_nChamb.C($i, \"HR1n\")"
    fi

    if [ "$i" -ge $HR2n_start ] && [ "$i" -le $HR2n_stop ]; then
        echo "Running analysis_nChamb.C($i, \"HR2n\")"
        root -l -q "analysis_nChamb.C($i, \"HR2n\")"
    fi

    if [ "$reaction" == "206Pbdp" ]; then
    if [ "$i" -ge $HR3n_start ] && [ "$i" -le $HR3n_stop ]; then
            echo "Running analysis_nChamb.C($i, \"HR3n\")"
            root -l -q "analysis_nChamb.C($i, \"HR3n\")"
        fi
    fi

done
