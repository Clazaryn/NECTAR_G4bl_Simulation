#!/bin/bash

reaction=$(grep '^reaction' reac_info.txt | awk -F'=' '{gsub(/^ +| +$/,"",$2); print $2}' | awk '{print $1}')

echo "   ####################################################   "
echo "                                                          "
echo "     SIMULATIONS FOR NEW NECTAR REACTION CHAMBER          "
echo "     Analysis of Detector output for  HR modes            "
echo "     Reaction: $reaction - check scripts!                 "
echo "                                                          "
echo "   ####################################################   "

# Read excitation energies from reac_info.txt
excitation_Ens=($(grep '^recoil_excEns' reac_info.txt | awk -F'=' '{print $2}' | tr ',' ' '))
echo "Excitation energies: ${excitation_Ens[@]}"

# Calculate HR channel ranges based on separation energies (3 MeV beyond next channel opening)
eval $(python3 UtilityScripts/calc_hr_ranges.py)
echo "HR channel ranges:"
echo "  HRg: indices $HRg_start to $HRg_stop (${excitation_Ens[$HRg_start]} to ${excitation_Ens[$HRg_stop]} MeV)"
echo "  HR1n: indices $HR1n_start to $HR1n_stop (${excitation_Ens[$HR1n_start]} to ${excitation_Ens[$HR1n_stop]} MeV)"
echo "  HR2n: indices $HR2n_start to $HR2n_stop (${excitation_Ens[$HR2n_start]} to ${excitation_Ens[$HR2n_stop]} MeV)"
echo "  HR3n: indices $HR3n_start to $HR3n_stop (${excitation_Ens[$HR3n_start]} to ${excitation_Ens[$HR3n_stop]} MeV)"
if [ "$HR4n_start" -le "$HR4n_stop" ] && [ "$HR4n_stop" -ge 0 ]; then
  echo "  HR4n: indices $HR4n_start to $HR4n_stop (${excitation_Ens[$HR4n_start]} to ${excitation_Ens[$HR4n_stop]} MeV)"
fi

# Calculate total number of iterations across all channels
total_iterations=0
total_iterations=$((total_iterations + HRg_stop - HRg_start + 1))
total_iterations=$((total_iterations + HR1n_stop - HR1n_start + 1))
if [ "$HR2n_start" -le "$HR2n_stop" ] && [ "$HR2n_stop" -ge 0 ]; then
  total_iterations=$((total_iterations + HR2n_stop - HR2n_start + 1))
fi
if [ "$HR3n_start" -le "$HR3n_stop" ] && [ "$HR3n_stop" -ge 0 ]; then
  total_iterations=$((total_iterations + HR3n_stop - HR3n_start + 1))
fi
if [ "$HR4n_start" -le "$HR4n_stop" ] && [ "$HR4n_stop" -ge 0 ]; then
  total_iterations=$((total_iterations + HR4n_stop - HR4n_start + 1))
fi

completed_iterations=0  # Track completed iterations

# Function to display progress bar
show_progress() {
  local completed=$1
  local total=$2
  local percent=$((completed * 100 / total))
  local bar_length=50
  local filled=$((completed * bar_length / total))
  local empty=$((bar_length - filled))
  
  # Build progress bar string
  local bar=""
  for ((i=0; i<filled; i++)); do
    bar="${bar}█"
  done
  for ((i=0; i<empty; i++)); do
    bar="${bar}░"
  done
  
  # Print progress bar (overwrite same line)
  printf "\rProgress: [%s] %d/%d (%d%%)" "$bar" "$completed" "$total" "$percent"
}

echo ""  # New line before progress bar starts

# Loop over excitation energies and run all applicable HR modes sequentially
for i in $(seq 0 $((${#excitation_Ens[@]}-1))); do

    en=${excitation_Ens[$i]}
    lbl=$(echo "$en" | awk '{printf "%02dMeV", int($1 + 0.5)}')  # Format: XXMeV (2 sig figs with leading zero, rounded to nearest integer)

    if [ "$i" -ge $HRg_start ] && [ "$i" -le $HRg_stop ]; then
        echo "Running analysis_nChamb.C($i, \"HRg\")"
        root -l -q "analysis_nChamb.C($i, \"HRg\")"
        completed_iterations=$((completed_iterations + 1))
        show_progress $completed_iterations $total_iterations
    fi

    if [ "$i" -ge $HR1n_start ] && [ "$i" -le $HR1n_stop ]; then
        echo "Running analysis_nChamb.C($i, \"HR1n\")"
        root -l -q "analysis_nChamb.C($i, \"HR1n\")"
        completed_iterations=$((completed_iterations + 1))
        show_progress $completed_iterations $total_iterations
    fi

    if [ "$i" -ge $HR2n_start ] && [ "$i" -le $HR2n_stop ]; then
        echo "Running analysis_nChamb.C($i, \"HR2n\")"
        root -l -q "analysis_nChamb.C($i, \"HR2n\")"
        completed_iterations=$((completed_iterations + 1))
        show_progress $completed_iterations $total_iterations
    fi

    if [ "$HR3n_start" -le "$HR3n_stop" ] && [ "$HR3n_stop" -ge 0 ]; then
    if [ "$i" -ge $HR3n_start ] && [ "$i" -le $HR3n_stop ]; then
            echo "Running analysis_nChamb.C($i, \"HR3n\")"
            root -l -q "analysis_nChamb.C($i, \"HR3n\")"
        completed_iterations=$((completed_iterations + 1))
        show_progress $completed_iterations $total_iterations
      fi
    fi

    if [ "$HR4n_start" -le "$HR4n_stop" ] && [ "$HR4n_stop" -ge 0 ]; then
      if [ "$i" -ge $HR4n_start ] && [ "$i" -le $HR4n_stop ]; then
        echo "Running analysis_nChamb.C($i, \"HR4n\")"
        root -l -q "analysis_nChamb.C($i, \"HR4n\")"
        completed_iterations=$((completed_iterations + 1))
        show_progress $completed_iterations $total_iterations
        fi
    fi

done

echo ""  # New line after progress bar completes
