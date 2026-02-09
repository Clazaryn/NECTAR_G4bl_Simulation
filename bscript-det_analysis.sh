#!/bin/bash

reaction=$(grep '^reaction' reac_info.txt | awk -F'=' '{gsub(/^ +| +$/,"",$2); print $2}' | awk '{print $1}')

echo "   ####################################################   "
echo "                                                          "
echo "     SIMULATIONS FOR NEW NECTAR REACTION CHAMBER          "
echo "     Detector analysis - TTree output for HR modes       "
echo "     Reaction: $reaction - check scripts!                 "
echo "                                                          "
echo "   ####################################################   "

# Precompile the analysis library
echo "Precompiling analysis library..."
g++ -shared -fPIC -std=c++11 -I. \
  UtilityScripts/det_analysis_impl.cxx det_analysis.cxx \
  -o libDetAnalysis.so \
  $(root-config --cflags --libs) 2>&1

if [ $? -ne 0 ]; then
  echo "Error: Failed to compile analysis library"
  exit 1
fi

echo "Analysis library compiled successfully"
echo ""

### RUNS IN PARALLEL over the excitation energy range defined in reac_info.txt ###
N=12  # number of concurrent jobs

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
echo "  HR4n: indices $HR4n_start to $HR4n_stop (${excitation_Ens[$HR4n_start]} to ${excitation_Ens[$HR4n_stop]} MeV)"
echo "  Max excitation energy: $max_excEn MeV"
echo ""

# Create output directory
mkdir -p ../${reaction}_sim/Det_analysis

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

jobs_running=0  # Initialize job count
job_pids=()     # Track job PIDs

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

# Function to wait for jobs when limit is reached
# If max_jobs is 0, waits for all jobs to complete
wait_for_jobs() {
  local max_jobs=${1:-$N}  # Default to N, but can be overridden (use 0 for all jobs)
  
  while [ "$jobs_running" -ge "$max_jobs" ] && [ ${#job_pids[@]} -gt 0 ]; do
    # Wait for the oldest job (first in array) to finish
    wait ${job_pids[0]} 2>/dev/null
    # Remove the finished job from the array
    job_pids=("${job_pids[@]:1}")
    jobs_running=$((jobs_running - 1))
    completed_iterations=$((completed_iterations + 1))
    show_progress $completed_iterations $total_iterations
  done
}

# Function to run ROOT analysis in background: 
# root loads the library and the header file, then executes the det_analysis function
run_analysis() {
  local exc_index=$1
  local rec_type=$2
  
  # Get script directory to ensure ROOT runs from correct location
  local script_dir=$(cd "$(dirname "$0")" && pwd)
  local lib_path="${script_dir}/libDetAnalysis.so"
  
  (
    # Change to script directory so relative paths work correctly
    cd "$script_dir"
    # Run ROOT - redirect both stdout and stderr to log file
    root -l -b <<ROOT_EOF &>../${reaction}_sim/Det_analysis/root_output_excEn${exc_index}_${rec_type}.log
gSystem->Load("$lib_path");
.L UtilityScripts/det_analysis.h
det_analysis($exc_index,"$rec_type");
.q
ROOT_EOF
  ) &  # Background job
  
  job_pids+=($!)  # Capture the PID of the current background job
  jobs_running=$((jobs_running + 1))  # Increment job count
  wait_for_jobs  # Wait if we've reached the limit
}

# HRg block
for i in $(seq $HRg_start $HRg_stop); do
  run_analysis $i "HRg"
done

# HR1n block
for i in $(seq $HR1n_start $HR1n_stop); do
  run_analysis $i "HR1n"
done

# HR2n block (only if range is valid)
if [ "$HR2n_start" -le "$HR2n_stop" ] && [ "$HR2n_stop" -ge 0 ]; then
  for i in $(seq $HR2n_start $HR2n_stop); do
    run_analysis $i "HR2n"
  done
fi

# HR3n block (only if range is valid)
if [ "$HR3n_start" -le "$HR3n_stop" ] && [ "$HR3n_stop" -ge 0 ]; then
  for i in $(seq $HR3n_start $HR3n_stop); do
    run_analysis $i "HR3n"
  done
fi

# HR4n block (only if range is valid)
if [ "$HR4n_start" -le "$HR4n_stop" ] && [ "$HR4n_stop" -ge 0 ]; then
  for i in $(seq $HR4n_start $HR4n_stop); do
    run_analysis $i "HR4n"
  done
fi

# Wait for all remaining jobs to finish and update progress
wait_for_jobs 0  # Wait until 0 jobs running (i.e., all done)

echo ""  # New line after progress bar completes

# Cleanup compiled library
if [ -f libDetAnalysis.so ]; then
    rm libDetAnalysis.so
    echo "Cleaned up compiled library"
fi
