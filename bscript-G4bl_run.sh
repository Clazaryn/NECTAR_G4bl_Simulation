#!/bin/bash

# get reaction from info file
reaction=$(grep '^reaction' reac_info.txt | awk -F'=' '{gsub(/^ +| +$/,"",$2); print $2}' | awk '{print $1}')

echo "   ####################################################   "
echo "                                                          "
echo "     SIMULATIONS FOR NEW NECTAR REACTION CHAMBER          "
echo "     G4beamline simulation for different HR modes         "
echo "     Reaction: $reaction - check scripts!                 "
echo "                                                          "
echo "   ####################################################   "

# --------------- Choose particle type ---------------
#loop_list=("ejectile")
#loop_list=("recoil")
loop_list=("ejectile" "recoil")

# RUNS IN PARALLEL over the excitation energy range defined in reac_info.txt
# Choose the range i and number of cores N to run over
N=15  # number of concurrent jobs

# get ring parameters from info file
det_setup=$(grep '^det_setup' reac_info.txt | awk -F'=' '{gsub(/^ +| +$/,"",$2); print $2}' | awk '{print $1}')
Bfield=$(grep '^dip_Bfld' reac_info.txt | awk -F'=' '{gsub(/^ +| +$/,"",$2); print $2}' | awk '{print $1}')
Q1=$(grep '^quad_Q1' reac_info.txt | awk -F'=' '{gsub(/^ +| +$/,"",$2); print $2}' | awk '{print $1}')
Q2=$(grep '^quad_Q2' reac_info.txt | awk -F'=' '{gsub(/^ +| +$/,"",$2); print $2}' | awk '{print $1}')
Q3=$(grep '^quad_Q3' reac_info.txt | awk -F'=' '{gsub(/^ +| +$/,"",$2); print $2}' | awk '{print $1}')
Q4=$(grep '^quad_Q4' reac_info.txt | awk -F'=' '{gsub(/^ +| +$/,"",$2); print $2}' | awk '{print $1}')
Q5=$(grep '^quad_Q5' reac_info.txt | awk -F'=' '{gsub(/^ +| +$/,"",$2); print $2}' | awk '{print $1}')

# choose detector setup based on det_setup in reac_info.txt
if [ "$det_setup" == "new" ]; then
  echo "Using new ESR detector setup"
  input_file="ESR_ring-New_Detectors.in"
elif [ "$det_setup" == "PoP" ]; then
  echo "Using PoP ESR detector setup"
  input_file="ESR_ring-PoP_Detectors.in"
  echo "Magnetic field scaling set manually: check value is correct in input file"
else
  echo "Invalid detector setup"
  exit 1
fi

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

jobs_running=0	# Initialize job count	
job_pids=()		# Track job PIDs

# Calculate total number of iterations per particle type across all channels
iterations_per_particle=0
iterations_per_particle=$((iterations_per_particle + HRg_stop - HRg_start + 1))
iterations_per_particle=$((iterations_per_particle + HR1n_stop - HR1n_start + 1))
if [ "$HR2n_start" -le "$HR2n_stop" ] && [ "$HR2n_stop" -ge 0 ]; then
  iterations_per_particle=$((iterations_per_particle + HR2n_stop - HR2n_start + 1))
fi
if [ "$HR3n_start" -le "$HR3n_stop" ] && [ "$HR3n_stop" -ge 0 ]; then
  iterations_per_particle=$((iterations_per_particle + HR3n_stop - HR3n_start + 1))
fi
if [ "$HR4n_start" -le "$HR4n_stop" ] && [ "$HR4n_stop" -ge 0 ]; then
  iterations_per_particle=$((iterations_per_particle + HR4n_stop - HR4n_start + 1))
fi

# Total iterations = iterations per particle type × number of particle types
num_particle_types=${#loop_list[@]}
total_iterations=$((iterations_per_particle * num_particle_types))

# Separate counters for each particle type
completed_ejectile=0
completed_recoil=0

# Function to display progress bar for a specific particle type
show_progress() {
  local particle_type=$1
  local completed=$2
  local total=$3
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
  
  # Print progress bar with particle type label
  printf "\r%s Progress: [%s] %d/%d (%d%%)" "$particle_type" "$bar" "$completed" "$total" "$percent"
}

# Function to wait for jobs and update progress
# If max_jobs is 0, waits for all jobs to complete
# Keeps particle-specific logic for timing scripts and progress tracking
wait_for_jobs() {
  local particle_type=$1
  local max_jobs=${2:-$N}  # Default to N, but can be overridden (use 0 for all jobs)
  
  while [ "$jobs_running" -ge "$max_jobs" ] && [ ${#job_pids[@]} -gt 0 ]; do
    # Wait for the oldest job (first in array) to finish
    local finished_pid=${job_pids[0]}
    wait $finished_pid 2>/dev/null
    
    # Remove the finished job from the array
    job_pids=("${job_pids[@]:1}")
    jobs_running=$((jobs_running - 1))
    
    # If this was the first job, stop the timing script
    if [ "${particle_type}" == "ejectile" ] && [ "$finished_pid" == "$first_ejectile_job_pid" ] && [ -n "$timing_ejectile_pid" ]; then
      # Check if process still exists before killing to avoid error messages
      if kill -0 $timing_ejectile_pid 2>/dev/null; then
        kill $timing_ejectile_pid >/dev/null 2>&1
        wait $timing_ejectile_pid 2>/dev/null
      fi
      timing_ejectile_pid=""
    elif [ "${particle_type}" == "recoil" ] && [ "$finished_pid" == "$first_recoil_job_pid" ] && [ -n "$timing_recoil_pid" ]; then
      # Check if process still exists before killing to avoid error messages
      if kill -0 $timing_recoil_pid 2>/dev/null; then
        kill $timing_recoil_pid >/dev/null 2>&1
        wait $timing_recoil_pid 2>/dev/null
      fi
      timing_recoil_pid=""
    fi
    
    # Update appropriate counter based on particle type
    if [ "${particle_type}" == "recoil" ]; then
      completed_recoil=$((completed_recoil + 1))
      show_progress "Recoil" $completed_recoil $iterations_per_particle
      recoil_remaining=$((recoil_remaining - 1))
    else
      completed_ejectile=$((completed_ejectile + 1))
      show_progress "Ejectile" $completed_ejectile $iterations_per_particle
    fi
  done
}

# Function to run G4beamline simulation for a single excitation energy
# Parameters: exc_index, rec_type (HRg, HR1n, etc.), particle (ejectile/recoil)
run_G4bl_sim() {
  local exc_index=$1
  local rec_type=$2
  local particle=$3
  
  en=${excitation_Ens[$exc_index]}
  lbl=$(echo "$en" | awk '{printf "%02dMeV", int($1 + 0.5)}')  # Format: XXMeV (2 sig figs with leading zero, rounded to nearest integer)
  
  # Format the file names in bash (outside subshell for HRg so they're available for Python script)
  histoFile="Detectors_${reaction}_${rec_type}_excEn${lbl}_${particle}.root"
  outFile="g4bl_${reaction}_${rec_type}_excEn${lbl}_${particle}.out"
  eventFile="../${reaction}_sim/Event_output/output_event_generator_${reaction}_${rec_type}_excEn${lbl}_${particle}.txt"
  
  (
    # run G4beamline simulation
    g4bl $input_file histoFile=$histoFile eventFile=$eventFile Bfield=$Bfield Q1=$Q1 Q2=$Q2 Q3=$Q3 Q4=$Q4 Q5=$Q5 > $outFile 2>&1
  ) &  # Background job
  
  current_pid=$!
  job_pids+=($current_pid)  # Capture the PID of the current background job
  
  # Track first job for timing (separate for ejectile and recoil) - only for first block (HRg)
  if [ "${rec_type}" == "HRg" ]; then
    if [ "${particle}" == "ejectile" ] && [ -z "$first_ejectile_job_pid" ]; then
      first_ejectile_job_pid=$current_pid
      first_ejectile_start_time=$(date +%s)
      remaining_iters=$((ejectile_remaining - 1))
      # Start Python timing estimation in background and track its PID
      python3.5 -u UtilityScripts/estimate_recoil_timing.py "$outFile" "$eventFile" "$remaining_iters" "$N" &
      timing_ejectile_pid=$!
    elif [ "${particle}" == "recoil" ] && [ -z "$first_recoil_job_pid" ]; then
      first_recoil_job_pid=$current_pid
      first_recoil_start_time=$(date +%s)
      remaining_iters=$((recoil_remaining - 1))
      # Start Python timing estimation in background and track its PID
      python3.5 -u UtilityScripts/estimate_recoil_timing.py "$outFile" "$eventFile" "$remaining_iters" "$N" &
      timing_recoil_pid=$!
    fi
  fi
  
  jobs_running=$((jobs_running + 1))  # Increment job count
  
  # Wait for jobs if we've reached the limit
  wait_for_jobs "$particle"
}

# Loop over all particle types in loop_list
for particle in ${loop_list[@]}; do

if [ ${#loop_list[@]} -eq 1 ]; then
  echo "You have chosen to simulate only the ${particle}s. Check if you already have corresponding recoils!"
fi

if [ ${particle} == "ejectile" ]; then
  echo ""
  echo "Starting simulations for ${particle} particles (fast)..."
  # Initialize timing tracking variables for ejectile
  first_ejectile_job_pid=""
  first_ejectile_start_time=""
  timing_ejectile_pid=""
  ejectile_remaining=$iterations_per_particle
elif [ ${particle} == "recoil" ]; then
  echo ""
  echo "Starting simulations for ${particle} particles (quite slow due to magnetic fields)..."
  # Initialize timing tracking variables for recoil
  first_recoil_job_pid=""
  first_recoil_start_time=""
  timing_recoil_pid=""
  # Calculate remaining iterations for recoil particles
  recoil_remaining=$iterations_per_particle
fi

# HRg block
for i in $(seq $HRg_start $HRg_stop); do
  run_G4bl_sim $i "HRg" "$particle"
done

# HR1n block
for i in $(seq $HR1n_start $HR1n_stop); do
  run_G4bl_sim $i "HR1n" "$particle"
done

# HR2n block (only if range is valid)
if [ "$HR2n_start" -le "$HR2n_stop" ] && [ "$HR2n_stop" -ge 0 ]; then
  for i in $(seq $HR2n_start $HR2n_stop); do
    run_G4bl_sim $i "HR2n" "$particle"
  done
fi

# HR3n block (only if range is valid)
if [ "$HR3n_start" -le "$HR3n_stop" ] && [ "$HR3n_stop" -ge 0 ]; then
  for i in $(seq $HR3n_start $HR3n_stop); do
    run_G4bl_sim $i "HR3n" "$particle"
  done
fi

# HR4n block (only if range is valid)
if [ "$HR4n_start" -le "$HR4n_stop" ] && [ "$HR4n_stop" -ge 0 ]; then
  for i in $(seq $HR4n_start $HR4n_stop); do
    run_G4bl_sim $i "HR4n" "$particle"
  done
fi

# Wait for all remaining jobs for this particle type to finish and update progress
wait_for_jobs "$particle" 0  # Wait until 0 jobs running (i.e., all done)

# Reset job tracking variables for the next particle type
jobs_running=0
job_pids=()

echo ""
echo "Completed all simulations for ${particle} particles."

done  # End of particle type loop

echo ""  # New line after progress bar completes

# Create output directory if it doesn't exist
mkdir -p ../${reaction}_sim/Detector_output/

# Moves the output of G4beamline to the right folder and renames it
# Use shopt to handle cases where no files match the pattern
shopt -s nullglob
out_files=(g4bl*.out)
root_files=(Detectors*.root)
shopt -u nullglob

# Move .out files if any exist
if [ ${#out_files[@]} -gt 0 ]; then
  mv "${out_files[@]}" ../${reaction}_sim/Detector_output/ 2>/dev/null || true
fi

# Move .root files if any exist
if [ ${#root_files[@]} -gt 0 ]; then
  mv "${root_files[@]}" ../${reaction}_sim/Detector_output/ 2>/dev/null || true
fi
