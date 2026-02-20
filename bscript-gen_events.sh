#!/bin/bash

# Read reaction from reac_info.txt
reaction=$(grep '^reaction' reac_info.txt | awk -F'=' '{gsub(/^ +| +$/,"",$2); print $2}' | awk '{print $1}')
# Read recoil Z and A too, needed for GEF input
recoil_A=$(grep '^recoil_A' reac_info.txt | awk -F'=' '{gsub(/^ +| +$/,"",$2); print $2}' | awk '{print $1}')
recoil_Z=$(grep '^recoil_Z' reac_info.txt | awk -F'=' '{gsub(/^ +| +$/,"",$2); print $2}' | awk '{print $1}')


echo "   ####################################################   "
echo "                                                          "
echo "     NECTAR SIMULATION SUITE - Event Generator            "
echo "     Python Event Generator for different HR modes        "
echo "     Reaction: $reaction - check scripts!                 "
echo "                                                          "
echo "   ####################################################   "

### User input parameters ###
read -p "Enter number of ejectile events: " nevents < /dev/tty
	
### RUNS IN PARALLEL over the excitation energy range defined in reac_info.txt ###
N=12  # number of concurrent jobs

# Check if verbose = True in event_generator.py and exit if so
if grep -q "verbose\s*=\s*True" event_generator.py; then
  echo "ERROR: verbose = True in event_generator.py. Please set verbose = False before running batch script."
  exit 1
fi

# Read excitation energies from reac_info.txt
excitation_Ens=($(grep '^recoil_excEns' reac_info.txt | awk -F'=' '{print $2}' | tr ',' ' '))
echo "Excitation energies: ${excitation_Ens[@]}"

#Get the path to the fission barrier
export BJORNHOLM_BARRIER_FILE="NuclideDataMaster/Bjornholm_Fission_Barrier.txt"

# Calculate HR channel ranges based on separation energies (3 MeV beyond next channel opening)
hr_vars="$(python3 UtilityScripts/calc_hr_ranges.py)" || exit 1
eval "$hr_vars"
echo "HR channel ranges:"
echo "  HRg: indices $HRg_start to $HRg_stop (${excitation_Ens[$HRg_start]} to ${excitation_Ens[$HRg_stop]} MeV)"
echo "  HR1n: indices $HR1n_start to $HR1n_stop (${excitation_Ens[$HR1n_start]} to ${excitation_Ens[$HR1n_stop]} MeV)"
echo "  HR2n: indices $HR2n_start to $HR2n_stop (${excitation_Ens[$HR2n_start]} to ${excitation_Ens[$HR2n_stop]} MeV)"
echo "  HR3n: indices $HR3n_start to $HR3n_stop (${excitation_Ens[$HR3n_start]} to ${excitation_Ens[$HR3n_stop]} MeV)"
echo "  HR4n: indices $HR4n_start to $HR4n_stop (${excitation_Ens[$HR4n_start]} to ${excitation_Ens[$HR4n_stop]} MeV)"
echo "  Z = $recoil_Z and A = $recoil_A,  Fission barrier Ua: $Ua MeV"
echo "  HRf: indices $HRf_start to $HRf_stop (${excitation_Ens[$HRf_start]} to ${excitation_Ens[$HRf_stop]} MeV)"
echo "  Max excitation energy: $max_excEn MeV"

jobs_running=0	# Initialize job count	
job_pids=()		# Track job PIDs

# Calculate total number of iterations across all channels for progress bar
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

# Create reaction directory structure
mkdir -p ../${reaction}_sim/Event_output

# HERE COMES THE FISSION BLOCK ____ 
# TEST___BLOCK__SO__FAR____________

#for i in $(seq $HRf_start $HRf_start); do		# HRf block
for i in $(seq 6 6); do		# HRf block
  (
    en=${excitation_Ens[$i]}
    lbl=$(echo "$en" | awk '{printf "%02dMeV", int($1 + 0.5)}')  # Format: XXMeV (2 sig figs with leading zero, rounded to nearest integer)
    
    #__FIRST__Make__the_GEF_.lmd_file
    #__(SECOND__Make the GEF_tree)__if_one_want_to_investigate_it
    enhancement_factor=$(( (nevents + 99999) / 100000 ))
    ./GEFbashscript.sh "$recoil_Z" "$recoil_A" "$en" "$enhancement_factor"

    #__THIRD___generate events !!___
    python3 event_generator.py "$nevents" "HRf" "$en" "$lbl" "$enhancement_factor"
  ) &  # Background job
  
  job_pids+=($!)		# Capture the PID of the current background job

  jobs_running=$((jobs_running + 1))		# Throttle the number of parallel jobs

  # Wait for a job to finish if we've reached the limit
  while [ "$jobs_running" -ge "$N" ]; do
    # Wait for the oldest job (first in array) to finish
    wait ${job_pids[0]} 2>/dev/null
    # Remove the finished job from the array
    job_pids=("${job_pids[@]:1}")
    jobs_running=$((jobs_running - 1))
    completed_iterations=$((completed_iterations + 1))
    show_progress $completed_iterations $total_iterations
  done
done


# END__OF_THE_TEST_BLOCK____________


for i in $(seq $HRg_start $HRg_stop); do		# HRg block
  (
    en=${excitation_Ens[$i]}
    lbl=$(echo "$en" | awk '{printf "%02dMeV", int($1 + 0.5)}')  # Format: XXMeV (2 sig figs with leading zero, rounded to nearest integer)
    python3 event_generator.py "$nevents" "HRg" "$en" "$lbl"
  ) &  # Background job
  
  job_pids+=($!)		# Capture the PID of the current background job

  jobs_running=$((jobs_running + 1))		# Throttle the number of parallel jobs

  # Wait for a job to finish if we've reached the limit
  while [ "$jobs_running" -ge "$N" ]; do
    # Wait for the oldest job (first in array) to finish
    wait ${job_pids[0]} 2>/dev/null
    # Remove the finished job from the array
    job_pids=("${job_pids[@]:1}")
    jobs_running=$((jobs_running - 1))
    completed_iterations=$((completed_iterations + 1))
    show_progress $completed_iterations $total_iterations
  done
done

for i in $(seq $HR1n_start $HR1n_stop); do		# HR1n block
  (
    en=${excitation_Ens[$i]}
    lbl=$(echo "$en" | awk '{printf "%02dMeV", int($1 + 0.5)}')  # Format: XXMeV (2 sig figs with leading zero, rounded to nearest integer)
	  python3 event_generator.py "$nevents" "HR1n" "$en" "$lbl"
  ) &  # Background job
  
  job_pids+=($!)		# Capture the PID of the current background job

  jobs_running=$((jobs_running + 1))		# Throttle the number of parallel jobs

  # Wait for a job to finish if we've reached the limit
  while [ "$jobs_running" -ge "$N" ]; do
    # Wait for the oldest job (first in array) to finish
    wait ${job_pids[0]} 2>/dev/null
    # Remove the finished job from the array
    job_pids=("${job_pids[@]:1}")
    jobs_running=$((jobs_running - 1))
    completed_iterations=$((completed_iterations + 1))
    show_progress $completed_iterations $total_iterations
  done
done

# HR2n block (only if range is valid)
if [ "$HR2n_start" -le "$HR2n_stop" ] && [ "$HR2n_stop" -ge 0 ]; then
  for i in $(seq $HR2n_start $HR2n_stop); do		# HR2n block
    (
      en=${excitation_Ens[$i]}
      lbl=$(echo "$en" | awk '{printf "%02dMeV", int($1 + 0.5)}')  # Format: XXMeV (2 sig figs with leading zero, rounded to nearest integer)
  	  python3 event_generator.py "$nevents" "HR2n" "$en" "$lbl"
    ) &  # Background job

    job_pids+=($!)		# Capture the PID of the current background job

    jobs_running=$((jobs_running + 1))		# Throttle the number of parallel jobs

    # Wait for a job to finish if we've reached the limit
    while [ "$jobs_running" -ge "$N" ]; do
      # Wait for the oldest job (first in array) to finish
      wait ${job_pids[0]} 2>/dev/null
      # Remove the finished job from the array
      job_pids=("${job_pids[@]:1}")
      jobs_running=$((jobs_running - 1))
      completed_iterations=$((completed_iterations + 1))
      show_progress $completed_iterations $total_iterations
    done
  done
fi

# HR3n block (only if range is valid)
if [ "$HR3n_start" -le "$HR3n_stop" ] && [ "$HR3n_stop" -ge 0 ]; then
  for i in $(seq $HR3n_start $HR3n_stop); do
    (
      en=${excitation_Ens[$i]}
      lbl=$(echo "$en" | awk '{printf "%02dMeV", int($1 + 0.5)}')  # Format: XXMeV (2 sig figs with leading zero, rounded to nearest integer)
      python3 event_generator.py "$nevents" "HR3n" "$en" "$lbl"
    ) &  # Background job

    job_pids+=($!)		# Capture the PID of the current background job

    jobs_running=$((jobs_running + 1))		# Throttle the number of parallel jobs

    # Wait for a job to finish if we've reached the limit
    while [ "$jobs_running" -ge "$N" ]; do
      # Wait for the oldest job (first in array) to finish
      wait ${job_pids[0]} 2>/dev/null
      # Remove the finished job from the array
      job_pids=("${job_pids[@]:1}")
      jobs_running=$((jobs_running - 1))
      completed_iterations=$((completed_iterations + 1))
      show_progress $completed_iterations $total_iterations
    done
  done
fi

# HR4n block (only if range is valid and Sn_4nDght is defined)
if [ "$HR4n_start" -le "$HR4n_stop" ] && [ "$HR4n_stop" -ge 0 ]; then
  for i in $(seq $HR4n_start $HR4n_stop); do
    (
      en=${excitation_Ens[$i]}
      lbl=$(echo "$en" | awk '{printf "%02dMeV", int($1 + 0.5)}')  # Format: XXMeV (2 sig figs with leading zero, rounded to nearest integer)
      python3 event_generator.py "$nevents" "HR4n" "$en" "$lbl"
    ) &  # Background job

    job_pids+=($!)		# Capture the PID of the current background job

    jobs_running=$((jobs_running + 1))		# Throttle the number of parallel jobs

    # Wait for a job to finish if we've reached the limit
    while [ "$jobs_running" -ge "$N" ]; do
      # Wait for the oldest job (first in array) to finish
      wait ${job_pids[0]} 2>/dev/null
      # Remove the finished job from the array
      job_pids=("${job_pids[@]:1}")
      jobs_running=$((jobs_running - 1))
      completed_iterations=$((completed_iterations + 1))
      show_progress $completed_iterations $total_iterations
    done
  done
fi

# Wait for all remaining jobs to finish and update progress
while [ ${#job_pids[@]} -gt 0 ]; do
  wait ${job_pids[0]}
  job_pids=("${job_pids[@]:1}")
  completed_iterations=$((completed_iterations + 1))
  show_progress $completed_iterations $total_iterations
done

echo ""  # New line after progress bar completes
