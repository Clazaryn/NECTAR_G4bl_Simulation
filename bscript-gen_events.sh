#!/bin/bash
### --------------------------------------------------------------------------- ###
### ================= Setup for script and parallel analysis ================== ###

# Check if running on correct machine
if [ "$(hostname)" != "borlin305.cenbg.in2p3.fr" ]; then
    echo "###\!/\!/\!/\!/\!/\!/ Caution \!/\!/\!/\!/\!/\!###"
    echo "Error: This script was designed to run on borlin305.cenbg.in2p3.fr. Current hostname: $(hostname)"
    echo "###\!/\!/\!/\!/\!/\!/ Caution \!/\!/\!/\!/\!/\!###"
fi
source /usr/local/Modules/3.2.10/init/bash  # Initialize environment modules for bash
module load python/3.5.2    # Load Python 3.5.2

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
N=14  # number of concurrent jobs

# Check if verbose = True in event_generator.py and exit if so
if grep -q "verbose\s*=\s*True" event_generator.py; then
  echo "ERROR: verbose = True in event_generator.py. Please set verbose = False before running batch script."
  exit 1
fi

# Calculate HR channel ranges and get excitation energies from calc_hr_ranges.py
eval $(python3.5 UtilityScripts/calc_hr_ranges.py)
echo "Excitation energies: ${excitation_Ens[@]}"
echo "Running HR modes: ${run_HR_modes[@]}"

# Build arrays to store mode information for reuse (prevents calculation inconsistencies)
declare -a mode_names=()
declare -a mode_start_indices=()
declare -a mode_stop_indices=()

# Process each mode once and store the information
for mode in "${run_HR_modes[@]}"; do
    mode_clean=$(echo "$mode" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    start_var="${mode_clean}_start"
    stop_var="${mode_clean}_stop"
    if [ -n "${!start_var}" ] && [ -n "${!stop_var}" ]; then
        start_idx=${!start_var}
        stop_idx=${!stop_var}
        mode_names+=("$mode_clean")
        mode_start_indices+=("$start_idx")
        mode_stop_indices+=("$stop_idx")
    fi
done

# Display HR channel ranges using stored information
echo "HR channel ranges:"
for i in "${!mode_names[@]}"; do
    mode_clean=${mode_names[$i]}
    start_idx=${mode_start_indices[$i]}
    stop_idx=${mode_stop_indices[$i]}
    start_energy=${excitation_Ens[$start_idx]}
    stop_energy=${excitation_Ens[$stop_idx]}
    echo "  $mode_clean: indices $start_idx to $stop_idx ($start_energy to $stop_energy MeV)"
done

# Calculate total number of iterations across all channels for progress bar
total_iterations=0
completed_iterations=0
for i in "${!mode_names[@]}"; do
    start_idx=${mode_start_indices[$i]}
    stop_idx=${mode_stop_indices[$i]}
    if [ "$start_idx" -le "$stop_idx" ] && [ "$stop_idx" -ge 0 ]; then
        total_iterations=$((total_iterations + stop_idx - start_idx + 1))
    fi
done

jobs_running=0	# Initialize job count	
job_pids=()		# Track job PIDs


### --------------------------------------------------------------------------- ###
### =============== Function definitions for parallel analysis ================ ###

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

# Function to wait for jobs when limit is reached. If max_jobs is 0, waits for all jobs to complete
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

# Function to check if an excitation energy is in an overlap region between modes
# Overlap occurs at the boundary between consecutive modes (3.0 MeV overlap_range)
is_in_overlap() {
    local exc_index=$1
    local current_mode=$2
    local exc_en=${excitation_Ens[$exc_index]}
    
    # Check if this excitation energy appears in multiple mode ranges
    local count=0
    for i in "${!mode_names[@]}"; do
        local mode_name=${mode_names[$i]}
        local mode_start=${mode_start_indices[$i]}
        local mode_stop=${mode_stop_indices[$i]}
        
        # Check if exc_index is within this mode's range
        if [ "$exc_index" -ge "$mode_start" ] && [ "$exc_index" -le "$mode_stop" ]; then
            count=$((count + 1))
        fi
    done
    
    [ $count -gt 1 ]    # Return true if the excitation energy is in an overlap region
}

# Function to run event generation in background
run_event_generation() {
  local exc_index=$1
  local rec_type=$2
  
  (
    en=${excitation_Ens[$exc_index]}
    lbl=$(echo "$en" | awk '{printf "%04.1fMeV", $1}')  # Format: XX.XMeV (1 decimal place with leading zero)
    
    # Check if this excitation energy is in an overlap region. If so, use half the events to avoid double-counting
    local events_to_use=$nevents
    if is_in_overlap $exc_index "$rec_type"; then
        events_to_use=$((nevents / 2))
    fi
    
    if [ "$rec_type" == "HRf" ]; then
      #__FIRST__Make__the_GEF_.lmd_file
      #__(SECOND__Make the GEF_tree)__if_one_want_to_investigate_it
      enhancement_factor=$(( (nevents + 99999) / 100000 ))
      ./GEFbashscript.sh "$recoil_Z" "$recoil_A" "$en" "$enhancement_factor"

      #__THIRD___generate events !!___
      python3.5 event_generator.py "$nevents" "HRf" "$en" "$lbl" "$enhancement_factor" 
    else
      python3.5 event_generator.py "$events_to_use" "$rec_type" "$en" "$lbl"
    fi
  ) &  # Background job
  
  job_pids+=($!)  # Capture the PID of the current background job
  jobs_running=$((jobs_running + 1))  # Increment job count
  wait_for_jobs  # Wait if we've reached the limit
}

echo ""  # New line before progress bar starts

### --------------------------------------------------------------------------- ###
### ========== Main script execution - executed jobs run in parallel ========== ###

# Create reaction directory structure
mkdir -p ./${reaction}_results/Event_output

# Loop over each HR mode using stored information
for i in "${!mode_names[@]}"; do
    mode_clean=${mode_names[$i]}
    start_idx=${mode_start_indices[$i]}
    stop_idx=${mode_stop_indices[$i]}
    if [ "$start_idx" -le "$stop_idx" ] && [ "$stop_idx" -ge 0 ]; then
        for j in $(seq $start_idx $stop_idx); do
            run_event_generation $j "$mode_clean"
        done
    fi
done

# Wait for all remaining jobs to finish and update progress
wait_for_jobs 0  # Wait until 0 jobs running (i.e., all done)

echo ""  # New line after progress bar completes
