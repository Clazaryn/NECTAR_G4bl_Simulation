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

# Loop over all particle types in loop_list
for particle in ${loop_list[@]}; do

if [ ${#loop_list[@]} -eq 1 ]; then
  echo "You have chosen to simulate only the ${particle}s. Check if you already have corresponding recoils!"
fi

if [ ${particle} == "ejectile" ]; then
  echo "Starting simulations for ${particle} particles (fast)..."
elif [ ${particle} == "recoil" ]; then
  echo "Starting simulations for ${particle} particles (quite slow due to magnetic fields)..."
fi

for i in $(seq $HRg_start $HRg_stop); do		# HRg block
  (
    en=${excitation_Ens[$i]}
    lbl=$(echo "$en" | awk '{printf "%02dMeV", int($1 + 0.5)}')  # Format: XXMeV (2 sig figs with leading zero, rounded to nearest integer)

    # Format the file names in bash
    histoFile="Detectors_${reaction}_HRg_excEn${lbl}_${particle}.root"
    outFile="g4bl_${reaction}_HRg_excEn${lbl}_${particle}.out"
    eventFile="../${reaction}_sim/Event_output/output_event_generator_${reaction}_HRg_excEn${lbl}_${particle}.txt"

    g4bl $input_file histoFile=$histoFile eventFile=$eventFile Bfield=$Bfield Q1=$Q1 Q2=$Q2 Q3=$Q3 Q4=$Q4 Q5=$Q5 > $outFile 2>&1
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
    
    # Format the file names in bash
    histoFile="Detectors_${reaction}_HR1n_excEn${lbl}_${particle}.root"
    outFile="g4bl_${reaction}_HR1n_excEn${lbl}_${particle}.out"
    eventFile="../${reaction}_sim/Event_output/output_event_generator_${reaction}_HR1n_excEn${lbl}_${particle}.txt"

    g4bl $input_file histoFile=$histoFile eventFile=$eventFile Bfield=$Bfield Q1=$Q1 Q2=$Q2 Q3=$Q3 Q4=$Q4 Q5=$Q5 > $outFile 2>&1
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
      
      # Format the file names in bash
      histoFile="Detectors_${reaction}_HR2n_excEn${lbl}_${particle}.root"
      outFile="g4bl_${reaction}_HR2n_excEn${lbl}_${particle}.out"
      eventFile="../${reaction}_sim/Event_output/output_event_generator_${reaction}_HR2n_excEn${lbl}_${particle}.txt"

      g4bl $input_file histoFile=$histoFile eventFile=$eventFile Bfield=$Bfield Q1=$Q1 Q2=$Q2 Q3=$Q3 Q4=$Q4 Q5=$Q5 > $outFile 2>&1
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
  for i in $(seq $HR3n_start $HR3n_stop); do		# HR3n block
    (
      en=${excitation_Ens[$i]}
      lbl=$(echo "$en" | awk '{printf "%02dMeV", int($1 + 0.5)}')  # Format: XXMeV (2 sig figs with leading zero, rounded to nearest integer)
      
      # Format the file names in bash
      histoFile="Detectors_${reaction}_HR3n_excEn${lbl}_${particle}.root"
      outFile="g4bl_${reaction}_HR3n_excEn${lbl}_${particle}.out"
      eventFile="../${reaction}_sim/Event_output/output_event_generator_${reaction}_HR3n_excEn${lbl}_${particle}.txt"

      g4bl $input_file histoFile=$histoFile eventFile=$eventFile Bfield=$Bfield Q1=$Q1 Q2=$Q2 Q3=$Q3 Q4=$Q4 Q5=$Q5 > $outFile 2>&1
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
      
      # Format the file names in bash
      histoFile="Detectors_${reaction}_HR4n_excEn${lbl}_${particle}.root"
      outFile="g4bl_${reaction}_HR4n_excEn${lbl}_${particle}.out"
      eventFile="../${reaction}_sim/Event_output/output_event_generator_${reaction}_HR4n_excEn${lbl}_${particle}.txt"

      g4bl $input_file histoFile=$histoFile eventFile=$eventFile Bfield=$Bfield Q1=$Q1 Q2=$Q2 Q3=$Q3 Q4=$Q4 Q5=$Q5 > $outFile 2>&1
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

# Wait for all remaining jobs for this particle type to finish and update progress
while [ ${#job_pids[@]} -gt 0 ]; do
  wait ${job_pids[0]} 2>/dev/null
  job_pids=("${job_pids[@]:1}")
  completed_iterations=$((completed_iterations + 1))
  show_progress $completed_iterations $total_iterations
done

# Reset job tracking variables for the next particle type
jobs_running=0
job_pids=()

echo "Completed all simulations for ${particle} particles."

done  # End of particle type loop

echo ""  # New line after progress bar completes

# Create output directory if it doesn't exist
mkdir -p ../${reaction}_sim/Detector_output/

# Moves the output of G4beamline to the right folder and renames it
mv g4bl*.out ../${reaction}_sim/Detector_output/
mv Detectors*.root ../${reaction}_sim/Detector_output/
