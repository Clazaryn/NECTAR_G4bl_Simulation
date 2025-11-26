#!/bin/bash

reaction=$(grep '^reaction' reac_info_nChamb.txt | awk -F'=' '{gsub(/^ +| +$/,"",$2); print $2}' | awk '{print $1}')

echo "   ####################################################   "
echo "                                                          "
echo "     SIMULATIONS FOR NEW NECTAR REACTION CHAMBER          "
echo "     Python Event Generator for different HR modes        "
echo "     Reaction: $reaction - check scripts!                 "
echo "                                                          "
echo "   ####################################################   "

### User input parameters ###
read -p "Enter number of ejectile events: " nevents < /dev/tty
	
# RUNS IN PARALLEL over the excitation energy range defined in test_event_Eexc.py and TreeAnalysis_eexc.C
# Choose the range i and number of cores N to run over
N=12  # number of concurrent jobs

# Read excitation energies from reac_info_nChamb.txt
excitation_Ens=($(grep '^recoil_excEns' reac_info_nChamb.txt | awk -F'=' '{print $2}' | tr ',' ' '))
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

jobs_running=0	# Initialize job count	
job_pids=()		# Track job PIDs

for i in $(seq $HRg_start $HRg_stop); do		# HRg block
  (
    en=${excitation_Ens[$i]}
    lbl=$(printf "%02d" $(echo "$en" | awk '{print int($1 + 0.5)}'))  # Round to nearest integer
    mkdir -p ../Event_output/$reaction/excEn_$lbl
    python3.5 event_gen_nChamb.py "$nevents" "HRg" "$en" "$lbl"
  ) &  # Background job
  
  job_pids+=($!)		# Capture the PID of the current background job

  jobs_running=$((jobs_running + 1))		# Throttle the number of parallel jobs

  if [ "$jobs_running" -ge "$N" ]; then
    wait ${job_pids[0]}						# Wait for the first job to finish (blocking wait for one job)
    job_pids=("${job_pids[@]:1}")  			# Remove the finished job's PID from the array
    jobs_running=$((jobs_running - 1))  	# Decrease the job count after one finishes
  fi
done

for i in $(seq $HR1n_start $HR1n_stop); do		# HR1n block
  (
    en=${excitation_Ens[$i]}
    lbl=$(printf "%02d" $(echo "$en" | awk '{print int($1 + 0.5)}'))  # Round to nearest integer
    mkdir -p ../Event_output/$reaction/excEn_$lbl
	  python3.5 event_gen_nChamb.py "$nevents" "HR1n" "$en" "$lbl"
  ) &  # Background job
  
  job_pids+=($!)		# Capture the PID of the current background job

  jobs_running=$((jobs_running + 1))		# Throttle the number of parallel jobs

  if [ "$jobs_running" -ge "$N" ]; then
    wait ${job_pids[0]}						# Wait for the first job to finish (blocking wait for one job)
    job_pids=("${job_pids[@]:1}")  			# Remove the finished job's PID from the array
    jobs_running=$((jobs_running - 1))  	# Decrease the job count after one finishes
  fi
done

for i in $(seq $HR2n_start $HR2n_stop); do		# HR2n block
  (
    en=${excitation_Ens[$i]}
    lbl=$(printf "%02d" $(echo "$en" | awk '{print int($1 + 0.5)}'))  # Round to nearest integer
    mkdir -p ../Event_output/$reaction/excEn_$lbl
	  python3.5 event_gen_nChamb.py "$nevents" "HR2n" "$en" "$lbl"
  ) &  # Background job
  
  job_pids+=($!)		# Capture the PID of the current background job

  jobs_running=$((jobs_running + 1))		# Throttle the number of parallel jobs

  if [ "$jobs_running" -ge "$N" ]; then
    wait ${job_pids[0]}						# Wait for the first job to finish (blocking wait for one job)
    job_pids=("${job_pids[@]:1}")  			# Remove the finished job's PID from the array
    jobs_running=$((jobs_running - 1))  	# Decrease the job count after one finishes
  fi
done

if [ "$reaction" == "206Pbdp" ]; then
  for i in $(seq $HR3n_start $HR3n_stop); do		# HR3n block
    (
      en=${excitation_Ens[$i]}
      lbl=$(printf "%02d" $(echo "$en" | awk '{print int($1 + 0.5)}'))  # Round to nearest integer
      mkdir -p ../Event_output/$reaction/excEn_$lbl
  	  python3.5 event_gen_nChamb.py "$nevents" "HR3n" "$en" "$lbl"
    ) &  # Background job

    job_pids+=($!)		# Capture the PID of the current background job

    jobs_running=$((jobs_running + 1))		# Throttle the number of parallel jobs

    if [ "$jobs_running" -ge "$N" ]; then
      wait ${job_pids[0]}						# Wait for the first job to finish (blocking wait for one job)
      job_pids=("${job_pids[@]:1}")  			# Remove the finished job's PID from the array
      jobs_running=$((jobs_running - 1))  	# Decrease the job count after one finishes
    fi
  done
fi

wait	# Wait for all remaining jobs to finish
