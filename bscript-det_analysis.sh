#!/bin/bash
### --------------------------------------------------------------------------- ###
### ================= Setup for script and parallel analysis ================== ###

reaction=$(grep '^reaction' reac_info.txt | awk -F'=' '{gsub(/^ +| +$/,"",$2); print $2}' | awk '{print $1}')

echo "   ####################################################   "
echo "                                                          "
echo "     SIMULATIONS FOR NEW NECTAR REACTION CHAMBER          "
echo "     Detector analysis - TTree output for HR modes       "
echo "     Reaction: $reaction - check scripts!                 "
echo "                                                          "
echo "   ####################################################   "

### RUNS IN PARALLEL over the excitation energy range defined in reac_info.txt ###
N=14  # number of concurrent jobs

# Calculate HR channel ranges and get excitation energies from calc_hr_ranges.py
eval $(python3 UtilityScripts/calc_hr_ranges.py)
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
echo ""

# Calculate total number of iterations across all channels
total_iterations=0
for i in "${!mode_names[@]}"; do
    start_idx=${mode_start_indices[$i]}
    stop_idx=${mode_stop_indices[$i]}
    if [ "$start_idx" -le "$stop_idx" ] && [ "$stop_idx" -ge 0 ]; then
        total_iterations=$((total_iterations + stop_idx - start_idx + 1))
    fi
done

completed_iterations=0  # Track completed iterations

jobs_running=0  # Initialize job count
job_pids=()     # Track job PIDs

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

# Function to wait for jobs when limit is reached
# If max_jobs is 0, waits for all jobs to complete
# Exits the entire script immediately if any job fails
wait_for_jobs() {
  local max_jobs=${1:-$N}  # Default to N, but can be overridden (use 0 for all jobs)
  
  while [ "$jobs_running" -ge "$max_jobs" ] && [ ${#job_pids[@]} -gt 0 ]; do
    # Wait for the oldest job (first in array) to finish
    wait ${job_pids[0]} 2>/dev/null
    local exit_status=$?
    # Remove the finished job from the array
    job_pids=("${job_pids[@]:1}")
    jobs_running=$((jobs_running - 1))
    completed_iterations=$((completed_iterations + 1))
    show_progress $completed_iterations $total_iterations
    
    # Check if job failed - exit immediately if so
    if [ $exit_status -ne 0 ]; then
      echo ""
      echo "Error: Analysis job failed - check log file for details. Exiting script."
      exit 1
    fi
  done
}

# Helper function to compile analysis and plotting libraries only if source files are newer than the output file
compile_library_if_needed() {
    local output_file=$1
    shift  # Remove first argument, rest are arguments (source files and linker flags)
    local compile_args=("$@")
    
    # Separate source files from linker flags
    local source_files=()
    local linker_flags=()
    for arg in "${compile_args[@]}"; do
        if [[ "$arg" == -* ]]; then
            linker_flags+=("$arg")
        else
            source_files+=("$arg")
        fi
    done
    
    local needs_compile=0
    
    if [ ! -f "$output_file" ]; then    # if no source file, compile
        needs_compile=1
    else
        # Check if any source file is newer than output
        for source_file in "${source_files[@]}"; do
            if [ "$source_file" -nt "$output_file" ]; then
                needs_compile=1
                break
            fi
        done
        # Check if libDetAnalysis.so is newer (for libMakePlots.so dependency)
        if [[ "$output_file" == *"MakePlots"* ]] && [ -f "libDetAnalysis.so" ] && [ "libDetAnalysis.so" -nt "$output_file" ]; then
            needs_compile=1
        fi
    fi
    
    if [ $needs_compile -eq 1 ]; then
        echo "Compiling -> $output_file (g++)..."
        g++ -shared -fPIC -std=c++11 -I. \
          -Wl,-rpath,'$ORIGIN' \
          "${source_files[@]}" \
          -o "$output_file" \
          "${linker_flags[@]}" \
          $(root-config --cflags --libs) 2>&1
        # Compile with g++:
        # -shared: Create a shared library (.so file) instead of an executable
        # -fPIC: Generate Position-Independent Code (required for shared libraries)
        # -std=c++11: Use C++11 language standard
        # -I.: Add current directory to include path for finding header files
        # -Wl,-rpath,'$ORIGIN': Set runtime library search path to directory where .so file is located
        # "${source_files[@]}": Array of source files (.cxx files) to compile
        # -o "$output_file": Specify output file name
        # "${linker_flags[@]}": Additional linker flags (e.g., -L. -lDetAnalysis)
        # $(root-config --cflags --libs): Get ROOT compilation flags and library paths
        # 2>&1: Redirect stderr to stdout so compilation errors appear in the same stream
        
        if [ $? -ne 0 ]; then
            echo "Error: Failed to compile $output_file"
            exit 1
        fi
        echo "Library compiled successfully"
    else
        echo "Skipping compilation of $output_file (up to date)"
    fi
}

# Function to run ROOT analysis in background: 
# root loads the library and the header file, then executes the det_analysis function
run_analysis() {
  local exc_index=$1
  local rec_type=$2
  
  # Get excitation energy value and format as label (matches det_analysis.cxx format: XX.XMeV)
  local en=${excitation_Ens[$exc_index]}
  local lbl=$(echo "$en" | awk '{printf "%04.1fMeV", $1}')  # Format: XX.XMeV (1 decimal place with leading zero)
  
  (
    # Change to script directory so relative paths work correctly
    # Run ROOT - redirect both stdout and stderr to single log file
    root -l -b <<ROOT_EOF &>../${reaction}_sim/Det_analysis/root_output_excEn${lbl}_${rec_type}.log
gSystem->Load("libDetAnalysis.so");
.L UtilityScripts/det_analysis.h
det_analysis($exc_index,"$rec_type");
.q
ROOT_EOF
    # Check if ROOT exited with error (non-zero exit code)
    if [ $? -ne 0 ]; then
        echo "Error: Analysis failed for excEn${lbl}_${rec_type}" >&2
        exit 1
    fi
  ) &  # Background job
  
  job_pids+=($!)  # Capture the PID of the current background job
  jobs_running=$((jobs_running + 1))  # Increment job count
  wait_for_jobs  # Wait if we've reached the limit
}

### --------------------------------------------------------------------------- ###
### ========== Main script execution - executed jobs run in parallel ========== ###

# Create output directory
mkdir -p ../${reaction}_sim/Det_analysis

# Precompile the analysis library
compile_library_if_needed "libDetAnalysis.so" \
  "UtilityScripts/det_analysis_impl.cxx" "det_analysis.cxx"
echo ""
echo "Starting analysis loop over HR modes..."

# Loop over each HR mode using stored information
for i in "${!mode_names[@]}"; do
    mode_clean=${mode_names[$i]}
    start_idx=${mode_start_indices[$i]}
    stop_idx=${mode_stop_indices[$i]}
    if [ "$start_idx" -le "$stop_idx" ] && [ "$stop_idx" -ge 0 ]; then
        for j in $(seq $start_idx $stop_idx); do
            run_analysis $j "$mode_clean"
        done
    fi
done

# Wait for all remaining jobs to finish and update progress
wait_for_jobs 0  # Wait until 0 jobs running (i.e., all done)
# Note: wait_for_jobs will exit the script immediately if any job fails

echo ""  # New line after progress bar completes

### Now that analysis is complete, run the plotting script ###
echo ""
echo "Analysis complete. Running plotting script automatically..."

### --------------------------------------------------------------------------- ###
### ================ Default plots are generated automatically ================ ###

# Compile plotting library (links against libDetAnalysis.so)
compile_library_if_needed "libMakePlots.so" \
  "UtilityScripts/make_plots_impl.cxx" \
  "UtilityScripts/banana_plots_impl.cxx" \
  "make_plots.cxx" \
  -L. -lDetAnalysis

# Load compiled library and run (cd to script directory to ensure libraries are found)
echo "Running plotting script..."
root -l -b <<ROOT_EOF
gSystem->Load("libDetAnalysis.so");
gSystem->Load("libMakePlots.so");
.L UtilityScripts/make_plots.h
make_plots();
.q
ROOT_EOF

# Cleanup compiled libraries (optional - comment out to keep for faster reruns)
# rm libDetAnalysis.so
# rm libMakePlots.so
# echo "Cleaned up compiled libraries"