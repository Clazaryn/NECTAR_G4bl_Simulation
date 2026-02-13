#!/usr/bin/env python3.5
"""
Estimate remaining time for recoil simulations by monitoring G4bl output files.
Reads the output file every second and calculates time estimates based on event completion.
"""

import sys
import time
import os
import re
from datetime import datetime, timedelta

def is_simulation_complete(output_file):
    """
    Check if the G4bl simulation has completed by looking for completion indicators.
    Returns True if simulation appears to be finished, False otherwise.
    """
    # Convert to absolute path to handle relative paths correctly when run in background
    if not os.path.isabs(output_file):
        output_file = os.path.abspath(output_file)
    
    if not os.path.exists(output_file):
        return False
    
    try:
        with open(output_file, 'r') as f:
            content = f.read()
        
        # Look for completion indicators in G4bl output
        # Common patterns: "Simulation completed", "Finished", "Total events", etc.
        completion_patterns = [
            r'Simulation\s+completed',
            r'Simulation\s+finished',
            r'Finished\s+simulation',
            r'Total\s+events\s+processed',
            r'All\s+events\s+completed',
            r'Run\s+completed',
        ]
        
        for pattern in completion_patterns:
            if re.search(pattern, content, re.IGNORECASE):
                return True
        
        return False
    except Exception as e:
        return False

def count_events_in_output(output_file):
    """
    Count completed events in G4bl output file.
    G4bl typically outputs progress like "Event N completed" or similar.
    We'll look for the highest event number mentioned.
    """
    # Convert to absolute path to handle relative paths correctly when run in background
    if not os.path.isabs(output_file):
        output_file = os.path.abspath(output_file)
    
    if not os.path.exists(output_file):
        return 0
    
    try:
        with open(output_file, 'r') as f:
            content = f.read()
        
        # Look for event completion patterns in G4bl output
        # Common patterns: "Event N", "Processing event N", "Event N completed", etc.
        event_patterns = [
            r'Event\s+(\d+)\s+completed',
            r'Processing\s+event\s+(\d+)',
            r'Event\s+(\d+)\s+of',
            r'Event\s+(\d+)[\s:]',
            r'event\s+(\d+)\s+completed',
            r'Completed\s+event\s+(\d+)',
        ]
        
        max_event = 0
        for pattern in event_patterns:
            matches = re.findall(pattern, content, re.IGNORECASE)
            for match in matches:
                try:
                    event_num = int(match)
                    if event_num > max_event:
                        max_event = event_num
                except ValueError:
                    continue
        
        # If no pattern matches, try a simpler approach: look for "Event" followed by numbers
        if max_event == 0:
            # Look for "Event" or "event" followed by a number
            simple_pattern = r'(?:Event|event)[\s:]+(\d+)'
            matches = re.findall(simple_pattern, content, re.IGNORECASE)
            for match in matches:
                try:
                    event_num = int(match)
                    if event_num > max_event:
                        max_event = event_num
                except ValueError:
                    continue
        
        # Return max_event + 1 if we found events (since events are 0-indexed or 1-indexed)
        # But actually, if max_event is 5, we've completed 6 events (0-5) or events 1-5
        # To be safe, return max_event as the count (assuming 1-indexed)
        # If max_event is 0, return 0
        return max_event if max_event > 0 else 0
    except Exception as e:
        return 0

def count_total_events(input_file):
    """Count total number of events in the input event file."""
    # Convert to absolute path to handle relative paths correctly when run in background
    if not os.path.isabs(input_file):
        input_file = os.path.abspath(input_file)
    
    if not os.path.exists(input_file):
        return 0
    
    try:
        with open(input_file, 'r') as f:
            lines = f.readlines()
        
        # Skip header lines (lines starting with #)
        event_lines = [line for line in lines if not line.strip().startswith('#') and line.strip()]
        return len(event_lines)
    except Exception as e:
        return 0

def format_time(seconds):
    """Format seconds into hours, minutes, seconds."""
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    secs = int(seconds % 60)
    return hours, minutes, secs

def estimate_recoil_timing(output_file, input_file, total_iterations, N_parallel):
    """
    Monitor G4bl output and estimate remaining time.
    
    Args:
        output_file: Path to G4bl output file (e.g., g4bl_238U_dp_HRg_excEn00MeV_recoil.out)
        input_file: Path to input event file (e.g., ./238U_dp_results/Event_output/...recoil.txt)
        total_iterations: Total number of iterations (excitation energies) remaining
        N_parallel: Number of parallel jobs
    """
    # Convert input file to absolute path for reliable access when run in background
    if not os.path.isabs(input_file):
        input_file = os.path.abspath(input_file)
    
    # Count total events needed
    total_events = count_total_events(input_file)
    
    # If less than 100 events, skip timing (since 100 recoils take ~17s)
    if total_events < 100:
        # Debug output to help diagnose issues
        if total_events == 0:
            print("Warning: Could not read events from file: {}".format(input_file), file=sys.stderr, flush=True)
            print("File exists: {}".format(os.path.exists(input_file)), file=sys.stderr, flush=True)
        print("Timing suppressed (input file has < 100 events, found: {})".format(total_events), file=sys.stderr, flush=True)
        return
    
    # Print initial message to stderr so it's visible even when run in background
    # Use newline to avoid overwriting progress bar
    print("Starting timing estimation (monitoring {} events, {} iterations remaining)...".format(total_events, total_iterations), file=sys.stderr, flush=True)
    
    start_time = time.time()
    last_update_time = start_time
    update_interval = 10  # Update every 10 seconds
    max_duration = 90  # Run for 90 seconds max
    first_update = True  # Flag to ensure first update is printed
    last_events_completed = 0  # Track previous event count to only update when it changes
    
    while True:
        current_time = time.time()
        elapsed = current_time - start_time
        
        # Check if simulation has completed - if so, exit early
        if is_simulation_complete(output_file):
            print("Timing estimation completed (simulation finished)".format(int(elapsed)), file=sys.stderr, flush=True)
            break
        
        # Stop after max_duration
        if elapsed >= max_duration:
            print("Timing estimation completed (monitored for {} seconds)".format(int(elapsed)), file=sys.stderr, flush=True)
            break
        
        # Count events completed
        events_completed = count_events_in_output(output_file)
        
        # Only update display if events have changed AND it's time for an update
        events_changed = (events_completed != last_events_completed)
        time_for_update = first_update or (current_time - last_update_time >= update_interval)
        
        if events_changed and time_for_update:
            first_update = False
            last_events_completed = events_completed  # Update tracked count
            if events_completed > 0:
                # Calculate time per event
                time_per_event = elapsed / events_completed
                
                # Estimate time for all remaining iterations
                # Account for parallelism: remaining iterations will run in batches of N
                if total_iterations > 0:
                    # Time per iteration = time per event * total events
                    time_per_iteration = time_per_event * total_events
                    # Total time for remaining iterations (accounting for parallelism)
                    estimated_remaining_seconds = (time_per_iteration * total_iterations) / N_parallel
                    
                    # Calculate predicted finish time
                    finish_time = datetime.now() + timedelta(seconds=estimated_remaining_seconds)
                    finish_time_str = finish_time.strftime("%Y-%m-%d %H:%M:%S")
                    
                    # Print on new line to avoid overwriting progress bar
                    # Use newline so it appears below the progress bar
                    print("Predicted finish time: {} "
                          "(Event {}/{} completed, "
                          "{} iterations remaining)".format(finish_time_str, events_completed, total_events, total_iterations), file=sys.stderr, flush=True)
                else:
                    # No iterations remaining, just show current progress
                    print("Event {}/{} completed".format(events_completed, total_events), file=sys.stderr, flush=True)
            else:
                # No events completed yet
                print("Waiting for events to complete... (elapsed: {}s)".format(int(elapsed)), file=sys.stderr, flush=True)
            
            last_update_time = current_time
        elif events_changed:
            # Events changed but not time for update yet - just update the tracked count
            last_events_completed = events_completed
        
        # Sleep for 1 second before next check
        time.sleep(1)
    
    print("", file=sys.stderr)  # Final newline

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: estimate_recoil_timing.py <output_file> <input_file> <total_iterations> <N_parallel>", file=sys.stderr)
        sys.exit(1)
    
    try:
        output_file = sys.argv[1]
        input_file = sys.argv[2]
        total_iterations = int(sys.argv[3])
        N_parallel = int(sys.argv[4])
        
        # Validate that input_file is not empty
        if not input_file or input_file.strip() == "":
            print("Error: input_file argument is empty", file=sys.stderr)
            sys.exit(1)
        
        # Validate that output_file is not empty
        if not output_file or output_file.strip() == "":
            print("Error: output_file argument is empty", file=sys.stderr)
            sys.exit(1)
        
        estimate_recoil_timing(output_file, input_file, total_iterations, N_parallel)
    except Exception as e:
        # Print errors to stderr so they're visible even if stdout is redirected
        print("Error in estimate_recoil_timing: {}".format(e), file=sys.stderr)
        import traceback
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)
