#!/bin/bash

# Output file to store timing results
output_file="timing_results.txt"

# Clear the output file if it exists
> "$output_file"

# Directory containing input files
input_dir="input_data"

# Loop through the specific input files by expanding the brace expression manually
for input_file in "$input_dir"/ellipse_N_03000.gal "$input_dir"/ellipse_N_04000.gal "$input_dir"/ellipse_N_05000.gal "$input_dir"/ellipse_N_06000.gal "$input_dir"/ellipse_N_07000.gal "$input_dir"/ellipse_N_08000.gal "$input_dir"/ellipse_N_09000.gal "$input_dir"/ellipse_N_10000.gal; do
    # Extract the number of particles (N) from the filename
    N=$(basename "$input_file" | grep -oE '[0-9]+')

    # Loop through the number of threads (1 to 8)
    for threads in {7..14}; do
        echo "Running simulation for N = $N, Threads = $threads, Input file = $input_file"

        # Prepare the command to run
        command="./A4/OpenMP/galsim $N $input_file 100 1e-5 0 $threads"

        # Print the command to check it
        #echo "Executing command: $command"

        # Run the command with gtime and capture the output
        gtime_output=$(gtime -v $command 2>&1)

        # Extract the elapsed time from the gtime output
        elapsed_time=$(echo "$gtime_output" | grep 'Elapsed (wall clock) time (h:mm:ss or m:ss):' | awk '{print $8}')

        # Append the results to the output file
        echo "N = $N, Threads = $threads, Input file = $input_file, Elapsed time = $elapsed_time" >> "$output_file"
    done
done

echo "Timing results saved to $output_file"