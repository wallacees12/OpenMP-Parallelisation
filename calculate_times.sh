#!/bin/bash

# Define parameters
input_1000="input_data/ellipse_N_01000.gal"
input_3000="input_data/ellipse_N_03000.gal"
output_dir="OpenMP/result.gal"
ref_output_3000="ref_output_data/ellipse_N_03000_after100steps.gal"
ref_output_1000="ref_output_data/ellipse_N_01000_after200steps.gal"

# Define output file
final_output="benchmark_results.txt"

# Clear previous output file
> $final_output

# Function to extract elapsed time from gtime output
extract_time() {
    grep "Elapsed (wall clock) time" | awk -F ': ' '{print $2}'
}

# Function to extract max difference from compare output
extract_max_diff() {
    grep "pos_maxdiff" | awk -F '=' '{print $2}' | tr -d ' '  # Extracts the number after '=' and removes spaces
}

# Run benchmarks for N=1000 (NO compare needed)
echo "Running OpenMP Galsim for N=1000..."
for threads in {1..11}; do
    cmd="./OpenMP/galsim 1000 $input_1000 200 0.00001 0 $threads"
    echo "Running: gtime -v $cmd" | tee -a $final_output
    time_result=$(gtime -v $cmd 2>&1 | extract_time)

    # Run compare to calculate the error
    compare_cmd="./compare_gal_files/compare 1000 $output_dir $ref_output_1000"
    echo "Running: $compare_cmd" | tee -a $final_output
    error_result=$($compare_cmd 2>&1 | extract_max_diff)

    # Handle case where no error value is found
    if [[ -z "$error_result" ]]; then
        error_result="N/A"
    fi

    # Print formatted result
    echo "$cmd | $time_result | $error_result" | tee -a $final_output
    echo "" >> $final_output
done

# Run benchmarks for N=3000 and calculate errors
echo "Running OpenMP Galsim for N=3000..."
for threads in {1..19}; do
    cmd="./OpenMP/galsim 3000 $input_3000 200 0.00001 0 $threads"
    echo "Running: gtime -v $cmd" | tee -a $final_output
    time_result=$(gtime -v $cmd 2>&1 | extract_time)

    # Run compare to calculate the error
    compare_cmd="./compare_gal_files/compare 3000 $output_dir $ref_output_3000"
    echo "Running: $compare_cmd" | tee -a $final_output
    error_result=$($compare_cmd 2>&1 | extract_max_diff)

    # Handle case where no error value is found
    if [[ -z "$error_result" ]]; then
        error_result="N/A"
    fi

    # Print formatted result
    echo "$cmd | $time_result | $error_result" | tee -a $final_output
    echo "" >> $final_output
done

echo "Benchmarking and error calculations completed!"