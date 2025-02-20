#!/bin/bash

# Define parameters
input_1000="input_data/ellipse_N_01000.gal"
input_3000="input_data/ellipse_N_03000.gal"

# Define reference output files for normal tests (steps=200)
ref_output_1000="ref_output_data/ellipse_N_01000_after200steps.gal"
ref_output_3000="ref_output_data/ellipse_N_03000_after100steps.gal"

# Define reference output file for final test (steps=400)
ref_output_3000_final="ref_output_data/ellipse_N_03000_after200steps.gal"

# Define variants (directories) to test
variants=("Openmp" "Vanilla" "Pthreads")

# Define output file and clear previous content
final_output="benchmark_results.txt"
> "$final_output"

# Function to extract elapsed time from gtime output
extract_time() {
    grep "Elapsed (wall clock) time" | awk -F ': ' '{print $2}'
}

# Function to extract max difference from compare output
extract_max_diff() {
    grep "pos_maxdiff" | awk -F '=' '{print $2}' | tr -d ' '
}

#############################################
# Normal Tests: Using 200 steps
#############################################
for variant in "${variants[@]}"; do
    output_dir="${variant}/result.gal"
    echo "Variant: $variant - Normal Test (steps=200)" >> "$final_output"
    
    # For N=1000
    if [ "$variant" = "Vanilla" ]; then
        cmd="./${variant}/galsim 1000 $input_1000 200 0.00001 0"
        time_result=$(gtime -v $cmd 2>&1 | extract_time)
        # Print executed command and time to console
        echo "Executed: $cmd | Time: $time_result"
        compare_cmd="./compare_gal_files/compare 1000 $output_dir $ref_output_1000"
        error_result=$($compare_cmd 2>&1 | extract_max_diff)
        [[ -z "$error_result" ]] && error_result="N/A"
        echo "$cmd | $time_result | $error_result" >> "$final_output"
    else
        for threads in {1..14}; do
            cmd="./${variant}/galsim 1000 $input_1000 200 0.00001 0 $threads"
            time_result=$(gtime -v $cmd 2>&1 | extract_time)
            # Print executed command and time to console
            echo "Executed: $cmd | Time: $time_result"
            compare_cmd="./compare_gal_files/compare 1000 $output_dir $ref_output_1000"
            error_result=$($compare_cmd 2>&1 | extract_max_diff)
            [[ -z "$error_result" ]] && error_result="N/A"
            echo "$cmd | $time_result | $error_result" >> "$final_output"
        done
    fi

    # For N=3000
    if [ "$variant" = "Vanilla" ]; then
        cmd="./${variant}/galsim 3000 $input_3000 100 0.00001 0"
        time_result=$(gtime -v $cmd 2>&1 | extract_time)
        echo "Executed: $cmd | Time: $time_result"
        compare_cmd="./compare_gal_files/compare 3000 $output_dir $ref_output_3000"
        error_result=$($compare_cmd 2>&1 | extract_max_diff)
        [[ -z "$error_result" ]] && error_result="N/A"
        echo "$cmd | $time_result | $error_result" >> "$final_output"
    else
        for threads in {1..14}; do
            cmd="./${variant}/galsim 3000 $input_3000 100 0.00001 0 $threads"
            time_result=$(gtime -v $cmd 2>&1 | extract_time)
            echo "Executed: $cmd | Time: $time_result"
            compare_cmd="./compare_gal_files/compare 3000 $output_dir $ref_output_3000"
            error_result=$($compare_cmd 2>&1 | extract_max_diff)
            [[ -z "$error_result" ]] && error_result="N/A"
            echo "$cmd | $time_result | $error_result" >> "$final_output"
        done
    fi
done

#############################################
# Final Test: Using 400 steps for N=3000 only
#############################################
echo "Final Test: steps=400 (N=3000)" >> "$final_output"
for variant in "${variants[@]}"; do
    output_dir="${variant}/result.gal"
    echo "Variant: $variant - Final Test (steps=400)" >> "$final_output"
    if [ "$variant" = "Vanilla" ]; then
        cmd="./${variant}/galsim 3000 $input_3000 400 0.00001 0"
        time_result=$(gtime -v $cmd 2>&1 | extract_time)
        echo "Executed: $cmd | Time: $time_result"
        compare_cmd="./compare_gal_files/compare 3000 $output_dir $ref_output_3000_final"
        error_result=$($compare_cmd 2>&1 | extract_max_diff)
        [[ -z "$error_result" ]] && error_result="N/A"
        echo "$cmd | $time_result | $error_result" >> "$final_output"
    else
        for threads in {1..14}; do
            cmd="./${variant}/galsim 3000 $input_3000 400 0.00001 0 $threads"
            time_result=$(gtime -v $cmd 2>&1 | extract_time)
            echo "Executed: $cmd | Time: $time_result"
            compare_cmd="./compare_gal_files/compare 3000 $output_dir $ref_output_3000_final"
            error_result=$($compare_cmd 2>&1 | extract_max_diff)
            [[ -z "$error_result" ]] && error_result="N/A"
            echo "$cmd | $time_result | $error_result" >> "$final_output"
        done
    fi
done

echo "Benchmarking and error calculations completed!"