#!/bin/bash

# Define the folder containing the input files
folder="/home/fluffymelon/PANDA_exp/scaling/production/calcite_decane_tip4p_sc_1.5_offset_0.2"

# Loop over begin times and reverse options
for begin_time in {0..20000..2000}; do
    for reverse in True False; do
        # Determine the output directory and file
        direction=$( [ "$reverse" == "True" ] && echo "reverse" || echo "forward" )
        output_dir="stat_ineff/$direction/$begin_time"
        mkdir -p "$output_dir"
        output_file="$output_dir/output.txt"

        # Run the Python script with the specified arguments
        python3 stat_ineff.py \
            --trajectory_file "$folder/cal_dec_tip4p.xtc" \
            --topology_file "$folder/cal_dec_tip4p.gro" \
            --output_file "$output_file" \
            --rho_bulk 30.896 \
            --l 2.000444 \
            --phi 0.5 \
            --H 9 \
            --interface_type "roll" \
            --residue "DECAN" \
            --sl 200 \
            --max_block_length 5000 \
            --begin_time "$begin_time" \
            --time 40000 \
            --timestep 2 \
            --chunk_length 1000 \
            --units "ps" \
            --reverse "$reverse" \
            --display False

    done
done
