#!/bin/bash

# Loop over begin times and reverse options
for begin_time in {0..20000..2000}; do
    for reverse in True False; do
    	for sc in 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5; do
			# Determine the output directory and file
			folder="/home/fluffymelon/PANDA_exp/scaling/production/calcite_decane_tip4p_sc_${sc}_offset_0.2"
			direction=$( [ "$reverse" == "True" ] && echo "reverse" || echo "forward" )
        	output_dir="${folder}/stat_ineff/$direction/$begin_time"
        	mv "$output_dir" "${output_dir}_backup"
			mkdir -p "$output_dir"
        	output_file="$output_dir/output.txt"

        	# Run the Python script with the specified arguments
        	python3 ~/PANDA/src/utils_py/stat_ineff.py \
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
            --display False &
		done
    done
done
