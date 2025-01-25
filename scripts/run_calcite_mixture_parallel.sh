# List of parameter triples (l, H, phi) as tuples
# pairs=("1 11.6 0.05" "1 11.6 0.15" "1 11.6 0.25" "1 11.6 0.35" \
#        "2 5.8 0.05" "2 5.8 0.15" "2 5.8 0.25" "2 5.8 0.35" "2 5.8 0.45" "2 5.8 0.55" "2 5.8 0.65"
#        "3 3.87 0.05" "3 3.87 0.15" "3 3.87 0.25" "3 3.87 0.35" "3 3.87 0.45" "3 3.87 0.55" "3 3.87 0.65" "3 3.87 0.75" "3 3.87 0.85" \
#        "4 2.9 0.05" "4 2.9 0.15" "4 2.9 0.25" "4 2.9 0.35" "4 2.9 0.45" "4 2.9 0.55" "4 2.9 0.65" "4 2.9 0.75" "4 2.9 0.85" "4 2.9 0.95" "4 2.9 1.05")

pairs=("2 5.8 0.35" "2 5.8 0.45" "3 3.87 0.25")


# # Set batch size
# batch_size=10

# # Calculate number of batches
# num_batches=$(( ${#pairs[@]} / batch_size ))

# # Loop over batches
# for ((i = 0; i < num_batches; i++)); do
#     # Get start and end indices for this batch
#     start=$(( i * batch_size ))
#     end=$(( (i + 1) * batch_size - 1 ))

#     if (( end >= ${#pairs[@]} )); then
#         end=$(( ${#pairs[@]} - 1 ))
#     fi

#     # Print batch number
#     echo "Batch $((i+1)) of ${num_batches}"

# Loop over tasks in this batch and execute Python script
# for pair in "${pairs[@]:start:end}"; do
for pair in "${pairs[@]}"; do
    l=$(echo "$pair" | cut -d' ' -f1)
    H=$(echo "$pair" | cut -d' ' -f2)
    phi=$(echo "$pair" | cut -d' ' -f3)

    # Create unique folder name
    exp_folder="calcite_decane_water_120_${l}_${phi}"

    # Execute Python script with corresponding arguments
    python py/cal_script_mixture.py \
        --H $H \
        --phi $phi \
        --build true \
        --Lx 12 \
        --Ly 12 \
        --Lz 4 \
        --offset 0.2 \
        --unitcell substrates/calcite/calcite_104_unitcell.gro \
        --gen_substr false \
        --freeze_substr false \
        --scale 2.3 \
        --ansambel nvt \
        --nsteps 2500000 \
        --temp 300 \
        --gpu_id 0 \
        --n_mpi 12 \
        --init_core 0 \
        --node 1 \
        --exp_folder $exp_folder \
        --system_name cal_dec_tip4p \
        --server_folder PANDA_exp/stability/angle_120_iter_2 \
        --send_to_server false &

done

#     # Wait for all jobs in this batch to finish
#     wait

#     # Print new line to separate batches
#     echo ""
# done
