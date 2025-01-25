# List of parameter triples (l, H, phi) as tuples
pairs=("1 11.6 0.05" "1 11.6 0.15" "1 11.6 0.25" "1 11.6 0.35" "1 11.6 0.45" "1 11.6 0.55" "1 11.6 0.65" "1 11.6 0.75" "1 11.6 0.85" \
       "2 5.8 0.05" "2 5.8 0.15" "2 5.8 0.25" "2 5.8 0.35" "2 5.8 0.45" "2 5.8 0.55" "2 5.8 0.65" "2 5.8 0.75" "2 5.8 0.85" \
       "3 3.87 0.05" "3 3.87 0.15" "3 3.87 0.25" "3 3.87 0.35" "3 3.87 0.45" "3 3.87 0.55" "3 3.87 0.65" "3 3.87 0.75" "3 3.87 0.85" \
       "4 2.9 0.05" "4 2.9 0.15" "4 2.9 0.25" "4 2.9 0.35" "4 2.9 0.45" "4 2.9 0.55" "4 2.9 0.65" "4 2.9 0.75" "4 2.9 0.85")

node=2
init_core=0
gpu_id="0"
init_dependency=10245
dependency=-1
counter=0

for pair in "${pairs[@]}"; do
    l=$(echo "$pair" | cut -d' ' -f1)
    H=$(echo "$pair" | cut -d' ' -f2)
    phi=$(echo "$pair" | cut -d' ' -f3)

    # Create unique folder name
    exp_folder="calcite_decane_water_120_${l}_${phi}"

    echo "Node: $node, Init Core: $init_core, GPU ID: $gpu_id, Dependency: $dependency"

    # Execute Python script with corresponding arguments
    python py/cal_script_cube.py \
        --H $H \
        --phi $phi \
        --build false \
        --Lx 12 \
        --Ly 12 \
        --Lz 4 \
        --offset 1.0 \
        --unitcell substrates/calcite/calcite_104_unitcell.gro \
        --gen_substr false \
        --freeze_substr false \
        --scale 2.3 \
        --ansambel nvt \
        --nsteps 2500000 \
        --temp 300 \
        --gpu_id $gpu_id \
        --n_mpi 12 \
        --init_core $init_core \
        --node $node \
        --dependency $dependency \
        --exp_folder $exp_folder \
        --system_name cal_dec_tip4p \
        --server_folder PANDA_exp/stability/cube/angle_120_iter_3_offset_1.0 \
        --send_to_server false &

    init_core=$(( init_core + 12 ))
    if [ $init_core -eq 12 ]; then
        gpu_id="1"
    elif [ $init_core -eq 24 ]; then
        gpu_id="01"
    elif [ $init_core -eq 36 ]; then
        init_core=0
        gpu_id="0"
        # if [ $node -eq 1 ]; then
        #     node=2
        # elif [ $node -eq 2 ]; then
        #     node=1
        # fi
    fi

    counter=$(( counter + 1 ))
    if [ $((counter)) -eq 3 ]; then
        dependency=$init_dependency
    fi

    if [ $((counter)) -gt 3 ]; then
        dependency=$(( dependency + 1 ))
    fi
done
