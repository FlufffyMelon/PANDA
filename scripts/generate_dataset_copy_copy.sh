#!/bin/bash
#SBATCH --job-name=python  # Name of the job
#SBATCH -w node3
#SBATCH --ntasks=1                       # Number of tasks (1 task = 1 instance of the script)
#SBATCH --cpus-per-task=128               # Number of CPU cores per task (adjust based on your server)

# Define the shapes and their iterations
# shapes=("Droplet" "Doughnut" "Worm" "Roll" "Perforation" "Layer")
shapes=("Perforation")
# shapes_iter=(40 10 20 10 10 2)
shapes_iter=(10)

# Define the parameter values
L=(1 2 3 4)
PHI=(0.05 0.15 0.25 0.35 0.45 0.55 0.65 0.75 0.85 0.95)
# PHI=(0.15 0.45 0.95)
THETA=(91 100 110 120 130 140 150 160 170 180)
# THETA=(91 130 180)
DELTA=(0 0.05 0.1 0.15 0.2)
# DELTA=(0 0.1)

# Output file
output_file="parameters.txt"

# Clear the output file if it exists
> "$output_file"

# Function to add normally distributed noise to a single value
add_noise() {
    local point=$1  # Single point value
    local step=$(awk -v a=$2 'BEGIN {
        print a / 6
    }')           # Standard deviation of the noise
    local lower=$3  # Lower bound for generated values
    local upper=$4  # Upper bound for generated values

    # Ensure the generated value is not below the lower bound
    while true; do
        val=$(awk -v point=$point -v step=$step -v seed=$RANDOM '
            BEGIN {
                srand(seed)
                pi = 3.14159265358979323846
                u1 = rand()
                u2 = rand()
                z0 = sqrt(-2 * log(u1)) * cos(2 * pi * u2)  # Standard normal distribution
                val = point + z0 * step  # Scale by step (standard deviation)
                print val
            }
        ')

        # Check if the value is within the lower and upper bounds
        if (( $(echo "$val >= $lower && $val <= $upper" | ~/bc-1.03/bc -l) )); then
            break
        fi
    done

    # Return the noisy value
    echo "$val"
}

# Iterate over shapes
for ((i=0; i<${#shapes[@]}; i++)); do
    echo $i
    shape=${shapes[i]}
    iter=${shapes_iter[i]}

    # Iterate over the number of iterations for the current shape
    for ((j=0; j<iter; j++)); do
        # Iterate over the initial arrays and add noise to each value
        for l in ${L[@]}; do
            for phi in ${PHI[@]}; do
                for theta in ${THETA[@]}; do
                    for delta in ${DELTA[@]}; do
                        l_noisy=$(add_noise $l 1 1 10)  # Step = 1, lower bound = 1, upper bound = 10
                        # echo "l_noisy $l $l_noisy"
                        phi_noisy=$(add_noise $phi 0.1 0.05 1)  # Step = 0.1, lower bound = 0.05, upper bound = 1.0
                        # echo "phi_noisy $phi $phi_noisy"
                        theta_noisy=$(add_noise $theta 10 91 180)  # Step = 10, lower bound = 91, upper bound = 180
                        # echo "theta_noisy $theta $theta_noisy"
                        delta_noisy=$(add_noise $delta 0.05 0 0.5)  # Step = 0.05, lower bound = 0, upper bound = 0.5
                        # echo "delta_noisy $delta $delta_noisy"

                        # Write parameters to the output file, formatted to 2 decimal places
                        echo "$shape $l_noisy $phi_noisy $theta_noisy $delta_noisy"

                        .venv/bin/python syn_gen_alpha.py \
                        --WIDTH_X 10 \
                        --WIDTH_Y 10 \
                        --l $l_noisy \
                        --phi $phi_noisy \
                        --theta $theta_noisy \
                        --delta $delta_noisy \
                        --interface_type $shape \
                        --extention "alpha" \
                        --folder "data/Interface_test/" \
                        --iteration $j &
                    done
                done
                wait
            done
        done
    done
done
