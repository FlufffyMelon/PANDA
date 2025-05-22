# PANDA – Predicting Angle from Nanoscale Density Analysis

## Overview

This repository provides tools for calculating the contact angle of a droplet using a one-dimensional density profile.

## Requirements

### Input Files

To perform the contact angle calculations, you need to create a virtual Python environment. This can be done by running the following commands in your terminal:

```bash
python -m venv ./.venv  # Create virtual environment
source .venv/bin/activate  # Activate virtual environment
```

You'll also need to install the dependencies listed in `requirements.txt` using pip:
```bash
pip install -r requirements.txt
```

### Building the System

The C++ code for building the system is located in the `src/utils_cpp` directory. To compile the code, run the following command:

```bash
make  # Compile C++ code
```

After compiling the code, you can execute python script to build the system:

```bash
.venv/bin/python py/cal_script.py \
--H 9 \
--phi 0.5 \
--build false \
--Lx 20 \
--Ly 5 \
--Lz 4 \
--offset 0.2 \
--unitcell substrates/calcite/calcite_104_unitcell.gro \
--freeze_substr false \
--scale 2.5 \
--exp_folder calcite_decane_tip4p \
--system_name cal_dec_tip4p \
--gpu_id 0 \
--n_mpi 8 \
--init_core 0 \
--node 3 \
--server_folder PANDA_exp/example \
--ansambel nvt \
--nsteps 20000000 \
--temp 300
```

or simply execute the bash script

```bash
./scripts/run_calcite.sh  # Run bash script to build system
```

Python script has number of arguments. Here's a brief description of each argument:

*   `--H`: Specify the height of the pore (default=10)
*   `--phi`: Specify the volume fraction of non-wetting compount (in this case decane) (default=0.5)
*   `--build`: Specify whether to build the system (default=true)
*   `--Lx`, `--Ly`, `--Lz`: Specify the dimensions of the simulation box
*   `--offset`: Specify the offset value for non-wetting fraction region (default=0.2)
*   `--unitcell`: Specify the unit cell file for the substrate
*   `--freeze_substr`: Specify whether to freeze the substrate (default=false)
*   `--scale`: Specify the scaling factor for calcite-decane interaction
*   `--exp_folder`: Specify the folder where the simulation results will be saved
*   `--system_name`: Specify the name of the system being simulated
*   `--gpu_id`: Specify the GPU ID to use for the simulation (default=0)
*   `--n_mpi`, `--init_core`, `--node`: Specify the number of MPI processes, initial core, and node id for the simulation
*   `--server_folder`: Specify the folder on the server, where configuration folder will be copy
*   `--ansambel`: Specify the ansamble type for the simulation (default=nvt)
*   `--nsteps`: Specify the total number of steps in the simulation (default=1000000)

### Running the Simulation

After running the bash script, you need to execute the sbatch script from the directory where the system was saved **on the server**:

```bash
sbatch run.sh  # Run sbatch script to execute simulation
```

### Analyzing Results

The results of the simulation can be analyzed using the Jupyter Notebook `example.ipynb` in the `examples` folder.

# panda

A scientific package for ... (describe your package here)

## Installation

You can install directly from git:

```
pip install git+<link to this repo>
```

## Usage

```python
import panda
# use panda functions
```

## Structure
- All code is in the `panda` folder.

## License
MIT

