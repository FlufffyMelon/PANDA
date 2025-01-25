import numpy as np
import argparse
import os
import sys

from src.utils_py.utils import str2bool

parser = argparse.ArgumentParser()
# Paths and names
parser.add_argument("--HOME_DIR", type=str, default="", help="Path to home directory of the project")
parser.add_argument("--exp_folder", type=str, help="Name of folder for saving")
parser.add_argument("--system_name", type=str, help="Name of final structure")

# System generation flags
parser.add_argument("--H", type=float, help="Slot H")
parser.add_argument("--phi", type=float, help="Amoung of water")
parser.add_argument("--build", type=str2bool, help="Build system flag")
parser.add_argument("--Lx", type=float, help="Lenght of substrate in Ox direction")
parser.add_argument("--Ly", type=float, help="Lenght of substrate in Oy direction")
parser.add_argument("--Lz", type=float, help="Lenght of substrate in Oz direction")
parser.add_argument("--offset", type=float, help="Wetting layer thickness")
parser.add_argument("--unitcell", type=str, help="Path to unitcell gro file")
parser.add_argument("--gen_substr", type=str2bool, default=True, help="Generate substrates` gro and itp")

# Gromacs files configuration
parser.add_argument("--freeze_substr", type=str2bool, default=True, help="Freeze substrate flag")
parser.add_argument("--scale", type=float, default=1,help="Scaling factor of calcite–decane pairwise LJ")
parser.add_argument("--ansambel", type=str, choices=["nvt", "npt"], help="Ansambel of simulation")
parser.add_argument("--nsteps", type=int, default=5000000, help="Number of simulation steps")
parser.add_argument("--temp", type=int, default=393, help="Target temperature, K")
parser.add_argument("--press", type=int, default=1, help="Target pressure, bar")

# Sbatch script configuration
parser.add_argument("--gpu_id", type=str, default="0", choices=["0", "1", "01"], help="Id of GPU")
parser.add_argument("--n_mpi", type=int, default=8, help="Number of MPI processes")
parser.add_argument("--init_core", type=int, default=0, help="Id of initial CPU core to bind process")
parser.add_argument("--node", type=int, default=1, help="Number of node on softcluster")

# Remote server params
parser.add_argument("--send_to_server", type=str2bool, default=True, help="Send to remote server flag")
parser.add_argument("--server_folder", type=str, help="Folder on softcluster")


# Parse arguments
args = parser.parse_args()

# Argument validation
assert (
    not args.freeze_substr or args.ansambel != "npt"
), "Can not run npt with freeze substr"
assert os.path.isfile(args.unitcell), "No such unitcell file!"
if not args.build: args.gen_substr = args.build

# sys.path.append(args.HOME_DIR)  # Avoid error with importing of src
from src.utils_py.io.gro import read_gro, write_gro  # noqa: E402
from src.utils_py.geom.Box import Box  # noqa: E402
from src.utils_py.geom.AntiBox import AntiBox  # noqa: E402
from src.utils_py.assembler.build import build_system  # noqa: E402
from src.utils_py.substr import generate_substrate, generate_calcite_itp  # noqa: E402
from src.utils_py.parser import parse_C6_C12  # noqa: E402


# Generating substrate with target dimension if needed
substr_name = generate_substrate(
    args.unitcell, args.Lx, args.Ly, args.Lz, freeze_substr=args.freeze_substr, build=args.gen_substr
)

substr_folder, _ = os.path.split(args.unitcell)
substr_path = os.path.join(substr_folder, "gro", substr_name)
structure = read_gro(substr_path)

# Generating itp for substrate if needed
if not args.freeze_substr:
    substr_itp_name = generate_calcite_itp(substr_path, build=args.gen_substr)


N = structure.atoms[-1].mol_id
WIDTH_X, WIDTH_Y, HEIGHT = structure.box
H = args.H
delta_h = args.offset
frac = args.phi

ff_type = "trappe"
ff_path = os.path.join(args.HOME_DIR, "ff", ff_type)

""" DO NOT TOUCH!!! """
insertion_limit = int(1e5)
rotation_limit = 1000
package = 0.3
distance = {"min": 0.08**2, "opt": 0.12**2}

system_size = np.array([WIDTH_X, WIDTH_Y, HEIGHT + H])
points = structure.atoms_xyz
structure.box = system_size

names = ["decane", "tip4p"]
density = [3.0896, 33.3277]  # nm-3


# Creating regions for decane and water
decane_box = Box(
    center=[WIDTH_X * (1 - frac / 2), WIDTH_Y / 2, HEIGHT + H / 2],
    borders=[WIDTH_X * frac, WIDTH_Y, H],
)

water_box = Box(
    center=[WIDTH_X * (1 - frac) / 2, WIDTH_Y / 2, HEIGHT + H / 2],
    borders=[WIDTH_X * (1 - frac), WIDTH_Y, H],
)

insertation_box = Box(
    center=[WIDTH_X / 2, WIDTH_Y / 2, HEIGHT + H / 2],
    borders=[WIDTH_X, WIDTH_Y, H - 2 * args.offset],
)

insert_shapes = [insertation_box, insertation_box]
shapes = [decane_box, water_box]
numbers = list(
    np.round(
        np.array([shapes[i].get_volume() * density[i] for i in range(len(names))])
    ).astype(int)
)


# Building system
if args.build:
    structure = build_system(
        os.path.join(ff_path, "gro"),
        structure,
        names,
        numbers,
        insert_shapes,
        points,
        insertion_limit=insertion_limit,
        rotation_limit=rotation_limit,
        package=package,
        min_dist2=distance["min"],
    )


filename = args.system_name
exp_path = os.path.join("systems", args.exp_folder)
if not os.path.exists(exp_path):
    os.makedirs(exp_path)


# Generating system.itp
with open(os.path.join(exp_path, "system.itp"), "w") as f:
    for name in names:
        f.write(f'#include "{name}.itp"\n')

    if args.freeze_substr:
        f.write('#include "cal.itp"\n')
    else:
        f.write(f'#include "{substr_itp_name}"\n')

    f.write(f"\n[ system ]\n{filename}\n")
    f.write("\n[ molecules ]\n; molecule name\tnr.\n")
    f.write(f"CAL\t{N if args.freeze_substr else 1}\n")
    for i, name in enumerate(names):
        f.write(f"{name}\t{numbers[i]}\n")

if args.build:
    # Writing structure in gro fileformat
    with open(os.path.join(exp_path, filename + ".gro"), "w") as f:
        f.write(write_gro(structure))

    print("Mixing system")
    os.system(
        f"./{os.path.join(args.HOME_DIR, 'src', 'utils_cpp', 'mixer')} -f {os.path.join(exp_path, filename + '.gro')} -o {os.path.join(exp_path, filename + '.gro')} -mn2 {distance['min']} -opt2 {distance['opt']}"
    )


# Generating run.sh
with open(os.path.join(exp_path, "run.sh"), "w") as f:
    f.write(f"""#!/bin/bash
#SBATCH -J gromacs
#SBATCH -w node{args.node}
#SBATCH -N 1\t# Number of nodes requested
#SBATCH -n {args.n_mpi}\t# Total number of mpi tasks requested
#SBATCH --cpus-per-task 1\t# Total number of omp tasks requested

gmx grompp -f nvt_cal_steep.mdp -c {filename}_init.gro -p {ff_type}.top -o {filename} -maxwarn 10
mpirun -np {args.n_mpi} --cpu-set {args.init_core}-{args.init_core+args.n_mpi-1} --bind-to core gmx_mpi mdrun -s -o -x -c -e -g -v -deffnm {filename} -ntomp 1 -nb gpu -gpu_id {args.gpu_id}
rm ./*pdb

gmx grompp -f nvt_cal_short.mdp -c {filename}.gro -p {ff_type}.top -o {filename} -maxwarn 10
mpirun -np {args.n_mpi} --cpu-set {args.init_core}-{args.init_core+args.n_mpi-1} --bind-to core gmx_mpi mdrun -s -o -x -c -e -g -v -deffnm {filename} -ntomp 1 -nb gpu -gpu_id {args.gpu_id} -dlb yes
rm ./*pdb

gmx grompp -f {args.ansambel}_cal_run.mdp -c {filename}.gro -p {ff_type}.top -o {filename} -maxwarn 10
mpirun -np {args.n_mpi} --cpu-set {args.init_core}-{args.init_core+args.n_mpi-1} --bind-to core gmx_mpi mdrun -s -o -x -c -e -g -v -deffnm {filename} -ntomp 1 -nb gpu -gpu_id {args.gpu_id} -dlb yes
rm ./*pdb""")


# Generating mdp files from templates
files = ["nvt_cal_steep.mdp", "nvt_cal_short.mdp", "nvt_cal_run.mdp"]
with open(os.path.join(ff_path, "mdp", "nvt_cal_steep.mdp")) as f:
    steep_final_text = f.read().format(
        **{
            "freeze": "freezegrps\t\t\t=  CAL\nfreezedim\t\t\t=  Y Y Y\n"
            if args.freeze_substr
            else ""
        }
    )
with open(os.path.join(exp_path, "nvt_cal_steep.mdp"), "w") as f:
    f.write(steep_final_text)

with open(os.path.join(ff_path, "mdp", "nvt_cal_short.mdp")) as f:
    short_final_text = f.read().format(
        **{
            "freeze": "freezegrps\t\t\t=  CAL\nfreezedim\t\t\t=  Y Y Y\n"
            if args.freeze_substr
            else "",
            "temp": args.temp,
        }
    )
with open(os.path.join(exp_path, "nvt_cal_short.mdp"), "w") as f:
    f.write(short_final_text)

with open(os.path.join(ff_path, "mdp", f"{args.ansambel}_cal_run.mdp")) as f:
    run_final_text = f.read().format(
        **{
            "nsteps": args.nsteps,
            "freeze": "freezegrps\t\t\t=  CAL\nfreezedim\t\t\t=  Y Y Y\n"
            if args.freeze_substr
            else "",
            "temp": args.temp,
            "press": args.press,
        }
    )
with open(os.path.join(exp_path, f"{args.ansambel}_cal_run.mdp"), "w") as f:
    f.write(run_final_text)


# Generating topoloty from template
with open(os.path.join(ff_path, f"{ff_type}.top")) as f:
    topology_text = f.read()
calcite_params = parse_C6_C12(topology_text, ["CA", "OCA", "CCA"])
decane_params = parse_C6_C12(topology_text, ["CH2", "CH3"])
pairwise_params = {}
for i, Ci in enumerate(["C6", "C12"]):
    for dec_name, Ci_dec in decane_params.items():
        for cal_name, Ci_cal in calcite_params.items():
            pairwise_params["_".join([Ci, dec_name, cal_name])] = "{:.2e}".format(
                args.scale * np.sqrt(Ci_dec[i] * Ci_cal[i])
            )
topology_text = topology_text.format(**pairwise_params)
with open(os.path.join(exp_path, f"{ff_type}.top"), "w") as f:
    f.write(topology_text)


# Sending files to the server or just copy
server_name = "softcluster"
if args.send_to_server:
    print(f"Sending the files to the {server_name}")
else:
    print(f"Copying the files on the server")

server_path = os.path.join(args.server_folder, args.exp_folder)
if args.send_to_server:
    os.system(f"ssh {server_name} 'mkdir -p {'~/'+server_path}'")
else:
    os.system(f"mkdir -p {'~/'+server_path}")

for name in names:
    if args.send_to_server:
        os.system(f"scp {os.path.join(ff_path, 'itp', str(name)+'.itp')} {server_name}:{server_path}")
    else:
        os.system(f"cp {os.path.join(ff_path, 'itp', str(name)+'.itp')} {os.path.join('~', server_path, str(name)+'.itp')}")

if args.send_to_server:
    if args.freeze_substr:
        os.system(f"scp {os.path.join(ff_path, 'itp', 'cal.itp')} {server_name}:{server_path}")
    else:
        os.system(f"scp {os.path.join(substr_folder, 'itp', substr_itp_name)} {os.path.join('~', server_path, substr_itp_name)}")
else:
    if args.freeze_substr:
        os.system(f"cp {os.path.join(ff_path, 'itp', 'cal.itp')} {server_path}")
    else:
        os.system(f"cp {os.path.join(substr_folder, 'itp', substr_itp_name)} {os.path.join('~/', server_path, substr_itp_name)}")

if args.send_to_server:
    os.system(f"scp {os.path.join(exp_path, str(ff_type)+'.top')} {server_name}:{server_path}")
else:
    os.system(f"cp {os.path.join(exp_path, str(ff_type)+'.top')} {os.path.join('~/', server_path, str(ff_type)+'.top')}")

files = ["nvt_cal_steep.mdp", "nvt_cal_short.mdp", f"{args.ansambel}_cal_run.mdp"]
for file in files:
    if args.send_to_server:
        os.system(f"scp {os.path.join(exp_path, file)} {server_name}:{server_path}")
    else:
        os.system(f"cp {os.path.join(exp_path, file)} {os.path.join('~/', server_path, file)}")

files = [filename + ".gro", "system.itp", "run.sh"]
for file in files:
    if args.send_to_server:
        os.system(f"scp {os.path.join(exp_path, file)} {server_name}:{server_path}/")
    else:
        os.system(f"cp {os.path.join(exp_path, file)} {os.path.join('~/', server_path, file)}")

if args.send_to_server:
    os.system(f"scp {os.path.join(exp_path, filename+'.gro')} {server_name}:{server_path}/{filename}_init.gro")
else:
    os.system(f"cp {os.path.join(exp_path, filename+'.gro')} {os.path.join('~/', server_path, f'{filename}_init.gro')}")

if args.send_to_server:
    print(f"The transfer to `{server_name}` was completed!")
else:
    print("Copying completed!")
