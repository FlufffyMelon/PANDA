import numpy as np
import argparse
import os
import sys


def str2bool(s: str) -> bool:
    """Helper function to support boolean command line arguments."""
    if s.lower() in {"true", "t", "1"}:
        return True
    elif s.lower() in {"false", "f", "0"}:
        return False
    else:
        raise ValueError(f"Invalid boolean value: {s}")


parser = argparse.ArgumentParser()
parser.add_argument(
    "--HOME_DIR", type=str, default="", help="Path to home directory of the project"
)
parser.add_argument("--H", type=float, help="Slot H")
parser.add_argument("--phi", type=float, help="Amoung of water")
parser.add_argument("--build", type=str2bool, help="Build system flag")
parser.add_argument("--Lx", type=float, help="Lenght of substrate in Ox direction")
parser.add_argument("--Ly", type=float, help="Lenght of substrate in Oy direction")
parser.add_argument("--Lz", type=float, help="Lenght of substrate in Oz direction")
parser.add_argument("--unitcell", type=str, help="Path to unitcell gro file")
parser.add_argument(
    "--freeze_substr", type=str2bool, default=True, help="Freeze substrate flag"
)
parser.add_argument("--exp_folder", type=str, help="Name of folder for saving")
parser.add_argument("--system_name", type=str, help="Name of final structure")
parser.add_argument("--gpu_id", type=int, default=0, help="Id of GPU")
parser.add_argument("--n_mpi", type=int, default=8, help="Number of MPI processes")
parser.add_argument(
    "--init_core", type=int, default=0, help="Id of initial CPU core to bind process"
)
parser.add_argument("--node", type=int, default=1, help="Number of node on softcluster")
parser.add_argument("--server_folder", type=str, help="Folder on softcluster")
parser.add_argument(
    "--ansambel", type=str, choices=["nvt", "npt"], help="Ansambel of simulation"
)
parser.add_argument(
    "--nsteps", type=int, default=5000000, help="Number of simulation steps"
)
parser.add_argument("--temp", type=int, default=393, help="Target temperature, K")
parser.add_argument("--press", type=int, default=1, help="Target pressure, bar")

# Parse arguments
args = parser.parse_args()

# Argument validation
assert (
    not args.freeze_substr or args.ansambel != "npt"
), "Can not run npt with freeze substr"
assert os.path.isfile(args.unitcell), "No such unitcell file!"

sys.path.append(args.HOME_DIR)  # Avoid error with importing of src
from src.utils_py.io.gro import read_gro, write_gro
from src.utils_py.geom.Box import Box
from src.utils_py.assembler.build import build_system
from src.utils_py.substr import generate_substrate, generate_calcite_itp

# Generating substrate with target dimension
substr_name = generate_substrate(
    args.unitcell, args.Lx, args.Ly, args.Lz, freeze_substr=args.freeze_substr
)

substr_folder, _ = os.path.split(args.unitcell)
substr_path = os.path.join(substr_folder, "gro", substr_name)
structure = read_gro(substr_path)

# Generating itp for substrate if needed
if not args.freeze_substr:
    substr_itp_name = generate_calcite_itp(substr_path)

N = structure.atoms[-1].mol_id
WIDTH_X, WIDTH_Y, HEIGHT = structure.box
H = args.H
delta_h = 0.2
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

names = ["decane", "tip3p"]
density = [3.0896, 33.3277]  # nm-3

# Creating regions for decane and water
box_left = Box(
    center=[WIDTH_X * (1 - frac) / 2, WIDTH_Y / 2, HEIGHT + H / 2],
    borders=[WIDTH_X * (1 - frac), WIDTH_Y, H],
)

box_left_delta = Box(
    center=[WIDTH_X * (1 - frac) / 2, WIDTH_Y / 2, HEIGHT + H / 2],
    borders=[WIDTH_X * (1 - frac), WIDTH_Y, H - 2 * delta_h],
)

box_left = Box(
    center=[WIDTH_X * (1 - frac) / 2, WIDTH_Y / 2, HEIGHT + H / 2],
    borders=[WIDTH_X * (1 - frac), WIDTH_Y, H],
)

box_left_delta = Box(
    center=[WIDTH_X * (1 - frac) / 2, WIDTH_Y / 2, HEIGHT + H / 2],
    borders=[WIDTH_X * (1 - frac), WIDTH_Y, H - 2 * delta_h],
)

box_right = Box(
    center=[WIDTH_X * (1 - frac / 2), WIDTH_Y / 2, HEIGHT + H / 2],
    borders=[WIDTH_X * frac, WIDTH_Y, H],
)

box_right_delta = Box(
    center=[WIDTH_X * (1 - frac / 2), WIDTH_Y / 2, HEIGHT + H / 2],
    borders=[WIDTH_X * frac, WIDTH_Y, H - 2 * delta_h],
)


insert_shapes = [box_left_delta, box_right_delta]
shapes = [box_left, box_right]
numbers = list(
    np.round(
        np.array([shapes[i].get_volume() * density[i] for i in range(len(names))])
    ).astype(int)
)

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


# Writing down .gro и system.itp
filename = args.system_name

exp_path = os.path.join("systems", args.exp_folder)
if not os.path.exists(exp_path):
    os.makedirs(exp_path)

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

print("Writing .gro files.")

if args.build:
    with open(os.path.join(exp_path, filename + ".gro"), "w") as f:
        f.write(write_gro(structure))

    print("Mixing system")
    os.system(
        f"./{os.path.join(args.HOME_DIR, 'src', 'utils_cpp', 'mixer')} -f {os.path.join(exp_path, filename + '.gro')} -o {os.path.join(exp_path, filename + '.gro')} -mn2 {distance['min']} -opt2 {distance['opt']}"
    )

with open(os.path.join(exp_path, "run.sh"), "w") as f:
    f.write(f"""#!/bin/bash
#SBATCH -J gromacs
#SBATCH -w node{args.node}
#SBATCH -N 1\t# Number of nodes requested
#SBATCH -n {args.n_mpi}\t# Total number of mpi tasks requested
#SBATCH --cpus-per-task 1\t# Total number of omp tasks requested

gmx grompp -f nvt_cal_steep.mdp -c {filename}_init.gro -p {ff_type}.top -o {filename} -maxwarn 10
mpirun -np {args.n_mpi} --cpu-set {args.init_core}-{args.init_core+args.n_mpi-1} --bind-to core gmx_mpi mdrun -s -o -x -c -e -g -v -deffnm {filename} -ntomp 1 -nb gpu -gpu_id {args.gpu_id}
rm *pdb

gmx grompp -f nvt_cal_short.mdp -c {filename}.gro -p {ff_type}.top -o {filename} -maxwarn 10
mpirun -np {args.n_mpi} --cpu-set {args.init_core}-{args.init_core+args.n_mpi-1} --bind-to core gmx_mpi mdrun -s -o -x -c -e -g -v -deffnm {filename} -ntomp 1 -nb gpu -gpu_id {args.gpu_id} -dlb yes
rm *pdb

gmx grompp -f {args.ansambel}_cal_run.mdp -c {filename}.gro -p {ff_type}.top -o {filename} -maxwarn 10
mpirun -np {args.n_mpi} --cpu-set {args.init_core}-{args.init_core+args.n_mpi-1} --bind-to core gmx_mpi mdrun -s -o -x -c -e -g -v -deffnm {filename} -ntomp 1 -nb gpu -gpu_id {args.gpu_id} -dlb yes
rm *pdb""")

# Generating mdp from templates
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


# Sending files to the server
server_name = "softcluster"

server_path = os.path.join(args.server_folder, args.exp_folder)
files = list((server_path).split("/"))
for i in range(len(files)):
    os.system(f"ssh {server_name} 'mkdir {'~/'+'/'.join(files[:i+1])}'")

for name in names:
    os.system(
        f"scp {os.path.join(ff_path, 'itp', str(name)+'.itp')} {server_name}:{server_path}"
    )

if args.freeze_substr:
    os.system(
        f"scp {os.path.join(ff_path, 'itp', 'cal.itp')} {server_name}:{server_path}"
    )
else:
    os.system(
        f"scp {os.path.join(substr_folder, 'itp', substr_itp_name)} {server_name}:{server_path}"
    )

os.system(
    f"scp {os.path.join(ff_path, str(ff_type)+'.top')} {server_name}:{server_path}"
)

files = ["nvt_cal_steep.mdp", "nvt_cal_short.mdp", f"{args.ansambel}_cal_run.mdp"]
for file in files:
    os.system(f"scp {os.path.join(exp_path, file)} {server_name}:{server_path}")

files = [filename + ".gro", "system.itp", "run.sh"]
for file in files:
    os.system(f"scp {os.path.join(exp_path, file)} {server_name}:{server_path}")

os.system(
    f"scp {os.path.join(exp_path, filename+'.gro')} {server_name}:{server_path}/{filename}_init.gro"
)
