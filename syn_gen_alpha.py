import argparse
import importlib
import os

import numpy as np

from src.utils_py.assembler.build import build_system
from src.utils_py.gro.Structure import Structure
from src.utils_py.io.gro import write_gro

""" Don't change!!! """
insertion_limit = int(1e5)
rotation_limit = 1000
package = 0.3
distance = {"min": 0.08**2, "opt": 0.12**2}


def str2list(arg):
    return list(map(float, arg.split(",")))


parser = argparse.ArgumentParser()
parser.add_argument("--WIDTH_X", type=float, required=True)
parser.add_argument("--WIDTH_Y", type=float, required=True)
parser.add_argument("--l", type=float, default=None, required=False)
parser.add_argument("--H", type=float, default=None, required=False)
parser.add_argument("--phi", type=float, required=True)
parser.add_argument("--theta", type=float, required=True)
parser.add_argument("--delta", type=float, required=True)
parser.add_argument("--interface_type", type=str, required=True)
parser.add_argument(
    "--extention", type=str, choices=["theta", "delta", "alpha"], required=True
)
parser.add_argument("--folder", type=str, required=True)
parser.add_argument("--iteration", type=int, required=True)

args = parser.parse_args()

assert args.l is not None or args.H is not None, "Provide eiter l, or H"

# Dynamically import the module
module = importlib.import_module(
    f"src.utils_py.geom.interface_geom"
)
Interface = getattr(module, args.interface_type.title())

if args.l is None:
    l = args.WIDTH_X / args.H
    H = args.H

if args.H is None:
    l = args.l
    H = args.WIDTH_X / l

theta = np.deg2rad(args.theta)
delta = args.delta

# Generating of synthetic data
box = np.array([args.WIDTH_X, args.WIDTH_Y, H])

structure = Structure(
    title="Points",
    box=box,
    atoms=np.empty(0, dtype=object),
    atoms_xyz=np.zeros((0, 3)),
)

points = structure.atoms_xyz

names = ["point"]
density = [12]  # nm-3 Water like density

region = Interface(
    box / 2,
    box,
    l,
    args.phi,
    theta=theta,
    delta=delta,
    extention=args.extention,
)

if not region.check_existence():
    print("Shape can not exist")
    exit()

if args.interface_type.lower() in ["droplet", "worm", "antidroplet", "antiworm"]:
    region.center[2] -= box[2] / 2
    region.center[2] -= region.rd * np.cos(theta) * H
    region.center[2] += delta * box[2]

insert_shapes = [region]
shapes = [region]
numbers = list(
    np.round(
        np.array([shapes[i].get_volume() * density[i] for i in range(len(names))])
    ).astype(int)
)

structure = build_system(
    "./",
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

dir_path = os.path.join(
    args.folder,
    args.interface_type.lower(),
    f"iter_{args.iteration}",
)
if not os.path.exists(dir_path):
    os.makedirs(dir_path, exist_ok=True)

# with open(os.path.join(dir_path, f"{i}.gro"), "w") as f:
with open(
    os.path.join(
        dir_path,
        f"{args.interface_type.lower()}_{l:.2f}_{args.phi:.2f}_{args.theta:.2f}_{delta:.2f}.gro",
    ),
    "w",
) as f:
    f.write(write_gro(structure))
