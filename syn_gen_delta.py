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
parser.add_argument("--H", type=float, required=True)
parser.add_argument("--phi", type=float, required=True)
parser.add_argument("--delta_list", type=str2list, required=True)
parser.add_argument("--interface_type", type=str, required=True)
parser.add_argument("--extention", type=str, required=True)
parser.add_argument("--folder", type=str, required=True)
parser.add_argument("--iterations", type=int, required=True)

args = parser.parse_args()

# Dynamically import the module
module = importlib.import_module(f"src.utils_py.geom.{args.interface_type.title()}")
Interface = getattr(module, args.interface_type.title())

l = args.WIDTH_X / args.H

# Generating of synthetic data
for d in args.delta_list:
    for i in range(args.iterations):
        print(f"Creating d:{d} iter:{i}")
        box = np.array([args.WIDTH_X, args.WIDTH_Y, args.H])

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
            delta=d,
            extention=args.extention,
        )

        insert_shapes = [region]
        shapes = [region]
        numbers = list(
            np.round(
                np.array(
                    [shapes[i].get_volume() * density[i] for i in range(len(names))]
                )
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

        dir_path = os.path.join(args.folder, f"{d:.2f}", str(args.interface_type))
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

        with open(os.path.join(dir_path, f"{i}.gro"), "w") as f:
            f.write(write_gro(structure))
