import os

import numpy as np
import argparse

from src.utils_py.assembler.build import build_system
from utils_py.geom.interface_geom.Droplet import Droplet
from utils_py.geom.interface_geom.Doughnut import Doughnut
from utils_py.geom.interface_geom.Worm import Worm
from utils_py.geom.interface_geom.Roll import Roll
from utils_py.geom.interface_geom.Perforation import Perforation
from src.utils_py.gro.Structure import Structure
from src.utils_py.io.gro import read_gro, write_gro

""" Don't change!!! """
insertion_limit = int(1e5)
rotation_limit = 1000
package = 0.3
distance = {"min": 0.08**2, "opt": 0.12**2}

parser = argparse.ArgumentParser()
parser.add_argument("--WIDTH_X", type=float, required=True)
parser.add_argument("--WIDTH_Y", type=float, required=True)
parser.add_argument("--H", type=float, required=True)
# parser.add_argument("--l", type=float, required=True)
parser.add_argument("--phi", type=float, required=True)
parser.add_argument("--theta", type=float, required=True)
parser.add_argument("--delta", type=float, required=True)
parser.add_argument("--interface_type", type=str, required=True)
parser.add_argument("--extention", type=str, required=True)
parser.add_argument("--folder", type=str, required=True)

args = parser.parse_args()

# Генерация искуственных датасетов

regions = [Droplet, Doughnut, Worm, Roll, Perforation]
types = ["droplet", "doughnut", "worm", "roll", "perforation"]
phi = np.array([0.1, 0.2, 0.3, 0.5, 0.8])
WIDTH_X = np.array([15, 15, 15, 15, 15])
WIDTH_Y = np.array([15, 15, 15, 15, 15])
H = np.array([7.5, 7.5, 7.5, 5, 7.5])
c = H / 2
c[0] = 0
c[2] = 0
thetas = [180, 170, 160, 150, 140, 130, 120, 110, 100, 90.1]

iterations = 11

for th in thetas:
    for id in range(len(regions)):
        # for id in [1]:
        box = np.array([WIDTH_X[id], WIDTH_Y[id], H[id]])

        names = ["point"]
        density = [12]  # nm-3 Плотность, как у воды

        for i in range(iterations, 100):
            structure = Structure(
                title="Points",
                box=box,
                atoms=np.empty(0, dtype=object),
                atoms_xyz=np.zeros((0, 3)),
            )
            points = structure.atoms_xyz

            # XZ surface
            region = regions[id](
                center=np.array([WIDTH_X[id] / 2, WIDTH_Y[id] / 2, c[id]]),
                borders=np.array([WIDTH_X[id], WIDTH_Y[id], H[id]]),
                phi=phi[id],
                th=np.deg2rad(th),
                H=H[id],
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

            os.system(f"mkdir synthetic/{th}")
            os.system(f"mkdir synthetic/{th}/{types[id]}")

            with open(f"synthetic/{th}/{types[id]}/{i+1}.gro", "w") as f:
                f.write(write_gro(structure))
