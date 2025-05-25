import sys
import numpy as np
from scipy.spatial.transform import Rotation

from .grid import check_grid
from .pbc import distance_pbc, distance_pbc2


def insert_point_into_shape(
    shape,
    points,
    point_id,
    system_size,
    mol_size,
    mol_id,
    insertion_limit=int(1e5),
    package=0.4,
):
    if mol_size < 0.05:
        mol_size = 0.05

    overlap = True
    insertion_counter = 0

    if not mol_id:
        return shape.generate_point()

    while overlap and (insertion_counter < insertion_limit):
        insertion_counter += 1

        if 5 * insertion_counter % insertion_limit == 0:
            print(f"Please, wait ... mol #{mol_id} is tried to be inserted\n")

        new_point = shape.generate_point()
        overlap = np.sum(
            distance_pbc2(new_point, points[0:point_id, :], system_size)
            < ((2 - insertion_counter / insertion_limit) * package * mol_size)
        )

    if insertion_counter > insertion_limit:
        sys.exit("Decrease package")

    return new_point


def find_position(
    all_xyz,  # all coordinates in the system (N, 3)
    new_point,  # where to insert
    atom_id,  # index to insert at
    mol_xyz,  # coordinates of the molecule to insert (M, 3)
    mol_id,  # id of the molecule
    rotation_limit=10,
    min_dist2=0.08**2,
    system_size=None,
):
    overlap = True
    rotation_counter = 0

    if not mol_id:
        new_mol_xyz = rotate_random(mol_xyz) + new_point
        return new_mol_xyz

    while overlap and (rotation_counter < rotation_limit):
        rotation_counter += 1
        new_mol_xyz = rotate_random(mol_xyz) + new_point
        overlap = False
        for atom_xyz in new_mol_xyz:
            if np.any(
                distance_pbc2(atom_xyz, all_xyz[0:atom_id, :], system_size) < min_dist2
            ):
                overlap = True
                break
    return new_mol_xyz


def rotate_random(xyz):
    theta = np.random.random(3) * 2 * np.pi
    rot = Rotation.from_euler("xyz", theta)
    new_xyz = rot.apply(xyz.copy())
    return new_xyz
