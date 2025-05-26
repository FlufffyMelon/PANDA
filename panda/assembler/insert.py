import sys
import numpy as np
from scipy.spatial.transform import Rotation
from .utils import distance_pbc, distance_pbc2
from panda.utils import validate_positions_and_box


def insert_point(
    shape,
    points,
    point_id,
    point_grid,
    system_size,
    mol_size,
    mol_id,
    insertion_limit=int(1e5),
    package=0.4,
):
    if mol_size < 0.05:
        mol_size = 0.05

    insertion_counter = 0

    while insertion_counter < insertion_limit:
        insertion_counter += 1

        if 5 * insertion_counter % insertion_limit == 0:
            print(f"Please, wait ... mol #{mol_id} is tried to be inserted\n")

        new_point = shape.generate_point()
        min_dist2 = ((2 - insertion_counter / insertion_limit) * package * mol_size)
        # overlap = np.sum(
        #     distance_pbc2(new_point, points[0:point_id, :], system_size)
        #     < ((2 - insertion_counter / insertion_limit) * package * mol_size)
        # )
        if not point_grid.check_collision(new_point, points, min_dist2):
            return new_point
        # overlap =  point_grid.check_collision(new_point, points, min_dist2)

    return None
    # sys.exit(f"Decrease package ({package}) or increase box size (insert_point)")


def insert_atom(
    all_xyz,  # all coordinates in the system (N, 3)
    new_point,  # where to insert
    atom_id,  # index to insert at
    mol_xyz,  # coordinates of the molecule to insert (M, 3)
    mol_id,  # id of the molecule
    atom_grid,
    rotation_limit=10,
    min_dist2=0.08**2,
    system_size=None,
):
    rotation_counter = 0

    while rotation_counter < rotation_limit:
        rotation_counter += 1
        new_mol_xyz = rotate_random(mol_xyz) + new_point
        collision = False
        for atom_xyz in new_mol_xyz:
            # if np.any(
            #     distance_pbc2(atom_xyz, all_xyz[0:atom_id, :], system_size) < min_dist2
            # ):
            if np.any(atom_grid.check_collision(atom_xyz, all_xyz, min_dist2)):
                collision = True
                break

        if not collision:
            return new_mol_xyz

    return None
    # sys.exit(
    #     f"Could not place molecule after rotation_limit ({rotation_limit}) attempts (insert_atom)"
    # )


def rotate_random(xyz):
    theta = np.random.random(3) * 2 * np.pi
    rot = Rotation.from_euler("xyz", theta)
    new_xyz = rot.apply(xyz.copy())
    return new_xyz
