import numpy as np
from scipy.spatial.transform import Rotation


def insert_point(
    shape,
    points,
    point_grid,
    mol_size,
    insertion_limit=int(1e5),
    package=0.4,
):
    if mol_size < 0.05:
        mol_size = 0.05

    insertion_counter = 0

    while insertion_counter < insertion_limit:
        insertion_counter += 1

        if 5 * insertion_counter % insertion_limit == 0:
            print(f"Please, wait ... mol is tried to be inserted\n")

        new_point = shape.generate_point()
        min_dist2 = ((2 - insertion_counter / insertion_limit) * package * mol_size)
        if not point_grid.check_collision(new_point, points, min_dist2):
            return new_point

    return None


def insert_atom(
    all_xyz,
    new_point,
    mol_xyz,
    atom_grid,
    rotation_limit=10,
    min_dist2=0.08**2,
):
    rotation_counter = 0

    while rotation_counter < rotation_limit:
        rotation_counter += 1
        new_mol_xyz = rotate_random(mol_xyz) + new_point
        collision = False
        for atom_xyz in new_mol_xyz:
            if np.any(atom_grid.check_collision(atom_xyz, all_xyz, min_dist2)):
                collision = True
                break

        if not collision:
            return new_mol_xyz

    return None


def rotate_random(xyz):
    theta = np.random.random(3) * 2 * np.pi
    rot = Rotation.from_euler("xyz", theta)
    new_xyz = rot.apply(xyz.copy())
    return new_xyz
