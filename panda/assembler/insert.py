import numpy as np
from scipy.spatial.transform import Rotation


def insert_point(
    shape,
    points,
    point_grid,
    mol_size,
    insertion_limit=int(1e5),
    min_dist2=None,
    package=0.4,
):
    if mol_size < 0.05:
        mol_size = 0.05

    single_atom = min_dist2 is None

    insertion_counter = 0

    while insertion_counter < insertion_limit:
        insertion_counter += 1

        if 5 * insertion_counter % insertion_limit == 0:
            print("    [Insert] Please, wait ... molecule is being inserted...")

        new_point = shape.generate_point()
        if single_atom:
            min_dist2 = (2 - insertion_counter / insertion_limit) * package * mol_size
        if not point_grid.check_collision(new_point, points, min_dist2):
            return new_point

    return None


def insert_atom(
    all_xyz,
    new_point,
    mol_xyz,
    atom_grid,
    rotation_limit=10,
    # single_atom=True,
    min_dist2=0.08**2,
):
    rotation_counter = 0

    while rotation_counter < rotation_limit:
        rotation_counter += 1

        # if single_atom:
        #     new_mol_xyz = mol_xyz + new_point
        # else:
        #     new_mol_xyz = rotate_random(mol_xyz) + new_point
        new_mol_xyz = rotate_random(mol_xyz) + new_point

        collision = False
        for atom_xyz in new_mol_xyz:
            if np.any(atom_grid.check_collision(atom_xyz, all_xyz, min_dist2)):
                collision = True
                break

        if not collision:
            return new_mol_xyz

        # if single_atom:
        #     break

    return None


def rotate_random(xyz):
    theta = np.random.random(3) * 2 * np.pi
    rot = Rotation.from_euler("xyz", theta)
    new_xyz = rot.apply(xyz.copy())
    return new_xyz
