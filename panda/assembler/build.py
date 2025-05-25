import sys
import numpy as np
from tqdm import tqdm
import mdtraj as md
from panda.assembler.insert import find_position, insert_point_into_shape
from .grid import add_mol_grid
from panda.utils import get_center_pbc


def build_system(
    dir,
    traj,  # mdtraj.Trajectory object
    names=None,
    numbers=None,
    shapes=None,
    substr_points=None,
    insertion_limit=int(1e5),
    rotation_limit=10,
    package=0.4,
    min_dist2=0.08**2,
):
    if names is None or numbers is None or shapes is None:
        sys.exit("Please, give all the parameters")

    if not (len(names) == len(numbers) == len(shapes)):
        sys.exit("Parameters have different lengths")

    system_size = traj.unitcell_lengths[0]
    print("Number of molecules:")
    for i, name in enumerate(names):
        print(name, "\t", numbers[i])

    atoms_number = []
    mol_trajs = []
    for i, name in enumerate(names):
        mol_traj = md.load(f"{dir}/{name}.gro")
        atoms_number.append(mol_traj.n_atoms)
        mol_trajs.append(mol_traj)
    atoms_number = np.array(atoms_number)

    # Prepare arrays for new atoms and coordinates
    total_new_atoms = np.sum(atoms_number * numbers)
    all_xyz = np.vstack((traj.xyz[0], np.zeros((total_new_atoms, 3))))
    points = np.vstack((substr_points, np.zeros((sum(numbers), 3))))
    point_id = substr_points.shape[0]
    atom_id = substr_points.shape[0]

    # Build new topology: flat, just residues and atoms
    new_top = md.Topology()
    # Copy substrate residues and atoms
    for res in traj.topology.residues:
        new_res = new_top.add_residue(res.name, new_top.add_chain())
        for atom in res.atoms:
            new_top.add_atom(atom.name, atom.element, new_res, serial=atom.serial)

    print("\nFilling system:")
    for i, name in enumerate(names):
        mol_traj = mol_trajs[i]
        mol_xyz = mol_traj.xyz[0]
        center = get_center_pbc(mol_xyz, system_size)
        mol_xyz -= center[0]
        mol_size = np.max(np.linalg.norm(mol_xyz, axis=1))
        print(name, mol_size)

        for mol_id in tqdm(range(numbers[i]), desc=name):
            new_point = insert_point_into_shape(
                shapes[i],
                points,
                point_id,
                system_size,
                mol_size,
                mol_id,
                insertion_limit,
                package,
            )
            points[point_id, :] = new_point
            point_id += 1

            # Find position and orientation for the molecule
            new_mol_xyz = find_position(
                all_xyz,
                new_point,
                atom_id,
                mol_xyz,
                mol_id,
                rotation_limit,
                min_dist2,
                system_size,
            )
            all_xyz[atom_id : atom_id + mol_xyz.shape[0], :] = new_mol_xyz

            # Add molecule residues and atoms
            for res in mol_traj.topology.residues:
                new_res = new_top.add_residue(res.name, new_top.add_chain())
                for atom in res.atoms:
                    new_top.add_atom(atom.name, atom.element, new_res)
            atom_id += mol_xyz.shape[0]

    # Finalize coordinates
    all_xyz = all_xyz[: len(list(new_top.atoms))]
    new_traj = md.Trajectory(
        xyz=all_xyz[np.newaxis, :, :],
        topology=new_top,
        unitcell_lengths=traj.unitcell_lengths,
        unitcell_angles=traj.unitcell_angles,
    )
    return new_traj
