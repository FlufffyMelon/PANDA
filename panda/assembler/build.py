import numpy as np
from tqdm import tqdm
import mdtraj as md
from .insert import insert_point, insert_atom
from .grid import Grid
from panda.utils import get_center_pbc


def build(
    mol_path,
    initial_traj,  # mdtraj.Trajectory object
    mol_names=None,
    mol_numbers=None,
    shapes=None,
    insertion_limit=int(1e5),
    rotation_limit=10,
    package=0.4,
    insertion_attempts=10,
    min_dist2=0.08**2,
):
    if mol_names is None or mol_numbers is None or shapes is None:
        print("\n[Error] Missing parameters: mol_names, mol_numbers, or shapes.\n")
        raise ValueError("Please, give all the parameters")

    if not (len(mol_names) == len(mol_numbers) == len(shapes)):
        print(
            "\n[Error] Parameters have different lengths: names=%d, numbers=%d, shapes=%d\n"
            % (len(mol_names), len(mol_numbers), len(shapes))
        )
        raise ValueError("Parameters have different lengths")

    system_size = initial_traj.unitcell_lengths[0]
    print(
        f"[Build] Number of molecules: "
        + ", ".join(f"{name}={mol_numbers[i]}" for i, name in enumerate(mol_names))
    )

    atoms_number = []
    mol_trajs = []
    mol_sizes = []
    for i, name in enumerate(mol_names):
        mol_traj = md.load(f"{mol_path}/{name}.gro")
        atoms_number.append(mol_traj.n_atoms)
        mol_trajs.append(mol_traj)

        mol_xyz = mol_traj.xyz
        mol_unitcell = mol_traj.unitcell_lengths
        center = get_center_pbc(mol_xyz, mol_unitcell)
        mol_xyz -= center
        mol_size = np.max(np.linalg.norm(mol_xyz, axis=2))
        mol_sizes.append(mol_size)
    atoms_number = np.array(atoms_number)
    max_diameter = 2 * max(mol_sizes)

    # Prepare arrays for new atoms and coordinates
    total_new_atoms = np.sum(atoms_number * mol_numbers)
    all_xyz = np.vstack((initial_traj.xyz[0], np.zeros((total_new_atoms, 3))))
    points = np.vstack((initial_traj.xyz[0], np.zeros((sum(mol_numbers), 3))))
    point_id = initial_traj.xyz[0].shape[0]
    atom_id = initial_traj.xyz[0].shape[0]

    # Build new topology: flat, just residues and atoms
    new_top = md.Topology()
    # Copy substrate residues and atoms
    for res in initial_traj.topology.residues:
        new_res = new_top.add_residue(res.name, new_top.add_chain())
        for atom in res.atoms:
            new_top.add_atom(atom.name, atom.element, new_res, serial=atom.serial)

    # Create grid for points and atoms
    point_grid = Grid(system_size, max_diameter, label="points")
    atom_grid = Grid(system_size, 3 * np.sqrt(min_dist2) / 2, label="atoms")
    for i in range(point_id):
        point_grid.add(points[i], i)
    for i in range(atom_id):
        atom_grid.add(all_xyz[i], i)
    print("[Build] Filling system...")

    for i, name in enumerate(mol_names):
        mol_traj = mol_trajs[i]
        mol_xyz = mol_traj.xyz[0]
        mol_size = mol_sizes[i]

        for _ in tqdm(range(mol_numbers[i]), desc=name):
            attempt = 0
            while attempt <= insertion_attempts:
                new_point = insert_point(
                    shapes[i], points, point_grid, mol_size, insertion_limit, package
                )
                if new_point is None:
                    attempt += 1
                    continue

                # Find position and orientation for the molecule
                new_mol_xyz = insert_atom(
                    all_xyz, new_point, mol_xyz, atom_grid, rotation_limit, min_dist2
                )
                if new_mol_xyz is None:
                    attempt += 1
                    continue

                break

            if new_point is None:
                print(
                    f"[Error] Could not insert molecule: decrease package ({package:.2f}) or increase box size\n"
                )
                raise RuntimeError(f"Decrease package ({package}) or increase box size")
            if new_mol_xyz is None:
                print(
                    f"[Error] Could not place molecule after rotation_limit ({rotation_limit}) attempts\n"
                )
                raise RuntimeError(
                    f"Could not place molecule after rotation_limit ({rotation_limit}) attempts"
                )

            # Adding new point
            points[point_id, :] = new_point
            point_grid.add(new_point, point_id)
            point_id += 1

            # Adding atoms of new molecular
            all_xyz[atom_id : atom_id + mol_xyz.shape[0], :] = new_mol_xyz
            for j in range(mol_xyz.shape[0]):
                atom_grid.add(new_mol_xyz[j], atom_id + j)

            # Add molecule residues and atoms
            for res in mol_traj.topology.residues:
                new_res = new_top.add_residue(res.name, new_top.add_chain())
                for atom in res.atoms:
                    new_top.add_atom(atom.name, atom.element, new_res)
            atom_id += mol_xyz.shape[0]

    # Finalize coordinates
    # all_xyz = all_xyz[: len(list(new_top.atoms))]
    new_traj = md.Trajectory(
        xyz=all_xyz[np.newaxis, :, :],
        topology=new_top,
        unitcell_lengths=initial_traj.unitcell_lengths,
        unitcell_angles=initial_traj.unitcell_angles,
    )
    return new_traj
