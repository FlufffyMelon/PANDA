import numpy as np
from tqdm import tqdm
from .utils import delta_pbc
from panda.utils import apply_pbc
from .grid import Grid


def mixer(
    traj, grid_size, min_dist2, opt_dist2, iteration_lim=10, check_same_molecule=False
):
    """
    Mixes the traj by pushing atoms apart if they are too close, using a Grid for fast neighbor search.
    Args:
        traj: mdtraj.Trajectory object
        grid_size: cell size for the grid
        min_dist2: minimum allowed squared distance between atoms of different molecules
        opt_dist2: maximum allowed squared distance for random displacement
        iteration_lim: maximum number of iterations
        check_same_molecule: If True, check and resolve collisions between atoms of the same molecule.
                             If False, only check and resolve collisions between atoms of different molecules (default).
    Returns:
        mdtraj.Trajectory: mixed traj
    """
    xyz = traj.xyz[0].copy()
    box = traj.unitcell_lengths[0]
    n_atoms = xyz.shape[0]

    # Build mol_ids: assign a unique id to each residue's atoms
    mol_ids = np.zeros(n_atoms, dtype=int)
    mol_counter = 0
    for res in traj.topology.residues:
        for atom in res.atoms:
            mol_ids[atom.index] = mol_counter
        mol_counter += 1

    min_dist = np.sqrt(min_dist2)
    opt_dist = np.sqrt(opt_dist2)

    # Build grid for initial positions
    grid = Grid(box, grid_size, label="mixing")
    for idx in range(n_atoms):
        grid.add(xyz[idx], idx)

    ins = 0
    overlap = True
    while overlap and (ins < iteration_lim):
        overlap = False
        overlap_counter = 0
        ins += 1
        print(f"[Mixing] Iteration {ins}")

        for i in tqdm(range(n_atoms), desc="Atoms"):
            # First check collision using the grid function
            if grid.check_collision(xyz[i], xyz, min_dist2):
                # If collision detected, loop over neighbor cells to find specific colliding atoms
                idx_i = grid.get_cell_index(xyz[i])
                for nidx in grid.neighbor_cell_indices[idx_i]:
                    # Get atom indices in the current neighbor cell, excluding atom i
                    if idx_i == nidx:
                        neighbor_atom_ids_in_cell = [
                            j for j in grid.grid[nidx] if j != i
                        ]
                    else:
                        neighbor_atom_ids_in_cell = grid.grid[nidx]
                    # neighbor_atom_ids_in_cell = [j for j in grid.grid[nidx] if j < i]

                    if not neighbor_atom_ids_in_cell:
                        continue

                    neighbor_xyz_in_cell = xyz[neighbor_atom_ids_in_cell]

                    # Calculate deltas and squared distances for atoms in this neighbor cell
                    deltas = delta_pbc(xyz[i], neighbor_xyz_in_cell, box)
                    d2s = np.sum(deltas**2, axis=1)

                    # Find indices of atoms in this neighbor cell that are too close
                    close_indices_in_cell = np.argwhere(d2s < min_dist2).flatten()

                    # Process colliding atoms found in this neighbor cell
                    for cell_idx in close_indices_in_cell:
                        j = neighbor_atom_ids_in_cell[cell_idx]

                        # Option to skip collisions within the same molecule
                        if not check_same_molecule and mol_ids[i] == mol_ids[j]:
                            continue  # skip same molecule

                        overlap = True
                        overlap_counter += 1
                        old_pos = xyz[j].copy()

                        # Get the specific delta for this colliding atom
                        delta = deltas[cell_idx]
                        dist = np.linalg.norm(delta)

                        # Push atom j apart
                        disp = delta / dist * np.random.uniform(min_dist, opt_dist)
                        xyz[j] += disp
                        # xyz[j] = xyz[j] % box
                        grid.move(old_pos, xyz[j], j)

        print(f"[Mixing] Detected {overlap_counter} overlaps")

    traj.xyz = apply_pbc(xyz, box)
    return traj
