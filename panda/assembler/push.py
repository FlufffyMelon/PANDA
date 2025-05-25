import numpy as np
from tqdm import tqdm
from .pbc import distance_pbc2, delta_pbc


def push_atoms_apart(xyz, box, mol_ids, min_dist2, max_dist2, iteration_lim=10):
    min_dist = np.sqrt(min_dist2)
    max_dist = np.sqrt(max_dist2)

    new_xyz = xyz.copy()
    ins = 0
    overlap = True
    while overlap and (ins < iteration_lim):
        overlap = False
        overlap_counter = 0
        ins += 1

        print(f"Iteration {ins}")
        for i in tqdm(range(len(new_xyz) - 1), desc="Atoms"):
            for j in range(i + 1, len(new_xyz)):
                if mol_ids[i] != mol_ids[j]:
                    delta = delta_pbc(new_xyz[i], [new_xyz[j]], box)[0]
                    if np.all(abs(delta) < min_dist):
                        overlap = True
                        overlap_counter += 1
                        new_xyz[j] += (
                            delta
                            / np.linalg.norm(delta)
                            * np.random.uniform(min_dist, max_dist)
                        )

        print(f"{overlap_counter} overlaps detected")

    # Apply periodic boundary conditions
    new_xyz = new_xyz % box
    return new_xyz
