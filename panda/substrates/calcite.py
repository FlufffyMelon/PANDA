import numpy as np
import string
import os
import mdtraj as md
from mdtraj.core.topology import Topology
from tqdm import tqdm
from itertools import combinations
from panda.utils import apply_pbc
from .utils import base62_encode, get_box_vectors, CustomSubstrate


class CalciteSubstrate(CustomSubstrate):
    def __init__(
        self, unitcell_path: str, Lx: float, Ly: float, Lz: float, build: bool = True
    ):
        self.unitcell_path = unitcell_path
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        self.build = build
        self.gro_path, self.itp_path, self.ndx_path = self._generate()

    def _generate(self):
        print("\nGenerating substrate...")
        substr_name = generate_calcite_substrate(
            self.unitcell_path,
            self.Lx,
            self.Ly,
            self.Lz,
            build=self.build,
        )
        substr_folder, _ = os.path.split(self.unitcell_path)
        substr_path = os.path.join(substr_folder, "gro", substr_name)

        substr_itp_name = generate_calcite_itp(substr_path, build=self.build)
        substr_itp_path = os.path.join(substr_folder, "itp", substr_itp_name)

        substr_ndx_name = generate_calcite_ndx(substr_path, build=self.build)
        substr_ndx_path = os.path.join(substr_folder, "ndx", substr_ndx_name)

        return substr_path, substr_itp_path, substr_ndx_path


# def generate_calcite_substrate(
#     unitcell_path: str,
#     Lx: float,
#     Ly: float,
#     Lz: float,
#     build: bool = True,
# ):
#     """
#     Generate a substrate from a given unitcell.

#     Parameters
#     ----------
#     unitcell_path : str
#         The path to the unitcell .gro file.
#     Lx : float
#         The desired length of the substrate in the x direction.
#     Ly : float
#         The desired length of the substrate in the y direction.
#     Lz : float
#         The desired length of the substrate in the z direction.
#     freeze_substr : bool
#         Whether or not to freeze the substrate. If True, the names of the atoms
#         will not be changed. Default is False.
#     build : bool
#         Whether to build the substrate if it does not exist. Default is True.

#     Returns
#     -------
#     str
#         The path to the file of the generated substrate.
#     """

#     unitcell = md.load(unitcell_path, top=unitcell_path)

#     unitcell_box = unitcell.unitcell_lengths[0]

#     Nx, Ny, Nz = np.round(
#         np.clip(np.array([Lx, Ly, Lz]) / unitcell_box[:3], 1, None)
#     ).astype(int)
#     N = Nx * Ny * Nz
#     N_atoms = unitcell.n_atoms
#     ex, ey, ez = get_box_vectors(unitcell_box)
#     print(
#         f"Creating substrate with dimensions {unitcell_box[0] * Nx:.1f}x{unitcell_box[1] * Ny:.1f}x{unitcell_box[2] * Nz:.1f}  ({Nx}x{Ny}x{Nz})..."
#     )

#     # Naming variables
#     folder, unitcell_filename = os.path.split(unitcell_path)
#     unitcell_name = os.path.splitext(unitcell_filename)[0]
#     filename = "_".join(unitcell_name.split("_")[:-1])
#     output_name = f"{filename}_{Nx}x{Ny}x{Nz}.gro"
#     output_path = os.path.join(folder, "gro", output_name)

#     if build:
#         if os.path.isfile(output_path):
#             print(f"A ready-made substrate is used from `{folder}/gro`")
#             return output_name
#         else:
#             print(
#                 "A substrate with this size does not exist. Forcibly generating it...."
#             )

#     substrate = generate_substrate(unitcell, Lx, Ly, Lz, 'CAL')
#     top = substrate.topology

#     # Rename atoms to be unique
#     counter_dict = {"C": 0, "O": 0, "Ca": 0}
#     for idx in range(N * N_atoms):
#         atom = top.atom(idx)

#         base = "".join([c for c in atom.name if not c.isdigit()])
#         counter_dict[base] += 1
#         atom.name = base + base62_encode(
#             counter_dict[base], length=(3 if base == "Ca" else 4)
#         )

#     # Build new topology
#     # top = Topology()
#     # residue = top.add_residue("CAL", top.add_chain())
#     # for i in range(N * N_atoms):
#     #     top.add_atom(atomnames[i], element=None, residue=residue, serial=i + 1)

#     # Set box
#     # box = unitcell_box[:3] * np.array([Nx, Ny, Nz])
#     # all_xyz = apply_pbc(all_xyz, box)
#     # traj = md.Trajectory(
#     #     all_xyz.reshape(1, -1, 3),
#     #     top,
#     #     unitcell_lengths=box.reshape(1, 3),
#     #     unitcell_angles=np.array([[90.0, 90.0, 90.0]]),
#     # )

#     # Write .gro file
#     os.makedirs(os.path.dirname(output_path), exist_ok=True)
#     substrate.save_gro(output_path)
#     print("Substrate successfully created!")

#     return output_name


def generate_calcite_substrate(
    unitcell_path: str,
    Lx: float,
    Ly: float,
    Lz: float,
    build: bool = True,
):
    """
    Generate a substrate from a given unitcell.

    Parameters
    ----------
    unitcell_path : str
        The path to the unitcell .gro file.
    Lx : float
        The desired length of the substrate in the x direction.
    Ly : float
        The desired length of the substrate in the y direction.
    Lz : float
        The desired length of the substrate in the z direction.
    freeze_substr : bool
        Whether or not to freeze the substrate. If True, the names of the atoms
        will not be changed. Default is False.
    build : bool
        Whether to build the substrate if it does not exist. Default is True.

    Returns
    -------
    str
        The path to the file of the generated substrate.
    """

    unitcell = md.load(unitcell_path, top=unitcell_path)
    unitcell_box = unitcell.unitcell_lengths[0]

    Nx, Ny, Nz = np.round(
        np.clip(np.array([Lx, Ly, Lz]) / unitcell_box[:3], 1, None)
    ).astype(int)
    N = Nx * Ny * Nz
    N_atoms = unitcell.n_atoms
    ex, ey, ez = get_box_vectors(unitcell_box)
    print(
        f"Creating substrate with dimensions {unitcell_box[0] * Nx:.1f}x{unitcell_box[1] * Ny:.1f}x{unitcell_box[2] * Nz:.1f}  ({Nx}x{Ny}x{Nz})..."
    )

    # Naming variables
    folder, unitcell_filename = os.path.split(unitcell_path)
    unitcell_name = os.path.splitext(unitcell_filename)[0]
    filename = "_".join(unitcell_name.split("_")[:-1])
    output_name = f"{filename}_{Nx}x{Ny}x{Nz}.gro"
    output_path = os.path.join(folder, "gro", output_name)

    if build:
        if os.path.isfile(output_path):
            print(f"A ready-made substrate is used from `{folder}/gro`")
            return output_name
        else:
            print(
                "A substrate with this size does not exist. Forcibly generating it...."
            )

    # Prepare arrays for new coordinates and atom info
    all_xyz = np.zeros((N * N_atoms, 3))
    atomnames = []
    mol_ids = []

    # Replicate unitcell
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                index = i * Ny * Nz + j * Nz + k
                offset = ex * i + ey * j + ez * k
                for a in range(N_atoms):
                    idx = index * N_atoms + a
                    all_xyz[idx, :] = unitcell.xyz[0, a, :] + offset
                    atom = unitcell.topology.atom(a)
                    atomnames.append(atom.name)
                    mol_ids.append(index + 1)

    # Rename atoms to be unique
    counter_dict = {"C": 0, "O": 0, "Ca": 0}
    for idx, name in enumerate(atomnames):
        base = "".join([c for c in name if not c.isdigit()])
        counter_dict[base] += 1
        atomnames[idx] = base + base62_encode(
            counter_dict[base], length=(3 if base == "Ca" else 4)
        )

    # Build new topology
    top = Topology()
    top.add_residue("", top.add_chain())
    residue = top.add_residue("CAL", top.add_chain())
    for i in range(N * N_atoms):
        top.add_atom(atomnames[i], element=None, residue=residue, serial=i + 1)

    # Set box
    box = unitcell_box[:3] * np.array([Nx, Ny, Nz])
    all_xyz = apply_pbc(all_xyz, box)
    traj = md.Trajectory(
        all_xyz.reshape(1, -1, 3),
        top,
        unitcell_lengths=box.reshape(1, 3),
        unitcell_angles=np.array([[90.0, 90.0, 90.0]]),
    )

    # Write .gro file
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    traj.save_gro(output_path)
    print("Substrate successfully created!")

    return output_name


def generate_calcite_ndx(substr_path: str, build: bool = True):
    """
    Generate .ndx file for a calcite substrate, based on the .gro file.
    Creates one group: CALCITE (all substrate atoms).

    Parameters
    ----------
    substr_path : str
        The path to the .gro file of the substrate.
    build : bool
        Whether to build the ndx if it does not exist. Default is True.

    Returns
    -------
    str
        The name of the generated .ndx file.
    """
    # Setup paths
    folder, substr_filename = os.path.split(substr_path)
    folder = os.path.split(folder)[0]
    substr_name = os.path.splitext(substr_filename)[0]
    output_name = substr_name + ".ndx"
    output_path = os.path.join(folder, "ndx", output_name)

    # Check if the NDX already exists
    if build and os.path.isfile(output_path):
        print(f"A ready-made ndx for substrate is used from `{folder}/ndx`")
        return output_name
    else:
        print(
            "The ndx does not exist for this substrate size. Forcibly generating it...."
        )

    # Load the substrate
    substr = md.load(substr_path)

    # All atoms in the substrate belong to CALCITE group
    calcite_atoms = list(range(1, substr.n_atoms + 1))  # 1-indexed for GROMACS

    # Write the index file
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        f.write("[ CALCITE ]\n")
        # Write atoms in groups of 15 per line
        for i in range(0, len(calcite_atoms), 15):
            f.write(" ".join(map(str, calcite_atoms[i : i + 15])) + "\n")

    print(f"NDX file with {len(calcite_atoms)} CALCITE atoms created.")

    return output_name


def generate_calcite_itp(substr_path: str, build: bool = True):
    """
    Generate .itp file for a substrate, based on the .gro file.

    Parameters
    ----------
    substr_path : str
        The path to the .gro file of the substrate.

    Returns
    -------
    str
        The path to the file of the generated .itp file.

    """
    substr = md.load(substr_path)
    # Check for orthorhombic box: first 3 components positive, rest (if present) zero
    box_lengths = substr.unitcell_lengths[0]
    assert np.all(box_lengths[:3] > 0), "Box lengths must be positive."
    if box_lengths.shape[0] > 3:
        assert np.allclose(box_lengths[3:], 0), (
            "Box should be orthorhombic (last 6 box components should be zero)."
        )

    folder, substr_filename = os.path.split(substr_path)
    folder = os.path.split(folder)[0]
    substr_name = os.path.splitext(substr_filename)[0]
    output_name = substr_name + ".itp"
    filename = "_".join(substr_name.split("_")[:-1])
    output_path = os.path.join(folder, "itp", output_name)

    if build:
        if os.path.isfile(output_path):
            print(f"A ready-made itp for substrate is used from `{folder}/itp`")
            return output_name
        else:
            print(
                "The itp does not exist for this substrate size. Forcibly generating it...."
            )

    neigh_dict = get_calcite_neighbors_list_numpy(substr)

    # Writing text of [ atoms ] section
    atoms_text = ""
    counter_dict = {"C": 0, "O": 0, "Ca": 0}
    metadata = {
        "C": ["CCA", 12.011, 0.999],
        "O": ["OCA", 15.999, -0.889],
        "Ca": ["CA", 40.078, 1.668],
    }

    atom_names = [a.name for a in substr.topology.atoms]
    for i, label in np.ndenumerate(atom_names):
        # name = "".join([i for i in label.name if not i.isdigit()])
        name = label[:2] if label[:2] == "Ca" else label[0]
        counter_dict[name] += 1
        type, mass, charge = metadata[name]

        atoms_text += "{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>10}{:>10}\n".format(
            i[0] + 1, type, 1, "CAL", label, 1, charge, mass
        )

    # Writing text of [ constraints ] section
    constraints_text = ""
    # Bond C-O
    for key, values in neigh_dict.items():
        for v in values:
            # constraints_text += "{:>8}{:>8}{:>8}{:>8.3f}\n".format(
            #     key + 1, v + 1, 1, np.linalg.norm(rij(key, v, substr))
            # )
            constraints_text += "{:>8}{:>8}{:>8}\n".format(key + 1, v + 1, 1)

    # Writing text of [ angles ] section
    angles_text = ""
    # Angles O-C-O
    for key, value in neigh_dict.items():
        for i, j in combinations(value, 2):
            # angles_text += "{:>8}{:>8}{:>8}{:>8}{:>8.1f}{:>8}\n".format(
            #     i + 1, key + 1, j + 1, 1, angle(i, key, j, substr), 1852.0
            # )
            angles_text += "{:>8}{:>8}{:>8}{:>8}\n".format(i + 1, key + 1, j + 1, 1)

    # Writing text of [ dihedrals ] section
    dihedrals_text = ""
    # dihedral O1-C-O2 O1-C-O3
    for key, value in neigh_dict.items():
        dihedrals_text += "{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}\n".format(
            key + 1,
            value[0] + 1,
            value[1] + 1,
            value[2] + 1,
            5,
            0.0,
            2 * 28.9,
            0.0,
            0.0,
        )

    # Combining everything into itp template
    with open(os.path.join(folder, filename + "_template" + ".itp")) as f:
        final_text = f.read().format(
            **{
                "atoms": atoms_text,
                "constraints": constraints_text,
                "angles": angles_text,
                "dihedrals": dihedrals_text,
            }
        )

    # Writing final itp file
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        f.write(final_text)

    return output_name


def get_calcite_neighbors_list_numpy(substr: md.Trajectory):
    """
    Generate a dictionary of calcite atom neighbors.

    For each carbon atom in the structure, find all neighboring oxygen atoms
    within a specified cutoff distance and store them in a dictionary.

    Parameters
    ----------
    substr : mdtraj.Trajectory
        A Trajectory object containing atoms and their positions.

    Returns
    -------
    dict
        A dictionary where keys are indices of carbon atoms and values are lists
        of indices of neighboring oxygen atoms.
    """
    assert len(substr.unitcell_lengths[0]) == 3, "Box should be orthogonal"

    # Cutoff distance for determining neighboring atoms
    l = 0.118 + 0.022
    neigh_dict = dict()

    xyz = substr.xyz[0]
    top = substr.topology
    atom_names = [a.name for a in top.atoms]

    # Initialize masks and lists for oxygen and carbon atom indices
    oxygen_mask = np.zeros(substr.n_atoms, dtype=bool)
    oxygen_real_ids = []
    carbon_ids = []

    # Identify indices of oxygen and carbon atoms
    for i in range(substr.n_atoms):
        name = atom_names[i][:2] if atom_names[i][:2] == "Ca" else atom_names[i][0]

        if name == "O":
            oxygen_mask[i] = True
            oxygen_real_ids.append(i)
        elif name == "C":
            carbon_ids.append(i)
    oxygen_real_ids = np.array(oxygen_real_ids)

    print("Generating neighbors list")
    # Iterate over carbon atoms to find their neighboring oxygen atoms
    for i in tqdm(carbon_ids):
        # Calculate relative positions of potential neighboring oxygen atoms
        rij_vecs = xyz[oxygen_mask] - xyz[i, :]
        mask = np.abs(rij_vecs) >= substr.unitcell_lengths[0] / 2
        rij_vecs -= substr.unitcell_lengths[0] * mask * np.sign(rij_vecs)

        # Determine indices of oxygen atoms within the cutoff distance
        neigh_oxygen = np.argwhere(np.linalg.norm(rij_vecs, axis=1) < l).ravel()
        O_neigh = list(oxygen_real_ids[neigh_oxygen])

        # Ensure exactly 3 neighboring oxygen atoms are found
        assert len(O_neigh) == 3, f"Incorrect number of neighbours ({len(O_neigh)})!!!"

        # Store the neighbors in the dictionary
        neigh_dict[i] = O_neigh.copy()
    return neigh_dict
