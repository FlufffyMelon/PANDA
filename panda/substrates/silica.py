import numpy as np
import os
import mdtraj as md
from mdtraj.core.topology import Topology
from panda.utils import apply_pbc
from .utils import base62_encode, CustomSubstrate, ase2mdtraj
from ase.io import read
from ase.build import surface


class SilicaSubstrate(CustomSubstrate):
    def __init__(
        self,
        unitcell_path: str,
        Lx: float,
        Ly: float,
        Lz: float,
        silanol_density: float = 0.0,
        build: bool = True,
    ):
        self.unitcell_path = unitcell_path
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        self.build = build
        self.silanol_density = silanol_density
        self.gro_path, self.itp_path, self.ndx_path = self._generate()

    def _generate(self):
        print("\nGenerating substrate...")
        substr_name = generate_silica_substrate(
            self.unitcell_path,
            self.Lx,
            self.Ly,
            self.Lz,
            self.silanol_density,
            build=self.build,
        )
        substr_folder, _ = os.path.split(self.unitcell_path)
        substr_path = os.path.join(substr_folder, "gro", substr_name)

        substr_itp_name = generate_silica_itp(substr_path, build=self.build)
        substr_itp_path = os.path.join(substr_folder, "itp", substr_itp_name)

        substr_ndx_name = generate_silica_ndx(substr_path, build=self.build)
        substr_ndx_path = os.path.join(substr_folder, "ndx", substr_ndx_name)

        return substr_path, substr_itp_path, substr_ndx_path


def generate_silica_substrate(
    unitcell_path: str,
    Lx: float,
    Ly: float,
    Lz: float,
    silanol_density: float,
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

    atoms = read(unitcell_path)
    unitcell_box = np.diag(atoms.cell) * 1e-1

    Nx, Ny, Nz = np.round(
        np.clip(np.array([Lx, Ly, Lz]) / unitcell_box, 1, None)
    ).astype(int)
    supercell = atoms.repeat((Nx, Ny, 1))
    slab = surface(supercell, (0, 0, 1), layers=Nz, vacuum=0.0)
    supercell = ase2mdtraj(slab)
    N_atoms = supercell.n_atoms
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

    # Rename atoms to be unique
    counter_dict = {"Si": 0, "O": 0}
    residue = supercell.topology.add_residue(
        "SIL", supercell.topology.add_chain(), resSeq=1
    )
    for idx in range(N_atoms):
        atom = supercell.topology.atom(idx)

        base = "".join([c for c in atom.name if not c.isdigit()])
        counter_dict[base] += 1
        atom.name = base + base62_encode(
            counter_dict[base], length=(3 if base == "Si" else 4)
        )
        atom.residue = residue

    # Set box
    box = np.diag(supercell.unitcell_vectors[0, :, :])
    shift = np.array([0.0, 0.0, -0.245])  # shift by -0.245 nm along z-axis
    all_xyz = supercell.xyz[0, :, :] + shift
    all_xyz = apply_pbc(all_xyz, box)

    # Removeing top layer of Si atoms
    offset = 0.1
    new_top = Topology()
    new_top.add_residue("", new_top.add_chain())
    residue = new_top.add_residue("SIL", new_top.add_chain())
    unique_id = 0
    for idx in range(N_atoms):
        if supercell.topology.atom(idx).name[:2] == "Si":
            if all_xyz[0, idx, 2] > box[2] - offset:
                continue

        new_top.add_atom(
            supercell.topology.atom(idx).name,
            element=None,
            residue=residue,
            serial=unique_id + 1,
        )
        unique_id += 1
    print(f"Removed {N_atoms - unique_id} Si atoms")
    new_xyz = all_xyz[:, all_xyz[0, :, 2] < box[2] - offset, :]

    traj = md.Trajectory(
        new_xyz.reshape(1, -1, 3),
        new_top,
        unitcell_lengths=box.reshape(1, 3),
        unitcell_angles=np.array([[90.0, 90.0, 90.0]]),
    )

    # Adding silanol atoms
    traj_silanol = add_silanol_atoms(traj, silanol_density)

    # Write .gro file
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    traj_silanol.save_gro(output_path)
    print("Substrate successfully created!")

    return output_name


def add_silanol_atoms(traj, silanol_density):
    box = traj.unitcell_lengths[0]
    area = box[0] * box[1]
    N_silanol = int(area * silanol_density)
    print(f"Adding {N_silanol} silanol atoms")
    top_ids = []
    bottom_ids = []
    for idx in range(traj.n_atoms):
        if traj.topology.atom(idx).name[0] == "O":
            if traj.xyz[0, idx, 2] > box[2] - 0.15:
                top_ids.append(idx)
            elif traj.xyz[0, idx, 2] < 0.05:
                bottom_ids.append(idx)

    print(f"Found {len(top_ids)} top O atoms and {len(bottom_ids)} bottom O atoms")

    assert N_silanol <= len(top_ids) and N_silanol <= len(bottom_ids), (
        "Not enough O atoms to add silanol atoms (is not implemented yet)"
    )

    top_O_ids = np.random.choice(top_ids, size=N_silanol, replace=False)
    bottom_O_ids = np.random.choice(bottom_ids, size=N_silanol, replace=False)

    silanol_length = np.array([0, 0, 0.095])
    xyz_silanol = np.zeros((2 * N_silanol, 3))
    new_top = traj.topology

    for i, idx in enumerate(top_O_ids):
        xyz_silanol[i, :] = traj.xyz[0, idx, :] + silanol_length
        name = "H" + base62_encode(i, length=4)
        new_top.add_atom(
            name,
            element=None,
            residue=new_top.atom(idx).residue,
            serial=traj.n_atoms + i + 1,
        )

    for i, idx in enumerate(bottom_O_ids):
        xyz_silanol[N_silanol + i, :] = traj.xyz[0, idx, :] - silanol_length
        name = "H" + base62_encode(N_silanol + i, length=4)
        new_top.add_atom(
            name,
            element=None,
            residue=new_top.atom(idx).residue,
            serial=traj.n_atoms + N_silanol + i + 1,
        )

    new_xyz = np.concatenate([traj.xyz[0, :, :], xyz_silanol], axis=0)
    new_xyz = new_xyz + np.array([0, 0, 0.095])  # shift by 0.095 nm along z-axis
    box[2] += 0.065  # increase the box size by 0.065 nm along z-axis
    traj_silanol = md.Trajectory(
        new_xyz.reshape(1, -1, 3),
        new_top,
        unitcell_lengths=box.reshape(1, 3),
        unitcell_angles=np.array([[90.0, 90.0, 90.0]]),
    )

    return traj_silanol


def classify_silica_atoms(substr: md.Trajectory):
    """
    Classify silica atoms and identify bonds and angles.

    This function analyzes a silica substrate to:
    1. Classify oxygen atoms as:
       - Ob: bridging oxygen (Si-O-Si)
       - Oh: hydroxyl oxygen (Si-O-H)
       - O: surface oxygen (Si-O only)
    2. Identify Si-O and O-H bonds
    3. Identify Si-O-Si angles

    Parameters
    ----------
    substr : mdtraj.Trajectory
        A Trajectory object containing silica atoms and their positions.

    Returns
    -------
    tuple
        (atom_types, bonds, angles)
        - atom_types: list of atom type strings
        - bonds: list of (i,j) tuples for bonds
        - angles: list of (i,j,k) tuples for angles
    """
    xyz = substr.xyz[0]
    box = substr.unitcell_lengths[0]
    n_atoms = substr.n_atoms
    top = substr.topology

    # --- Atom categorization ---
    # Identify atoms by element type using vectorized operations
    atom_names = [a.name for a in top.atoms]
    atom_names_array = np.array(atom_names)

    is_si = np.array([name.startswith("Si") for name in atom_names_array])
    is_h = np.array([name.startswith("H") for name in atom_names_array])
    is_o = ~(is_si | is_h)

    si_atoms = np.where(is_si)[0].tolist()
    o_atoms = np.where(is_o)[0].tolist()
    h_atoms = np.where(is_h)[0].tolist()

    # Initialize atom types
    atom_types = [""] * n_atoms
    for i in si_atoms:
        atom_types[i] = "Si"
    for i in o_atoms:
        atom_types[i] = "O"  # Will be refined later
    for i in h_atoms:
        atom_types[i] = "H"

    print(
        f"Found {len(si_atoms)} Si atoms, {len(o_atoms)} O atoms, {len(h_atoms)} H atoms"
    )

    # --- Bond detection ---
    print("Finding Si-O bonds...")
    bonds = []
    o_neighbors = {i: [] for i in o_atoms}
    si_o_cutoff = 0.2  # nm
    o_h_cutoff = 0.1  # nm

    # Si-O bonds with periodic boundary conditions in x,y
    if si_atoms and o_atoms:
        # Compute pairwise distances efficiently
        o_coords = xyz[o_atoms]
        si_coords = xyz[si_atoms]

        # Calculate displacement vectors with PBC
        rij_vectors = o_coords[:, np.newaxis, :] - si_coords[np.newaxis, :, :]
        for dim in [0, 1]:  # Apply PBC in x,y only
            mask = np.abs(rij_vectors[:, :, dim]) >= box[dim] / 2
            rij_vectors[:, :, dim][mask] -= box[dim] * np.sign(
                rij_vectors[:, :, dim][mask]
            )

        # Find pairs within cutoff
        distances = np.linalg.norm(rij_vectors, axis=2)
        within_cutoff = distances < si_o_cutoff
        o_indices, si_indices = np.where(within_cutoff)

        # Record bonds
        for idx in range(len(o_indices)):
            o_idx = o_atoms[o_indices[idx]]
            si_idx = si_atoms[si_indices[idx]]
            # bonds.append((si_idx, o_idx))
            o_neighbors[o_idx].append(si_idx)

        print(f"Found {len(o_indices)} Si-O bonds")

    # O-H bonds (no PBC needed as H atoms are on surface)
    print("Finding O-H bonds...")
    if h_atoms and o_atoms:
        o_coords = xyz[o_atoms]
        h_coords = xyz[h_atoms]

        # Calculate distances between all O-H pairs
        oh_vectors = o_coords[:, np.newaxis, :] - h_coords[np.newaxis, :, :]
        oh_distances = np.linalg.norm(oh_vectors, axis=2)
        oh_within_cutoff = oh_distances < o_h_cutoff
        o_indices, h_indices = np.where(oh_within_cutoff)

        # Record bonds
        for idx in range(len(o_indices)):
            o_idx = o_atoms[o_indices[idx]]
            h_idx = h_atoms[h_indices[idx]]
            bonds.append((o_idx, h_idx))
            o_neighbors[o_idx].append(h_idx)

        print(f"Found {len(o_indices)} O-H bonds")

    # --- Oxygen classification ---
    print("Classifying oxygen atoms...")
    si_set = set(si_atoms)
    h_set = set(h_atoms)

    # Classify each oxygen atom based on its neighbors
    for o_idx in o_atoms:
        neighbors = o_neighbors[o_idx]
        si_count = sum(1 for n in neighbors if n in si_set)
        h_count = sum(1 for n in neighbors if n in h_set)

        if si_count == 2 and h_count == 0:
            atom_types[o_idx] = "Ob"  # Bridging oxygen (Si-O-Si)
        elif si_count == 1 and h_count == 1:
            atom_types[o_idx] = "Oh"  # Hydroxyl oxygen (Si-O-H)
        elif si_count == 1 and h_count == 0:
            atom_types[o_idx] = "O"  # Surface oxygen (Si-O only)
        else:
            print(
                f"Warning: Unusual oxygen bonding pattern - O atom {o_idx} has {si_count} Si neighbors and {h_count} H neighbors"
            )

    # --- Angle detection ---
    print("Finding Si-O-Si angles...")
    angles = []

    # Find angles only for bridging oxygens
    bridging_o = [o_idx for o_idx in o_atoms if atom_types[o_idx] == "Ob"]
    for o_idx in bridging_o:
        si_neighbors = [n for n in o_neighbors[o_idx] if n in si_atoms]
        if len(si_neighbors) == 2:
            angles.append((si_neighbors[0], o_idx, si_neighbors[1]))

    # --- Summary ---
    atom_types_array = np.array(atom_types)
    ob_count = np.sum(atom_types_array == "Ob")
    oh_count = np.sum(atom_types_array == "Oh")
    o_count = np.sum(atom_types_array == "O")
    print(
        f"Classified O atoms: {ob_count} Ob (Si-O-Si), {oh_count} Oh (Si-O-H), {o_count} O (Si-O)"
    )
    print(f"Found {len(bonds)} bonds and {len(angles)} angles")

    return atom_types, bonds, angles


def generate_silica_ndx(substr_path: str, build: bool = True):
    """
    Generate .ndx file for a silica substrate, based on the .gro file.
    Creates two groups: SILICA (all substrate atoms) and OH (hydrogen atoms from silanol groups).

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

    # Identify all atoms and OH groups
    silica_atoms = []
    oh_atoms = []

    # Identify hydrogen atoms (silanol groups)
    for i, atom in enumerate(substr.topology.atoms):
        if atom.name.startswith("H"):
            oh_atoms.append(i + 1)  # 1-indexed for GROMACS

    # All atoms in the substrate belong to SILICA group
    silica_atoms = list(range(1, oh_atoms[0]))  # 1-indexed for GROMACS

    # Write the index file
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        f.write("[ SILICA ]\n")
        # Write atoms in groups of 15 per line
        for i in range(0, len(silica_atoms), 15):
            f.write(" ".join(map(str, silica_atoms[i : i + 15])) + "\n")

        f.write("\n[ OH ]\n")
        for i in range(0, len(oh_atoms), 15):
            f.write(" ".join(map(str, oh_atoms[i : i + 15])) + "\n")

    print(
        f"NDX file with {len(silica_atoms)} SILICA atoms and {len(oh_atoms)} OH atoms created."
    )

    return output_name


def generate_silica_itp(substr_path: str, build: bool = True):
    """
    Generate .itp file for a silica substrate, based on the .gro file.

    Parameters
    ----------
    substr_path : str
        The path to the .gro file of the substrate.
    build : bool
        Whether to build the itp if it does not exist. Default is True.

    Returns
    -------
    str
        The path to the file of the generated .itp file.
    """
    # --- Load and validate substrate ---
    substr = md.load(substr_path)
    box_lengths = substr.unitcell_lengths[0]

    # Validate box dimensions
    assert np.all(box_lengths[:3] > 0), "Box lengths must be positive."
    if box_lengths.shape[0] > 3:
        assert np.allclose(box_lengths[3:], 0), (
            "Box should be orthorhombic (last 6 box components should be zero)."
        )

    # --- Setup paths ---
    folder, substr_filename = os.path.split(substr_path)
    folder = os.path.split(folder)[0]
    substr_name = os.path.splitext(substr_filename)[0]
    output_name = substr_name + ".itp"
    filename = "_".join(substr_name.split("_")[:-1])
    output_path = os.path.join(folder, "itp", output_name)

    # Check if the ITP already exists
    if build and os.path.isfile(output_path):
        print(f"A ready-made itp for substrate is used from `{folder}/itp`")
        return output_name
    else:
        print(
            "The itp does not exist for this substrate size. Forcibly generating it...."
        )

    # --- Atom type definitions ---
    metadata = {
        "Si": ["Si", 28.086, 2.10],  # Silicon
        "Ob": ["Ob", 15.999, -1.050],  # Bridging oxygen (Si-O-Si)
        "Oh": ["Oh", 15.999, -0.950],  # Hydroxyl oxygen (Si-O-H)
        "O": ["O", 15.999, -0.525],  # Surface oxygen (Si-O)
        # "O": ["O", 15.999, -0.6336207],  # Surface oxygen (Si-O)
        "H": ["H", 1.008, 0.4250],  # Hydrogen
    }

    # --- Generate topology elements ---
    atom_types, bonds, angles = classify_silica_atoms(substr)

    # --- Format ITP file sections ---
    # Atoms section
    atoms_text = ""
    total_charge = 0.0

    for i, atom_type in enumerate(atom_types):
        atom_name = substr.topology.atom(i).name
        type_name, mass, charge = metadata[atom_type]
        total_charge += charge
        atoms_text += "{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>10}{:>10}\n".format(
            i + 1, type_name, 1, "SIL", atom_name, 1, charge, mass
        )

    # Print total charge for verification
    print(f"Total substrate charge: {total_charge:.6f}")

    # Bonds section
    bonds_text = ""
    for i, j in bonds:
        bonds_text += "{:>8}{:>8}{:>8}\n".format(i + 1, j + 1, 1)

    # Angles section
    angles_text = ""
    for i, j, k in angles:
        angles_text += "{:>8}{:>8}{:>8}{:>8}\n".format(i + 1, j + 1, k + 1, 1)

    # --- Generate final ITP file ---
    # Read template
    with open(os.path.join(folder, filename + "_template" + ".itp")) as f:
        final_text = f.read().format(
            atoms=atoms_text,
            bonds=bonds_text,
            # angles=angles_text,
        )

    # Write output file
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        f.write(final_text)

    return output_name
