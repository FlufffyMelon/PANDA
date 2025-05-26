import numpy as np
import string
import os
import mdtraj as md
from mdtraj.core.topology import Topology
from tqdm import tqdm
from itertools import combinations
from .utils import apply_pbc


def base62_encode(num, length=4):
    """
    Encode a given number in base62.

    Parameters
    ----------
    num : int
        The number to be encoded.
    length : int, optional
        The desired length of the encoded string. If the
        number is shorter than this, it will be padded with
        zeros. The default is 4.

    Returns
    -------
    str
        The encoded string.

    """
    chars = string.digits + string.ascii_uppercase + string.ascii_lowercase
    base = len(chars)
    encoded = []

    # Repeatedly divide the number by the base and use the remainder
    # to construct the encoded string.
    q, r = divmod(num, base)
    encoded.append(chars[r])
    while q > 0:
        q, r = divmod(q, base)
        encoded.append(chars[r])

    # Pad with zeros if necessary
    encoded.extend([chars[0]] * (length - len(encoded)))

    return "".join(reversed(encoded))


def generate_substrate(
    unitcell_path: str,
    Lx: float,
    Ly: float,
    Lz: float,
    freeze_substr: bool = False,
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
        The filename of the generated substrate.
    """

    unitcell = md.load(unitcell_path, top=unitcell_path)
    unitcell_box = unitcell.unitcell_lengths[0]

    Nx, Ny, Nz = np.round(
        np.clip(np.array([Lx, Ly, Lz]) / unitcell_box[:3], 1, None)
    ).astype(int)
    N = Nx * Ny * Nz
    N_atoms = unitcell.n_atoms
    ex, ey, ez = get_box_vectors(unitcell_box)

    # Naming variables
    folder, unitcell_filename = os.path.split(unitcell_path)
    unitcell_name = os.path.splitext(unitcell_filename)[0]
    filename = "_".join(unitcell_name.split("_")[:-1])
    output_name = f"{filename}_{Nx}x{Ny}x{Nz}.gro"
    output_path = os.path.join(folder, "gro", output_name)

    if not build:
        if os.path.isfile(output_path):
            print(f"A ready-made substrate is used from `{folder}/gro`")
            return output_name
        else:
            print(
                "A substrate with this size does not exist. Forcibly generating it...."
            )

    # Prepare arrays for new coordinates and atom info
    all_xyz = np.zeros((N * N_atoms, 3))
    resnames = []
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
                    resnames.append(atom.residue.name)
                    atomnames.append(atom.name)
                    mol_ids.append(index + 1)

    # Optionally rename atoms to be unique
    if not freeze_substr:
        counter_dict = {"C": 0, "O": 0, "Ca": 0}
        for idx, name in enumerate(atomnames):
            base = "".join([c for c in name if not c.isdigit()])
            counter_dict[base] += 1
            atomnames[idx] = base + base62_encode(
                counter_dict[base], length=(3 if base == "Ca" else 4)
            )

    # Build new topology
    top = Topology()
    chain = top.add_chain()
    residue_map = {}
    for i in range(N * N_atoms):
        resname = resnames[i]
        if (mol_ids[i], resname) not in residue_map:
            residue = top.add_residue(resname, chain)
            residue_map[(mol_ids[i], resname)] = residue
        else:
            residue = residue_map[(mol_ids[i], resname)]
        top.add_atom(atomnames[i], element=None, residue=residue)

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
    traj.save_gro(output_path)
    print("Substrate successfully created!")
    return output_name


def get_box_vectors(box: np.array):
    assert len(box) in {3, 9}

    if len(box) == 3:
        ex, ey, ez = np.diag(box)
    else:
        ex, ey, ez = box[:3], box[3:6], box[6:]

    return ex, ey, ez


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
        The filename of the generated .itp file.

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

    if not build:
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
    for i, label in np.ndenumerate(substr.atoms):
        # name = "".join([i for i in label.name if not i.isdigit()])
        name = label.name[:2] if label.name[:2] == "Ca" else label.name[0]
        counter_dict[name] += 1
        type, mass, charge = metadata[name]

        atoms_text += "{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>10}{:>10}\n".format(
            i[0] + 1, type, 1, "CAL", label.name, 1, charge, mass
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
    with open(output_path, "w") as f:
        f.write(final_text)

    return output_name


def get_calcite_neighbors_list(substr: md.Trajectory):
    """
    Generate a dictionary of neighbors for a given substrate using mdtraj.Trajectory.

    Parameters
    ----------
    substr : mdtraj.Trajectory
        The substrate structure as an mdtraj Trajectory object.

    Returns
    -------
    dict
        A dictionary of neighbors, where each key is a Carbon atom index and
        the value is a list of its Oxygen neighbors (indices).

    Notes
    -----
    The substrate structure should be a calcite crystal, with the atoms labeled
    as "Ca" and "C" for Calcium and Carbon, respectively. The Oxygen atoms are
    labeled as "O". The box vectors should be orthogonal.
    """
    assert len(substr.unitcell_lengths[0]) == 3, "Box should be orthogonal"
    l = 0.118 + 0.022
    neigh_dict = dict()
    xyz = substr.xyz[0]
    top = substr.topology
    atom_names = [a.name for a in top.atoms]
    # Identify indices of C and O atoms
    carbon_ids = [i for i, a in enumerate(atom_names) if a.startswith("C")]
    oxygen_ids = [i for i, a in enumerate(atom_names) if a.startswith("O")]
    print("Generating neighbors list")
    for i in tqdm(carbon_ids):
        O_neigh = []
        for j in oxygen_ids:
            if len(O_neigh) > 3:
                print("Too many neighbours...")
                break
            if np.linalg.norm(rij(i, j, xyz, substr.unitcell_lengths[0])) <= l:
                O_neigh.append(j)
        if len(O_neigh) != 3:
            print("Can`t find all neighbours!!!", len(O_neigh), i)
            exit()
        neigh_dict[i] = O_neigh.copy()
    return neigh_dict


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
    l = 0.118 + 0.022
    neigh_dict = dict()
    xyz = substr.xyz[0]
    top = substr.topology
    atom_names = [a.name for a in top.atoms]
    oxygen_mask = np.array([name.startswith("O") for name in atom_names])
    oxygen_real_ids = np.where(oxygen_mask)[0]
    carbon_ids = [i for i, name in enumerate(atom_names) if name.startswith("C")]
    print("Generating neighbors list")
    for i in tqdm(carbon_ids):
        rij_vecs = xyz[oxygen_mask] - xyz[i, :]
        mask = np.abs(rij_vecs) >= substr.unitcell_lengths[0] / 2
        rij_vecs -= substr.unitcell_lengths[0] * mask * np.sign(rij_vecs)
        neigh_oxygen = np.argwhere(np.linalg.norm(rij_vecs, axis=1) < l).ravel()
        O_neigh = list(oxygen_real_ids[neigh_oxygen])
        assert len(O_neigh) == 3, f"Incorrect number of neighbours ({len(O_neigh)})!!!"
        neigh_dict[i] = O_neigh.copy()
    return neigh_dict


def rij(i, j, xyz, box):
    """
    Calculate the relative position vector between atoms i and j taking into account PBC.

    Parameters
    ----------
    i : int
        Index of the first atom.
    j : int
        Index of the second atom.
    xyz : np.ndarray
        Array of atomic coordinates (shape: n_atoms x 3).
    box : np.ndarray
        Box dimensions (length 3).

    Returns
    -------
    numpy.array
        The relative position vector between atoms i and j.

    Notes
    -----
    The calculation takes into account the periodic boundary conditions of the
    simulation box. If the relative position vector is larger than half the box
    size in any dimension, the box size is subtracted from the relative position
    vector to 'wrap' it around to the other side of the box.
    """
    rij = xyz[j, :] - xyz[i, :]
    mask = np.abs(rij) >= box / 2
    rij -= box * mask * np.sign(rij)
    return rij


def angle(i, j, k, xyz, box):
    """
    Calculate the angle between the vectors defined by atoms i-j and j-k.

    Parameters
    ----------
    i : int
        Index of the first atom.
    j : int
        Index of the second atom.
    k : int
        Index of the third atom.
    xyz : np.ndarray
        Array of atomic coordinates (shape: n_atoms x 3).
    box : np.ndarray
        Box dimensions (length 3).

    Returns
    -------
    float
        The angle in degrees between the vectors defined by atoms i-j and j-k.

    Notes
    -----
    The calculation takes into account the periodic boundary conditions of the
    simulation box.
    """
    ji = rij(j, i, xyz, box)
    unit_ji = ji / np.linalg.norm(ji)
    jk = rij(j, k, xyz, box)
    unit_jk = jk / np.linalg.norm(jk)
    dot_product = np.dot(unit_ji, unit_jk)
    angle = np.arccos(dot_product)
    return np.rad2deg(angle)
