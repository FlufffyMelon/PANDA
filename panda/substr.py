import numpy as np
import string
import os
import sys

sys.path.append("..")  # Avoid error with importing of src
from src.utils_py.io.gro import read_gro, write_gro
from src.utils_py.gro.Structure import Structure

from tqdm import tqdm
from itertools import combinations


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
    unitcell_path: str, Lx: float, Ly: float, Lz: float, freeze_substr: bool = False, build: bool = True
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

    Returns
    -------
    str
        The filename of the generated substrate.

    """
    unitcell = read_gro(unitcell_path)
    unitcell_box = unitcell.box

    Nx, Ny, Nz = np.round(
        np.clip(np.array([Lx, Ly, Lz]) / unitcell_box[:3], 1, None)
    ).astype(int)
    N = Nx * Ny * Nz
    N_atoms = len(unitcell.atoms)
    ex, ey, ez = get_box_vectors(unitcell_box)

    structure = Structure(
        title=f"Substrate {Nx}x{Ny}x{Nz}",
        box=unitcell_box[:3] * np.array([Nx, Ny, Nz]),
        atoms=np.empty(N * N_atoms, dtype=object),
        atoms_xyz=np.zeros((N * N_atoms, 3)),
    )

    print("Possible substrate size {:.1f}x{:.1f}x{:.1f} ({}x{}x{})".format(*structure.box, Nx, Ny, Nz))

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
            print(f"A substrate with this size does not exist. Forcibly generating it....")

    # Replicate unitcell
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                index = i * Ny * Nz + j * Nz + k

                for a in range(N_atoms):
                    atom_label = unitcell.atoms[a].copy()
                    atom_label.mol_id = index + 1
                    atom_label.id = index * N_atoms + a + 1

                    structure.atoms[index * N_atoms + a] = atom_label
                    structure.atoms_xyz[index * N_atoms + a, :] = (
                        unitcell.atoms_xyz[a, :] + ex * i + ey * j + ez * k
                    )

    if not freeze_substr:
        # Rename the atoms to be unique
        counter_dict = {"C": 0, "O": 0, "Ca": 0}
        for i, label in np.ndenumerate(structure.atoms):
            name = "".join([i for i in label.name if not i.isdigit()])

            counter_dict[name] += 1

            structure.atoms[i].name = name + base62_encode(
                counter_dict[name], length=(3 if name == "Ca" else 4)
            )

    with open(output_path, "w") as f:
        f.write(write_gro(structure.apply_pbc()))

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
    substr = read_gro(substr_path)
    # substr.box = substr.box[:3]
    assert len(substr.box) == 3, "Box should be orthogonal"

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
            print(f"The itp does not exist for this substrate size. Forcibly generating it....")

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


def get_calcite_neighbors_list(substr: Structure):
    """
    Generate a dictionary of neighbors for a given substrate Structure.

    Parameters
    ----------
    substr : Structure
        The substrate structure.

    Returns
    -------
    dict
        A dictionary of neighbors, where each key is a Carbon atom index and
        the value is a list of its Oxygen neighbors.

    Notes
    -----
    The substrate structure should be a calcite crystal, with the atoms labeled
    as "Ca" and "C" for Calcium and Carbon, respectively. The Oxygen atoms are
    labeled as "O". The box vectors should be orthogonal.

    """
    # substr.box = substr.box[:3]
    assert len(substr.box) == 3, "Box should be orthogonal"

    # Do not touch!!!!!!
    l = 0.118 + 0.022
    neigh_dict = dict()

    print("Generating neighbors list")
    for i in tqdm(range(len(substr.atoms))):
        label_i = substr.atoms[i]
        name_i = label_i.name[:2] if label_i.name[:2] == "Ca" else label_i.name[0]

        if name_i == "C":
            O_neigh = []
            for j in range(len(substr.atoms)):
                # Check if the neighbor is an Oxygen atom
                if len(O_neigh) > 3:
                    print("Too many neighbours...")
                    break

                label_j = substr.atoms[j]
                name_j = (
                    label_j.name[:2] if label_j.name[:2] == "Ca" else label_j.name[0]
                )

                if name_j == "O":
                    # Check if the distance between the Carbon atom and the Oxygen atom is within the cutoff
                    if np.linalg.norm(rij(i, j, substr)) <= l:
                        O_neigh.append(j)

            if len(O_neigh) != 3:
                print("Can`t find all neighbours!!!", len(O_neigh), i, j)
                exit()

            neigh_dict[i] = O_neigh.copy()

    return neigh_dict


def get_calcite_neighbors_list_numpy(substr: Structure):
    """
    Generate a dictionary of calcite atom neighbors.

    For each carbon atom in the structure, find all neighboring oxygen atoms
    within a specified cutoff distance and store them in a dictionary.

    Parameters
    ----------
    substr : Structure
        A Structure object containing atoms and their positions.

    Returns
    -------
    dict
        A dictionary where keys are indices of carbon atoms and values are lists
        of indices of neighboring oxygen atoms.

    """
    assert len(substr.box) == 3, "Box should be orthogonal"

    # Cutoff distance for determining neighboring atoms
    l = 0.118 + 0.022
    neigh_dict = dict()

    # Initialize masks and lists for oxygen and carbon atom indices
    oxygen_mask = np.zeros(len(substr.atoms), dtype=bool)
    oxygen_real_ids = []
    carbon_ids = []

    # Identify indices of oxygen and carbon atoms
    for i in range(len(substr.atoms)):
        name = (
            substr.atoms[i].name[:2]
            if substr.atoms[i].name[:2] == "Ca"
            else substr.atoms[i].name[0]
        )

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
        rij = (substr.atoms_xyz - substr.atoms_xyz[i, :])[oxygen_mask]
        mask = np.abs(rij) >= substr.box / 2
        rij -= substr.box * mask * np.sign(rij)

        # Determine indices of oxygen atoms within the cutoff distance
        neigh_oxygen = np.argwhere(np.linalg.norm(rij, axis=1) < l).ravel()
        O_neigh = list(oxygen_real_ids[neigh_oxygen])

        # Ensure exactly 3 neighboring oxygen atoms are found
        assert len(O_neigh) == 3, f"Incorrect number of neighbours ({len(O_neigh)})!!!"

        # Store the neighbors in the dictionary
        neigh_dict[i] = O_neigh.copy()

    return neigh_dict


def rij(i, j, structure):
    """
    Calculate the relative position vector between atoms i and j taking into account PBC.

    Parameters
    ----------
    i : int
        Index of the first atom.
    j : int
        Index of the second atom.
    structure : Structure
        The structure object containing the atoms and their positions.

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
    rij = structure.atoms_xyz[j, :] - structure.atoms_xyz[i, :]
    mask = np.abs(rij) >= structure.box / 2
    rij -= structure.box * mask * np.sign(rij)

    return rij


def angle(i, j, k, structure):
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
    structure : Structure
        The structure object containing the atoms and their positions.

    Returns
    -------
    float
        The angle in degrees between the vectors defined by atoms i-j and j-k.

    Notes
    -----
    The calculation takes into account the periodic boundary conditions of the
    simulation box.

    """
    # Calculate the relative position vectors between atoms i-j and j-k
    ji = rij(j, i, structure)
    unit_ji = ji / np.linalg.norm(ji)
    jk = rij(j, k, structure)
    unit_jk = jk / np.linalg.norm(jk)

    # Calculate the dot product between the two vectors
    dot_product = np.dot(unit_ji, unit_jk)

    # Calculate the angle from the dot product
    angle = np.arccos(dot_product)

    return np.rad2deg(angle)
