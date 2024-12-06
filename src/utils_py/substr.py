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
    chars = string.digits + string.ascii_uppercase + string.ascii_lowercase
    base = len(chars)
    encoded = []

    while num > 0:
        encoded.append(chars[num % base])
        num //= base

    # Pad with zeros if necessary
    while len(encoded) < length:
        encoded.append(chars[0])

    return "".join(reversed(encoded))


def generate_substrate(
    unitcell_path: str, Lx: float, Ly: float, Lz: float, freeze_substr: bool = False
):
    unitcell = read_gro(unitcell_path)
    unitcell_box = unitcell.box

    Nx, Ny, Nz = np.clip(np.array([Lx, Ly, Lz]) / unitcell_box[:3], 1, None).astype(int)
    N = Nx * Ny * Nz
    N_atoms = len(unitcell.atoms)
    ex, ey, ez = get_box_vectors(unitcell_box)

    structure = Structure(
        title=f"Substrate {Nx}x{Ny}x{Nz}",
        box=unitcell_box[:3] * np.array([Nx, Ny, Nz]),
        atoms=np.empty(N * N_atoms, dtype=object),
        atoms_xyz=np.zeros((N * N_atoms, 3)),
    )

    print("Possible substrate size {:.1f}x{:.1f}x{:.1f}".format(*structure.box))

    # Naming variables
    folder, unitcell_filename = os.path.split(unitcell_path)
    unitcell_name = os.path.splitext(unitcell_filename)[0]
    filename = "_".join(unitcell_name.split("_")[:-1])
    output_name = f"{filename}_{Nx}x{Ny}x{Nz}.gro"
    output_path = os.path.join(folder, "gro", output_name)

    # if os.path.isfile(output_path):
    #     print(f"Such substrate already exist! Please check `{folder}/gro`")
    #     return output_name

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
    assert (len(box) == 3) or (len(box) == 9)

    if len(box) == 3:
        ex = np.array([box[0], 0, 0])
        ey = np.array([0, box[1], 0])
        ez = np.array([0, 0, box[2]])
    elif len(box) == 9:
        ex = np.array([box[0], box[3], box[4]])
        ey = np.array([box[5], box[1], box[6]])
        ez = np.array([box[7], box[8], box[2]])

    return (ex, ey, ez)


def generate_calcite_itp(substr_path: str):
    substr = read_gro(substr_path)
    # substr.box = substr.box[:3]
    assert len(substr.box) == 3, "Box should be orthogonal"

    folder, substr_filename = os.path.split(substr_path)
    folder = os.path.split(folder)[0]
    substr_name = os.path.splitext(substr_filename)[0]
    output_name = substr_name + ".itp"
    filename = "_".join(substr_name.split("_")[:-1])
    output_path = os.path.join(folder, "itp", output_name)

    # if os.path.isfile(output_path):
    #     print(f"Itp for substrate already exist! Please check `{folder}/itp`")
    #     return output_name

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
    with open(os.path.join(folder, "itp", substr_name + ".itp"), "w") as f:
        f.write(final_text)

    return output_name


def get_calcite_neighbors_list(substr: Structure):
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
                if len(O_neigh) > 3:
                    print("Too many neighbours...")
                    break

                label_j = substr.atoms[j]
                name_j = (
                    label_j.name[:2] if label_j.name[:2] == "Ca" else label_j.name[0]
                )

                if name_j == "O":
                    if np.linalg.norm(rij(i, j, substr)) <= l:
                        O_neigh.append(j)

            if len(O_neigh) != 3:
                print("Can`t find all neighbours!!!", len(O_neigh), i, j)
                exit()

            neigh_dict[i] = O_neigh.copy()

    return neigh_dict


def get_calcite_neighbors_list_numpy(substr: Structure):
    # substr.box = substr.box[:3]
    assert len(substr.box) == 3, "Box should be orthogonal"

    # Do not touch!!!!!!
    l = 0.118 + 0.022
    neigh_dict = dict()

    oxygen_mask = np.zeros(len(substr.atoms), dtype=bool)
    oxygen_real_ids = []
    carbon_ids = []
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
    for i in tqdm(carbon_ids):
        O_neigh = []

        rij = (substr.atoms_xyz - substr.atoms_xyz[i, :])[oxygen_mask]
        mask = np.abs(rij) >= substr.box / 2
        rij -= substr.box * mask * np.sign(rij)
        neigh_oxygen = np.argwhere(np.linalg.norm(rij, axis=1) < l).ravel()

        O_neigh = list(oxygen_real_ids[neigh_oxygen])

        if len(O_neigh) != 3:
            print("Incorrect number of neighbours!!!", len(O_neigh))
            exit()

        neigh_dict[i] = O_neigh.copy()

    return neigh_dict


def rij(i, j, structure):
    rij = structure.atoms_xyz[j, :] - structure.atoms_xyz[i, :]
    mask = np.abs(rij) >= structure.box / 2
    rij -= structure.box * mask * np.sign(rij)

    return rij


def angle(i, j, k, structure):
    ji = rij(j, i, structure)
    unit_ji = ji / np.linalg.norm(ji)
    jk = rij(j, k, structure)
    unit_jk = jk / np.linalg.norm(jk)

    dot_product = np.dot(unit_ji, unit_jk)
    angle = np.arccos(dot_product)

    return np.rad2deg(angle)
