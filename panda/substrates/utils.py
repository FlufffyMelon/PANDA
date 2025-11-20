import numpy as np
import string
import os
import ase
import mdtraj as md
from mdtraj.core.topology import Topology
from tqdm import tqdm
from itertools import combinations
from panda.utils import apply_pbc
from abc import ABC


class CustomSubstrate(ABC):
    """
    Abstract base class for substrate generators.

    Attributes
    ----------
    gro_path : str
        Path to the generated substrate gro file
    itp_path : str
        Path to the generated substrate itp file
    ndx_path : str
        Path to the generated substrate ndx file
    """

    pass


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


# def generate_substrate(
#     unitcell: md.Trajectory, Lx: float, Ly: float, Lz: float, residue_name: str
# ):
#     """
#     Generate a substrate from a given unitcell.

#     Parameters
#     ----------
#     unitcell : md.Trajectory
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

#     # unitcell = md.load(unitcell_path, top=unitcell_path)
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

#     # Prepare arrays for new coordinates and atom info
#     all_xyz = np.zeros((N * N_atoms, 3))
#     atomnames = []
#     mol_ids = []

#     # Replicate unitcell
#     for i in range(Nx):
#         for j in range(Ny):
#             for k in range(Nz):
#                 index = i * Ny * Nz + j * Nz + k
#                 offset = ex * i + ey * j + ez * k
#                 for a in range(N_atoms):
#                     idx = index * N_atoms + a
#                     all_xyz[idx, :] = unitcell.xyz[0, a, :] + offset
#                     atom = unitcell.topology.atom(a)
#                     atomnames.append(atom.name)
#                     mol_ids.append(index + 1)

#     # Build new topology
#     top = Topology()
#     residue = top.add_residue(residue_name, top.add_chain())
#     for i in range(N * N_atoms):
#         top.add_atom(atomnames[i], element=None, residue=residue, serial=i + 1)

#     # Set box
#     box = unitcell_box[:3] * np.array([Nx, Ny, Nz])
#     all_xyz = apply_pbc(all_xyz, box)
#     traj = md.Trajectory(
#         all_xyz.reshape(1, -1, 3),
#         top,
#         unitcell_lengths=box.reshape(1, 3),
#         unitcell_angles=np.array([[90.0, 90.0, 90.0]]),
#     )

#     print("Substrate successfully created!")
#     return traj


def get_box_vectors(box: np.array):
    assert len(box) in {3, 9}

    if len(box) == 3:
        ex, ey, ez = np.diag(box)
    else:
        ex, ey, ez = box[:3], box[3:6], box[6:]

    return ex, ey, ez


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


def ase2mdtraj(atoms: ase.Atoms, residue_name: str = "CRY"):
    """
    Convert ASE Atoms object to MDTraj Trajectory.

    Parameters
    ----------
    atoms : ase.Atoms
        Input ASE Atoms object.

    Returns
    -------
    traj : mdtraj.Trajectory
        MDTraj trajectory containing a single frame, with
        correct topology and unit cell.
    """
    # --- Build MDTraj Topology ---
    top = Topology()
    chain = top.add_chain()
    res = top.add_residue(residue_name, chain)  # single residue called CRY

    md_atoms = []
    for i, atom in enumerate(atoms):
        # use chemical symbol as element
        md_atom = top.add_atom(atom.symbol, element=None, residue=res, serial=i + 1)
        md_atoms.append(md_atom)

    # --- Coordinates ---
    # ASE positions are in Å; MDTraj expects nm
    coords = atoms.get_positions() / 10.0  # shape (N,3), units: nm
    coords = np.array([coords])  # add frame dimension

    # --- Unit cell ---
    # ASE gives cell in Å; MDTraj expects nm
    # cell = atoms.cell / 10.0  # a, b, c in nm
    cell = atoms.cell.lengths() / 10.0  # a, b, c in nm
    angles = atoms.cell.angles()  # alpha, beta, gamma in degrees

    # unitcell_lengths = np.array([cell[0, 0], cell[1, 1], cell[2, 2], cell[0, 1], cell[0, 2], cell[1, 0], cell[1, 2], cell[2, 0], cell[2, 1]])
    unitcell_lengths = np.array([cell])
    unitcell_angles = np.array([angles])

    # --- Build Trajectory ---
    traj = md.Trajectory(
        xyz=coords.reshape(1, -1, 3),
        topology=top,
        unitcell_lengths=unitcell_lengths,
        unitcell_angles=unitcell_angles,
        # unitcell_vectors=cell.reshape(1, 3, 3)
    )

    return traj
