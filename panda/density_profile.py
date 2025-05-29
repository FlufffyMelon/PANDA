import numpy as np
import mdtraj as md
from tqdm import tqdm
from panda.utils import (
    get_center_pbc,
    apply_pbc,
    validate_positions_and_box,
)


def get_numerical_density_profile(
    positions: np.array, box: np.array, sl: int, axis: str = "Z", center: bool = True
):
    # TODO: improve it
    assert axis == "Z", "Axis can only be 'X', 'Y' and 'Z'"
    positions, box = validate_positions_and_box(positions, box)

    axis, dz = np.linspace(
        0, np.mean(box, axis=0)[0, 2], sl + 1, endpoint=True, retstep=True
    )

    indexes = np.floor(positions[:, :, 2] / dz).astype(int)
    assert np.max(indexes) < sl, ""
    density_profiles = np.array([np.bincount(i, minlength=sl) for i in indexes])
    density_profiles = np.mean(density_profiles, axis=0)
    density_profiles /= np.mean(box, axis=0)[0, 0] * np.mean(box, axis=0)[0, 1] * dz

    if center:
        axis -= axis[-1] / 2

    return axis[:-1], density_profiles


def get_density_profile(
    trajectory_file: str,
    topology_file: str,
    residue: str,
    sl: int,
    chunk_length: float,
    begin_time: float,
    time: float,
    timestep: float,
    units: str = "ps",
):
    """
    Calculate all density profiles for given trajectory and topology files.

    Parameters
    ----------
    trajectory_file : str
        Path to the GROMACS trajectory (.xtc) file.
    topology_file : str
        Path to the GROMACS topology (.gro) file.
    residue : str
        Residue type of the droplet.
    sl : int
        Number of slices to divide the pore into.
    chunk_length : float
        Length of the chunk for dynamic load of trajectory.
    begin_time : float
        Starting time for calculations.
    time : float
        Total time of the simulation.
    timestep : float
        Time step of the simulation.
    units : str, optional
        Units of time in the simulation. Default is 'ps'.

    Returns
    -------
    axises : np.ndarray
        Array of axises of the profiles.
    denses : np.ndarray
        Array of density profiles.
    """

    # Validate inputs and initialize variables
    assert residue in ["DECAN"], "This residue type is not available"
    assert units.lower() in ["ns", "ps"], "Wrong units"
    u = 1000 if units.lower() == "ns" else 1

    # Calculate the number of frames and chunks
    N = int(time * u) // int(timestep * u) + 1
    L = int((time - begin_time) * u) // int(chunk_length * u)
    tau = int(chunk_length * u) // int(timestep * u)
    start_frame = N - tau * L

    # Pre-allocate arrays for results
    axises = np.empty((L * tau, sl))
    denses = np.empty((L * tau, sl))

    # Iterate over trajectory chunks
    chunk_iter = md.iterload(
        trajectory_file, top=topology_file, chunk=tau, skip=start_frame
    )
    for chunk_idx, chunk in enumerate(tqdm(chunk_iter, total=L, desc="Chunk")):
        assert int(timestep) == int(chunk.timestep), (
            f"The input timestep and the actual timestep do not match. Perhaps timestep = {int(chunk.timestep)}?"
        )

        # Select residue positions and apply periodic boundary conditions
        residue_mask = chunk.top.select(f"resname {residue}")
        positions = chunk.xyz[:, residue_mask, :]
        box = chunk.unitcell_lengths[:, np.newaxis, :]

        center = get_center_pbc(positions, box)
        positions -= center
        positions += box / 2
        positions = apply_pbc(positions, box)

        # Calculate density profiles for each frame
        for i in range(tau):
            axis_i, dens_i = get_numerical_density_profile(
                positions[i, :, :], box[i, :, :], sl, center=True
            )
            axises[chunk_idx * tau + i] = axis_i
            denses[chunk_idx * tau + i] = dens_i

    return axises, denses


def block_average_density_profile(
    axises: np.ndarray, denses: np.ndarray, block_size: int = 10
):
    """
    Compute block-averaged density profile and axis.

    Parameters
    ----------
    axises : np.ndarray
        Array of axis values (shape: [N, sl]), as from get_each_density_profile.
    denses : np.ndarray
        Array of density profiles (shape: [N, sl]), as from get_each_density_profile.
    block_size : int, optional
        Number of profiles to average in each block. Default is 10.

    Returns
    -------
    axis_avg : np.ndarray
        Block-averaged axis (shape: [n_blocks, sl]).
    dens_avg : np.ndarray
        Block-averaged density profile (shape: [n_blocks, sl]).
    """
    blocks_num = axises.shape[0] // block_size
    if blocks_num == 0:
        raise ValueError("Not enough profiles for the given block size.")

    mean_axises = np.zeros((blocks_num, axises.shape[1]))
    mean_denses = np.zeros((blocks_num, denses.shape[1]))

    for i in range(blocks_num):
        mean_axises[i, :] = np.mean(
            axises[(i * block_size) : ((i + 1) * block_size), :], axis=0
        )
        mean_denses[i, :] = np.mean(
            denses[(i * block_size) : ((i + 1) * block_size), :], axis=0
        )

    return mean_axises, mean_denses
