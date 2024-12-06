import mdtraj as md
import numpy as np
from tqdm import tqdm

from .profile_approx import _profile_approx_alpha_from_array
from .utils import get_numerical_density_profile, get_center_pbc, apply_pbc


def get_statistical_inefficiency_RCA(
    trajectory_file: str,
    topology_file: str,
    rho_bulk: float,
    l: float,
    phi: float,
    H: float,
    interface_type: str,
    residue: str,
    sl: int,
    chunk_length: float,
    begin_time: float,
    time: float,
    timestep: float,
    max_block_length: float,
    units: str = "ps",
    display: bool = True,
):
    """
    Calculate the statistical inefficiency of the trajectory as a function of block length
    using RCA method.

    Parameters
    ----------
    trajectory_file : str
        Path to the GROMACS trajectory (.xtc) file.
    topology_file : str
        Path to the GROMACS topology (.gro) file.
    rho_bulk : float
        Bulk density of the droplet.
    l : float
        Ratio of pore length to pore height.
    phi : float
        Volume fraction of the droplet.
    H : float
        Height of the pore.
    interface_type : str
        Type of surface to use in approximation, e.g., 'roll'.
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
    max_block_length : float
        Maximum block length to calculate the statistical inefficiency.
    units : str, optional
        Units of time in the simulation. Default is 'ps'.
    display : bool, optional
        If True, display optimization details. Default is True.

    Returns
    -------
    block_sizes : np.array
        Array of block sizes.
    si : np.array
        Array of scaled variances.
    """
    assert units.lower() in ["ns", "ps"], "Wrong units"
    u = 1000 if units.lower() == "ns" else 1

    L = int((time - begin_time) * u) // int(chunk_length * u)
    tau = int(chunk_length * u) // int(timestep * u)

    # Calculate all density profiles
    print("Calculating all density profiles")
    axises, denses = get_each_density_profiles(
        trajectory_file=trajectory_file,
        topology_file=topology_file,
        residue=residue,
        sl=sl,
        chunk_length=chunk_length,
        begin_time=begin_time,
        time=time,
        timestep=timestep,
        units=units,
    )

    # Calculate contact angles
    print("Calculating contact angles")
    assert len(denses) == L * tau, "Length of denses array does not match with L * tau"
    angles = np.zeros(L * tau)
    for i in tqdm(range(L * tau)):
        _, _, best_i = _profile_approx_alpha_from_array(
            denses[i, :], axises[i, :], rho_bulk, l, phi, H, interface_type, display
        )

        angles[i] = np.rad2deg(best_i["theta"])

    # Invert to calculate blocks in reverse order
    angles = np.flip(angles)

    # Compute variance, mean, and length
    v = np.var(angles)
    m = np.mean(angles)

    # Pre-allocate arrays
    block_sizes = np.zeros(max_block_length - 1, dtype=int)
    si = np.zeros(max_block_length - 1)

    # Calculate the scaled variance for each block length
    for t in range(2, max_block_length):
        nblocks = L * tau // t
        if nblocks > 0:
            xg = np.split(angles[: nblocks * t], nblocks)
            block_means = np.array([np.mean(block) for block in xg])
            v2 = np.sum((block_means - m) ** 2) / nblocks
            si[t - 2] = t * v2 / v
            block_sizes[t - 2] = nblocks
        else:
            break

    return block_sizes, si


def get_each_density_profiles(
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
    for chunk_idx, chunk in enumerate(tqdm(chunk_iter, total=L, desc="Chunk:")):
        assert (
            int(timestep) == int(chunk.timestep)
        ), f"The input timestep and the actual timestep do not match. Perhaps timestep = {int(chunk.timestep)}?"

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
