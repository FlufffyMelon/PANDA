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
    assert units.lower() in ["ns", "ps"], "Wrong units"
    u = 1000 if units.lower() == "ns" else 1
    
    L = int((time - begin_time) * u) // int(chunk_length * u)
    tau = int(chunk_length * u) // int(timestep * u)
    
    # Calculating all density profiles
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
    angles = []

    print("Calculating contact angles")
    assert len(denses) == L * tau, "Length of denses array does not match with L * tau"
    for i in tqdm(range(L * tau)):
        _, _, best_i = _profile_approx_alpha_from_array(
            denses[i, :], axises[i, :], rho_bulk, l, phi, H, interface_type, display
        )

        angles.append(np.rad2deg(best_i["theta"]))

    # Invert to calculate blocks in reverse order
    angles = np.array(angles[::-1])

    # Compute variance, mean, and length
    v = np.var(angles, ddof=1)  # ddof=1 for sample variance (like in R)
    m = np.mean(angles)
    n = len(angles)

    block_sizes = []
    si = []  # List to store the scaled variances
    for t in range(2, max_block_length):
        nblocks = n // t  # Number of blocks
        if nblocks > 0:
            # Create blocks of size `t`
            xg = np.split(angles[:nblocks * t], nblocks)
            # Calculate mean of each block
            block_means = np.array([np.mean(block) for block in xg])
            # Variance of the block means
            v2 = np.sum((block_means - m) ** 2) / nblocks
            # Compute scaled variance and append
            si.append(t * v2 / v)
            block_sizes.append(nblocks)
        else:
            break

    block_sizes = np.array(block_sizes)
    si = np.array(si)
    
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
    # TODO: Add better residue checker
    assert residue in ["DECAN"], "This residue type is not avaliable"
    assert units.lower() in ["ns", "ps"], "Wrong units"
    u = 1000 if units.lower() == "ns" else 1

    N = int(time * u) // int(timestep * u) + 1
    L = int((time - begin_time) * u) // int(chunk_length * u)
    tau = int(chunk_length * u) // int(timestep * u)
    start_frame = N - tau * L

    axises = []
    denses = []
    
    for chunk in tqdm(
        md.iterload(trajectory_file, top=topology_file, chunk=tau, skip=start_frame),
        total=L,
        desc="Chunk:"
    ):
        assert (
            int(timestep) == int(chunk.timestep)
        ), f"The input timestep and the actual timestep do not match. Perhaps timestep = {int(chunk.timestep)}?"

        residue_mask = chunk.top.select(f"resname {residue}")
        positions = chunk.xyz[:, residue_mask, :]
        box = chunk.unitcell_lengths[:, np.newaxis, :]

        # Calculating and centering aroud center mass relative to the PBC
        center = get_center_pbc(positions, box)
        positions -= center
        positions += box / 2
        positions = apply_pbc(positions, box)

        # Calculating profiles for each frame
        for i in range(tau):
            axis_i, dens_i = get_numerical_density_profile(
                positions[i, :, :], box[i, :, :], sl, center=True
            )

            axises.append(axis_i)
            denses.append(dens_i)

    axises = np.array(axises)
    denses = np.array(denses)

    return axises, denses
