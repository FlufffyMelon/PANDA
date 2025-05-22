import numpy as np
import argparse
import sys
from tqdm import tqdm

from panda import profile_approx_from_array, get_each_density_profile, str2bool


def get_statistical_inefficiency(
    trajectory_file: str,
    topology_file: str,
    rho_bulk: float,
    l: float,
    phi: float,
    H: float,
    interface_type: str,
    residue: str,
    sl: int,
    max_block_length: float,
    begin_time: float,
    time: float,
    timestep: float,
    chunk_length: float = 1000,
    units: str = "ps",
    reverse: bool = False,
    display: bool = True,
):
    """
    Calculate the statistical inefficiency of the trajectory as a function of block length.

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
    axises, denses = get_each_density_profile(
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
        _, _, best_i = profile_approx_from_array(
            denses[i, :], axises[i, :], rho_bulk, l, phi, H, interface_type, extention='alpha', display=display
        )

        angles[i] = np.rad2deg(best_i["theta"])

    # Invert to calculate blocks in reverse order
    if reverse:
        angles = np.flip(angles)

    # Compute variance and mean
    v = np.var(angles, ddof=1)
    m = np.mean(angles)

    # Pre-allocate arrays
    block_sizes = np.zeros(max_block_length - 1, dtype=int)
    si = np.zeros(max_block_length - 1)

    # Calculate the scaled variance for each block length
    for t in range(2, max_block_length):
        nblocks = L * tau // t
        if nblocks > 1:
            xg = np.split(angles[: nblocks * t], nblocks)
            block_means = np.array([np.mean(block) for block in xg])
            v2 = np.sum((block_means - m) ** 2) / (nblocks - 1)
            si[t - 2] = t * v2 / v
            block_sizes[t - 2] = t
        else:
            break

    return block_sizes, si


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--trajectory_file", type=str, help="Path to the GROMACS trajectory (.xtc) file.")
    parser.add_argument("--topology_file", type=str, help="Path to the GROMACS topology (.gro) file.")
    parser.add_argument("--output_file", type=str, help="Path to the output file where results will be saved.")
    parser.add_argument("--rho_bulk", type=float, help="Bulk density of the droplet.")
    parser.add_argument("--l", type=float, help="Ratio of pore length to pore height.")
    parser.add_argument("--phi", type=float, help="Volume fraction of the droplet.")
    parser.add_argument("--H", type=float, help="Height of the pore.")
    parser.add_argument("--interface_type", type=str, help="Type of surface to use in approximation, e.g., 'roll'.")
    parser.add_argument("--residue", type=str, help="Residue type of the droplet.")
    parser.add_argument("--sl", type=int, help="Number of slices to divide the pore into.")
    parser.add_argument("--max_block_length", type=int, help="Maximum block length to calculate the statistical inefficiency.")
    parser.add_argument("--begin_time", type=float, help="Starting time for calculations.")
    parser.add_argument("--time", type=float, help="Total time of the simulation.")
    parser.add_argument("--timestep", type=float, help="Time step of the simulation.")
    parser.add_argument("--chunk_length", type=float, default=1000, help="Length of the chunk for dynamic load of trajectory (default: 1000).")
    parser.add_argument("--units", type=str, default="ps", choices=["ps", "ns"], help="Units of time in the simulation (default: 'ps').")
    parser.add_argument("--reverse", type=str2bool, default=False, help="If True, process trajectory in reverse order (default: False).")
    parser.add_argument("--display", type=str2bool, default=True, help="If True, display optimization details (default: True).")

    args = parser.parse_args()

    block_sizes, si = get_statistical_inefficiency(
        trajectory_file=args.trajectory_file,
        topology_file=args.topology_file,
        rho_bulk=args.rho_bulk,
        l=args.l,
        phi=args.phi,
        H=args.H,
        interface_type=args.interface_type,
        residue=args.residue,
        sl=args.sl,
        max_block_length=args.max_block_length,
        begin_time=args.begin_time,
        time=args.time,
        timestep=args.timestep,
        chunk_length=args.chunk_length,
        units=args.units,
        reverse=args.reverse,
        display=args.display,
    )

    # Write results to the output file using numpy
    output_data = np.column_stack((block_sizes, si))
    header = "Block Sizes\tStatistical Inefficiencies"
    np.savetxt(args.output_file, output_data, delimiter="\t", header=header, comments="")
