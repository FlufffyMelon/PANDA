import mdtraj as md
import numpy as np
from tqdm import tqdm
import argparse
import sys

from profile_approx import _profile_approx_from_array
from utils import get_numerical_density_profile, get_center_pbc, apply_pbc, str2bool

def extract_density_profiles(
    trajectory_file: str,
    topology_file: str,
    rho_bulk: float,
    l: float,
    phi: float,
    H: float,
    interface_type: str,
    residue: str,
    sl: int,
    block_length: float,
    begin_time: float,
    time: float,
    timestep: float,
    units: str = "ps",
    samples: int = 3,
    display: bool = True,
):
    # Validate inputs and initialize variables
    assert residue in ["DECAN"], "This residue type is not available"
    assert units.lower() in ["ns", "ps"], "Wrong units"
    u = 1000 if units.lower() == "ns" else 1

    # Calculate the number of frames and chunks
    N = int(time * u) // int(timestep * u) + 1
    L = int((time - begin_time) * u) // int(block_length * u)
    tau = int(block_length * u) // int(timestep * u)
    start_frame = N - tau * L

    t = np.zeros(L)
    angles = np.zeros(L)
    delta = np.zeros(L)

    # Iterate over trajectory chunks
    chunk_iter = md.iterload(
        trajectory_file, top=topology_file, chunk=tau, skip=start_frame
    )
    for chunk_idx, chunk in enumerate(tqdm(chunk_iter, total=L, desc="Chunk")):
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

        axis, dens = get_numerical_density_profile(positions, box, sl, center=True)

        t[chunk_idx] = chunk.time[-1]
        angles[chunk_idx] = np.rad2deg(best_i["theta"])
        delta[chunk_idx] = best_i["delta"]

    return t, angles, delta


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # Adding arguments (based on the function signature)
    parser.add_argument("--trajectory_file", type=str, required=True, help="Path to the trajectory file (e.g., .xtc).")
    parser.add_argument("--topology_file", type=str, required=True, help="Path to the topology file (e.g., .gro).")
    parser.add_argument("--output_file", type=str, required=True, help="Path to the output file.")
    parser.add_argument("--rho_bulk", type=float, required=True, help="Bulk density of the system.")
    parser.add_argument("--l", type=float, required=True, help="Characteristic length scale.")
    parser.add_argument("--phi", type=float, required=True, help="Packing fraction.")
    parser.add_argument("--H", type=float, required=True, help="Height of the system.")
    parser.add_argument("--interface_type", type=str, required=True, help="Type of interface (e.g., roll, flat, etc.).")
    parser.add_argument("--residue", type=str, required=True, help="Residue name in the simulation.")
    parser.add_argument("--sl", type=int, required=True, help="Number of slices.")
    parser.add_argument("--block_length", type=float, required=True, help="Maximum block length to analyze.")
    parser.add_argument("--begin_time", type=float, required=True, help="Start time of the trajectory to analyze.")
    parser.add_argument("--time", type=float, required=True, help="Total simulation time.")
    parser.add_argument("--timestep", type=float, required=True, help="Time step of the simulation.")
    parser.add_argument("--units", type=str, default="ps", choices=["ps", "ns"], help="Units of time in the simulation (default: 'ps').")
    parser.add_argument("--samples", type=int, default=3, help="Number of points in a grid for minimization process")
    parser.add_argument("--display", type=str2bool, default=True, help="If True, display optimization details (default: True).")

    args = parser.parse_args()

    # Call the function
    time, angles, delta = get_angle_delta_blocks(
        args.trajectory_file,
        args.topology_file,
        args.rho_bulk,
        args.l,
        args.phi,
        args.H,
        args.interface_type,
        args.residue,
        args.sl,
        args.block_length,
        args.begin_time,
        args.time,
        args.timestep,
        args.units,
        args.samples,
        args.display
    )

    # Save the output to a file
    output_data = np.column_stack((time, angles, delta))
    header = "Block Size\tAngle\tDelta"
    np.savetxt(args.output_file, output_data, delimiter="\t", header=header, comments="")
