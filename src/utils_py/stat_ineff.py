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
    block_lenght: float,
    begin_time: float,
    time: float,
    timestep: float,
    units: str = "ps",
    display: bool = True,
):
    # TODO: Add better residue checker
    assert residue in ["DECAN"], "This residue type is not avaliable"
    assert units.lower() in ["ns", "ps"], "Wrong units format"
    units = 1000 if units.lower() == "ns" else 1

    N = int(time * units) // int(timestep * units) + 1
    L = int((time - begin_time) * units) // int(block_lenght * units)
    tau = int(block_lenght * units) // int(timestep * units)
    start_frame = N - tau * L

    X = []
    Y = []

    for chunk in tqdm(
        md.iterload(trajectory_file, top=topology_file, chunk=tau, skip=start_frame),
        total=L,
    ):
        assert (
            int(timestep) == int(chunk.timestep)
        ), f"The input timestep and the actual timestep do not match. Perhaps timestep = {int(chunk.timestep)}?"

        residue_mask = chunk.top.select(f"resname {residue}")
        positions = chunk.xyz[:, residue_mask, :]
        box = chunk.unitcell_lengths[:, np.newaxis, :]

        # Calculating and centering aroud center mass with relative to the PBC
        center = get_center_pbc(positions, box)
        positions -= center
        positions += box / 2
        positions = apply_pbc(positions, box)

        # Calculating contact angle for each frame
        for i in range(tau):
            axis_i, dens_i = get_numerical_density_profile(
                positions[i, :, :], box[i, :, :], sl, center=True
            )
            _, _, best_i = _profile_approx_alpha_from_array(
                dens_i, axis_i, rho_bulk, l, phi, H, interface_type, display
            )

            X.append(np.rad2deg(best_i["theta"]))

        # # Calculating contact angls for block average
        # axis, dens = get_numerical_density_profile(positions, box, sl, center=True)
        # _, _, best = _profile_approx_alpha_from_array(
        #     dens, axis, rho_bulk, l, phi, H, interface_type, display
        # )

        # Y.append(np.rad2deg(best["theta"]))

    X = np.array(X)
    Y = np.array(Y)
    print("X", X, np.std(X))
    print("Y", Y, np.std(Y))

    return tau * np.std(Y) / np.std(X)
