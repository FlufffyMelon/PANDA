import numpy as np
import sys

# sys.path.append("..")  # Avoid error with importing of src
# print(sys.path)
from src.utils_py.io.gro import read_gro
import src.utils_py.interface

# print("In module products sys.path[0], __package__ ==", sys.path[0], __package__)
# from .io.gro import read_gro
# import auxil
from scipy.optimize import minimize, Bounds


def L1(x, z, dens, l, phi, rho_func, grad_rho_func=None):
    """
    Function to minimize for finding the best angles theta and delta.

    Parameters:
        x (list): List of parameters to optimize, [theta, delta].

    Returns:
        float: Sum of absolute differences between the calculated and the experimental densities.
    """
    return np.sum(np.abs(rho_func(z, l, phi, *x) - dens))


def MAE(x, z, dens, l, phi, rho_func, grad_rho_func=None):
    return np.mean(np.abs(rho_func(z, l, phi, *x) - dens))


def RMSE(x, z, dens, l, phi, rho_func, grad_rho_func=None):
    return np.sqrt(np.mean((rho_func(z, l, phi, *x) - dens) ** 2))


def grad_L1(x, z, dens, l, phi, rho_func, grad_rho_func):
    """
    Gradient of a function to minimize for finding the best angles theta and delta.

    Parameters:
        x (list): List of parameters to optimize, [theta, delta].

    Returns:
        float: Sum of absolute differences between the calculated and the experimental densities.
    """
    return np.sum(
        np.sign(rho_func(z, l, phi, *x) - dens) * grad_rho_func(z, l, phi, *x), axis=1
    )


def profile_approx(
    dens_profile_file: str,
    rho_bulk: float,
    l: float,
    phi: float,
    H: float,
    interface_type: str,
    extention: str = "alpha",
    display: bool = True,
) -> tuple[np.array, np.array, dict]:
    """
    Approximate a density profile using a generalized algorithm taking into account theta and delta.

    Parameters:
        dens_profile_file (str): Path to the density profile file with an .xvg extension.
        rho_bulk (float): Bulk density of the droplet.
        l (float): Ratio of pore length to pore height.
        phi (float): Volume fraction of the droplet.
        H (float): Height of the pore.
        interface_type (str): Type of surface to use in approximation, e.g., 'roll'.
        display (bool, optional): If True, display optimization details. Default is True.

    Returns:
        tuple:
            - z (np.array): Normalized z-coordinates of the density profile.
            - dens (np.array): Normalized density profile values.
            - best (dict): Best-fit parameters, including 'theta' angle in radians.
    """

    # Check the file extension
    assert dens_profile_file.endswith("xvg"), (
        "The density profile file must have the extension .xvg"
    )

    # Read the density profile file
    z, dens = np.loadtxt(dens_profile_file, comments=["@", "#"], unpack=True)

    # Call the function to approximate the density profile
    return _profile_approx_from_array(
        dens,
        z,
        rho_bulk,
        l,
        phi,
        H,
        interface_type,
        extention=extention,
        display=display,
    )


def profile_approx_from_array(
    dens: np.array,
    z: np.array,
    rho_bulk: float,
    l: float,
    phi: float,
    H: float,
    interface_type: str,
    extention: str = "alpha",
    samples: int = 3,
    method="L-BFGS-B",
    x0: tuple = (3 * np.pi / 4, 0.1),
    callback=None,
    display: bool = True,
) -> tuple[np.array, np.array, dict]:
    """
    Approximate a density profile using a generalized algorithm considering theta and delta.

    Parameters:
        dens (np.array): Density profile values.
        z (np.array): z-coordinates of the density profile.
        rho_bulk (float): Bulk density of the droplet.
        l (float): Ratio of pore length to pore height.
        phi (float): Volume fraction of the droplet.
        H (float): Height of the pore.
        interface_type (str): Type of surface to use in approximation, e.g., 'roll'.
        display (bool, optional): If True, display optimization details. Default is True.

    Returns:
        tuple: Contains normalized z-coordinates, normalized density profile values,
               and best-fit parameters including 'theta' and 'delta'.
    """
    dens = dens.copy() / rho_bulk  # Normalize the density profile
    z = z.copy() / H  # Normalize the z-coordinates

    left = np.argmax(z > -0.5) - 1  # Find the left boundary of the droplet
    right = np.argmin(z < 0.5)  # Find the right boundary of the droplet

    interface_type = interface_type.lower()
    extention = extention.lower()
    # Create a function based on the surface type
    assert interface_type in [
        "droplet",
        "doughnut",
        "worm",
        "roll",
        "perforation",
        "layer",
    ], (
        "There is no such type of surface. Could be droplet, doughnut, worm, roll, perforation or layer"
    )

    assert extention in ["theta", "delta", "alpha"]
    rho = getattr(src.utils_py.interface, f"rho_{interface_type}_{extention}")
    try:
        grad_rho = getattr(
            src.utils_py.interface, f"grad_rho_{interface_type}_{extention}"
        )
    except AttributeError:
        # print("Gradient function does not exist. Using numerical approximation")
        grad_rho = None

    # Initial guess for theta and delta
    if samples == 1:
        if extention in ["alpha", "theta"]:
            theta0 = x0[0]
        if extention in ["alpha", "delta"]:
            delta0 = x0[1]
    else:
        if extention in ["alpha", "theta"]:
            theta0 = np.linspace(np.pi / 2 + 0.02, np.pi, samples)
        if extention in ["alpha", "delta"]:
            delta0 = np.linspace(0, 0.2, samples)

    if extention == "alpha":
        grid1, grid2 = np.meshgrid(theta0, delta0)
        x0_mesh = np.column_stack((grid1.ravel(), grid2.ravel()))
        bounds = ((np.pi / 2 + 0.01, np.pi), (0, 0.49))
    elif extention == "theta":
        x0_mesh = theta0
        bounds = ((np.pi / 2 + 0.01, np.pi),)
    elif extention == "delta":
        x0_mesh = delta0
        bounds = ((0, 0.49),)

    min_best = {
        "theta": np.pi,
        "delta": 0,
        "fun": np.inf,
        "mae": np.inf,
        "rmse": np.inf,
    }
    for x0i in x0_mesh:
        res = minimize(
            L1,
            x0i,
            method=method,
            # jac=None if grad_rho is None else grad_L1,
            bounds=bounds,
            args=(z[left:right], dens[left:right], l, phi, rho, grad_rho),
            options={"disp": display},
            callback=callback,
        )

        if extention == "alpha":
            best = {
                "theta": res.x[0],
                "delta": res.x[1],
                "fun": res.fun,
            }
        elif extention == "theta":
            best = {
                "theta": res.x[0],
                "fun": res.fun,
            }
        elif extention == "delta":
            best = {
                "delta": res.x[0],
                "fun": res.fun,
            }

        if best["fun"] < min_best["fun"]:
            mae = MAE(res.x, z[left:right], dens[left:right], l, phi, rho, grad_rho)
            rmse = RMSE(res.x, z[left:right], dens[left:right], l, phi, rho, grad_rho)

            best.update({"mae": mae, "rmse": rmse})

            min_best = best

    return z, dens, min_best


def profile_approx_modified(
    structure_file: str,
    dens_profile_file: str,
    rho_bulk: float,
    droplet_mols: list,
    interface_type: str,
    display: bool = True,
) -> tuple[np.array, np.array, float, float, float, dict]:
    """
    Modified approximation of a density profile, incorporating real pore height and volume fraction.

    Parameters:
    structure_file (str): Path to the structure file with a .gro extension.
    dens_profile_file (str): Path to the density profile file with an .xvg extension.
    rho_bulk (float): Bulk density of the droplet.
    droplet_mols (list): List of molecule identifiers in the droplet for processing.
    interface_type (str): Type of surface for approximation, e.g., 'roll'.
    display (bool, optional): If True, display optimization details. Default is True.

    Returns:
    tuple:
        - z (np.array): Normalized z-coordinates of the density profile.
        - dens (np.array): Normalized density profile values.
        - l (float): Adjusted pore length-to-height ratio.
        - real_phi (float): Real volume fraction calculated from density profile.
        - H (float): Real pore height based on droplet molecule positions.
        - best (dict): Best-fit parameters including 'theta' angle and 'offset'.
    """

    # Ensure the structure file has the correct extension
    assert structure_file.endswith("gro"), (
        "The structure file must have the extension .gro"
    )
    # Read the structure file
    structure = read_gro(structure_file)

    # Ensure the density profile file has the correct extension
    assert dens_profile_file.endswith("xvg"), (
        "The density profile file must have the extension .xvg"
    )
    # Read the density profile file
    z, dens = np.loadtxt(dens_profile_file, comments=["@", "#"], unpack=True)

    # Initialize z_min and z_max based on the structure box dimensions
    z_min, z_max = structure.box[2] / 2, structure.box[2] / 2
    # Process atoms of the droplet to find real pore height
    for i, atom_label in enumerate(structure.atoms):
        if atom_label.mol_name in droplet_mols:
            z_min = min(z_min, structure.atoms_xyz[i, 2])
            z_max = max(z_max, structure.atoms_xyz[i, 2])
    H = z_max - z_min
    print("real H:", H)

    # Calculate the adjusted pore length-to-height ratio
    l = structure.box[0] / H

    # Normalize the density profile
    dens /= rho_bulk
    z /= H

    # Integrate the density profile numerically to get the real volume fraction
    real_phi = np.trapz(dens, z)
    print("real phi:", real_phi)

    # Find the left non-zero bin
    left, right = 0, len(z)
    for i in range(len(z)):
        if z[i] > -0.5:
            left = i - 1
            break

    # Find the right non-zero bin
    for i in range(len(z) - 1, -1, -1):
        if z[i] < 0.5:
            right = i + 1
            break

    # Create a function based on the surface type
    assert interface_type.lower() in [
        "droplet",
        "doughnut",
        "worm",
        "roll",
        "perforation",
        "layer",
    ], (
        "There is no such type of surface. Could be droplet, doughnut, worm, roll, perforation or layer"
    )
    rho_theta = getattr(src.utils_py.auxil, f"rho_{interface_type.lower()}_theta")

    def L1(x, z, dens, l, phi):
        """
        Function to minimize for finding the best angle theta and offset.

        Parameters:
            x (list): List of parameters to optimize, [theta, offset].
            z (np.array): Normalized z-coordinates of the density profile.
            dens (np.array): Normalized density profile values.
            l (float): Ratio of pore length to pore height.
            phi (float): Volume fraction of the droplet.

        Returns:
            float: Sum of absolute differences between the calculated and the experimental densities.
        """
        return np.sum(np.abs(rho_theta(z, l, phi, x[0]) + x[1] - dens))

    # Minimization process
    x0 = [3 * np.pi / 4, 0]  # Initial guess for theta and offset
    res = minimize(
        L1,
        x0,
        (z[left:right], dens[left:right], l, real_phi),
        method="Nelder-Mead",
        bounds=((np.pi / 2, np.pi), (-0.5, 0.5)),
        options={"disp": display},
    )
    best = {"theta": res.x[0], "offset": res.x[1]}

    return z, dens, l, real_phi, H, best
