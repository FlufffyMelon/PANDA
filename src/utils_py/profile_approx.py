import numpy as np
import sys
sys.path.append('..') # Avoid error with importing of src
from src.utils_py.io.gro import read_gro
import src.utils_py.auxil
from scipy.optimize import minimize, Bounds


def profile_approx(
        dens_profile_file : str,
        rho_bulk : float,
        l: float,
        phi : float,
        H : float,
        surface_type : str,
        display : bool = True
    ) -> tuple[np.array, np.array, dict]:
    """
    Approximate a density profile using a classical algorithm from the article.

    Parameters:
    dens_profile_file (str): Path to the density profile file with an .xvg extension.
    rho_bulk (float): Bulk density of the droplet.
    l (float): Ratio of pore length to pore height.
    phi (float): Volume fraction of the droplet.
    H (float): Height of the pore.
    surface_type (str): Type of surface to use in approximation, e.g., 'roll'.
    display (bool, optional): If True, display optimization details. Default is True.

    Returns:
    tuple:
        - z (np.array): Normalized z-coordinates of the density profile.
        - dens (np.array): Normalized density profile values.
        - best (dict): Best-fit parameters, including 'theta' angle in radians.
    """

    assert dens_profile_file.endswith('xvg'), "The density profile file must have the extension .xvg"
    # Reading density profile file
    z, dens = np.loadtxt(dens_profile_file, comments=["@", "#"], unpack=True)

    dens /= rho_bulk
    z /= H

    left, right = 0, len(z)
    # Find left none zero bin
    for i in range(len(z)):
        if z[i] > -0.5:
            left = i-1
            break

    # Find right none zero bin
    for i in range(len(z)-1, -1, -1):
        if z[i] < 0.5:
            right = i+1
            break

    # Create a function from surface type
    assert surface_type.lower() in ['droplet', 'doughnut', 'worm', 'roll', 'perforation', 'layer'], \
    "There is no such type of surface. Could be droplet, doughnut, worm, roll, perforation or layer"
    rho_theta = getattr(src.utils_py.auxil, f'rho_{surface_type.lower()}_theta')
    def L1(x, z, dens, l, phi):
        return np.sum(np.abs(rho_theta(z, l, phi, x[0]) - dens))

    # Minimization process
    x0 = [3*np.pi/4]
    res = minimize(
        L1,
        x0,
        (z[left:right], dens[left:right], l, phi),
        method='Nelder-Mead',
        bounds=Bounds(np.pi / 2, np.pi),
        options={'disp': display}
    )
    best = {'theta': res.x[0]}

    return z, dens, best



def profile_approx_modified(
        structure_file : str,
        dens_profile_file : str,
        rho_bulk : float,
        surface_mols : list,
        surface_type : str,
        display : bool = True
    ) -> tuple[np.array, np.array, float, float, float, dict]:
    """
    Modified approximation of a density profile, incorporating real pore height and volume fraction.
    Some trick until formulas accounting for a non-zero wetting layer are obtained.

    Parameters:
    structure_file (str): Path to the structure file with a .gro extension.
    dens_profile_file (str): Path to the density profile file with an .xvg extension.
    rho_bulk (float): Bulk density of the droplet.
    surface_mols (list): List of molecule identifiers in the droplet for processing.
    surface_type (str): Type of surface for approximation, e.g., 'roll'.
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

    assert structure_file.endswith('gro'), "The structure file must have the extension .gro"
    # Reading structure file
    structure = read_gro(structure_file)

    assert dens_profile_file.endswith('xvg'), "The density profile file must have the extension .xvg"
    # Reading density profile file
    z, dens = np.loadtxt(dens_profile_file, comments=["@", "#"], unpack=True)

    z_min, z_max = structure.box[2] / 2, structure.box[2] / 2
    # Process atoms of droplet to find real pore height
    for i, atom_label in enumerate(structure.atoms):
        if atom_label.mol_name in surface_mols:
            z_min = min(z_min, structure.atoms_xyz[i, 2])
            z_max = max(z_max, structure.atoms_xyz[i, 2])
    H = z_max - z_min
    print('real H:', H)

    l = structure.box[0] / H

    dens /= rho_bulk
    z /= H

    # Integrate density profile numericaly to get real volume fraction
    real_phi = np.trapz(dens, z)
    print('real phi:', real_phi)

    left, right = 0, len(z)
    for i in range(len(z)):
        if z[i] > -0.5:
            left = i-1
            break

    for i in range(len(z)-1, -1, -1):
        if z[i] < 0.5:
            right = i+1
            break

    # Create a function from surface type
    assert surface_type.lower() in ['droplet', 'doughnut', 'worm', 'roll', 'perforation', 'layer'], \
    "There is no such type of surface. Could be droplet, doughnut, worm, roll, perforation or layer"
    rho_theta = getattr(src.utils_py.auxil, f'rho_{surface_type.lower()}_theta')
    def L1(x, z, dens, l, phi):
        return np.sum(np.abs(rho_theta(z, l, phi, x[0]) + x[1] - dens))

    # Minimization process
    x0 = [3*np.pi/4, 0]
    res = minimize(
        L1,
        x0,
        (z[left:right], dens[left:right], l, real_phi),
        method='Nelder-Mead',
        bounds=((np.pi / 2, np.pi), (-0.5, 0.5)),
        options={'disp': display}
    )
    best = {'theta': res.x[0], 'offset': res.x[1]}

    return z, dens, l, real_phi, H, best
