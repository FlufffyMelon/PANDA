import numpy as np
from numpy import pi, cos, sin, tan, sqrt, cbrt, abs

BIG_NUM = np.Inf

""" Droplet """


# Formulas in case of piº degree contact angle
def phi_min_droplet_pi(l):
    return 0


def phi_max_droplet_pi(l):
    return pi / (6 * l * l)


def S_droplet_pi(l, phi):
    phi_min = phi_min_droplet_pi(l)
    phi_max = phi_max_droplet_pi(l)

    if max(0, phi_min) < phi <= max(0, phi_max):
        return cbrt(36 * pi * phi * phi * l * l * l * l)

    return BIG_NUM


def rho_droplet_pi(z, l, phi):
    mask = abs(z) <= cbrt(6 * pi * pi * phi * l * l) / (2 * pi)
    z_cut = z * mask
    # Equation from article (incorrect)
    # rho = (cbrt(36 * pi * phi * phi * l * l * l * l) - pi * z_cut**2) / (l * l)
    rho = (cbrt(36 * pi * phi * phi * l * l * l * l) / 4 - pi * z_cut**2) / (l * l)

    return rho * mask


""" Doughnut """


# Formulas in case of piº degree contact angle
def phi_min_doughnut_pi(l):
    return pi / (6 * l * l)


def phi_max_doughnut_pi(l):
    return pi * (0.25 + (pi - 4) / (8 * l) + (10 - 3 * pi) / (24 * l * l))


def S_doughnut_pi(l, phi):
    phi_min = phi_min_doughnut_pi(l)
    phi_max = phi_max_doughnut_pi(l)

    if max(0, phi_min) < phi <= max(0, phi_max):
        a = -(pi**3) / 16 + 2 * pi / 3 + 2 * phi * l * l
        b = sqrt(3 * pi**3 * (3 * pi**3 - 32 * pi + 192 * phi * l * l))

        return a + b / 48

    return BIG_NUM


def rho_doughnut_pi(z, l, phi):
    # Equations from article (incorrect)
    # z_cut = z * (abs(z) < 0.5)
    # a = sqrt(pi * (3 * pi**3 + 192 * phi * l * l - 32 * pi) / 3) - pi * pi
    # b = 4 * pi * sqrt(1 - 4 * z_cut**2)

    # return (a + b) / (8 * l * l) * (abs(z) < 0.5)

    z_cut = z * (abs(z) < 0.5)
    a = (sqrt(pi * (3 * pi**3 + 192 * phi * l * l - 32 * pi) / 3) - pi * pi) / (4 * pi)
    b = sqrt(1 - 4 * z_cut**2)

    return 0.25 * pi * (a + b) ** 2 / (l * l) * (abs(z) < 0.5)


""" Worm """


# Formulas in case of piº degree contact angle
def phi_min_worm_pi(l):
    return 0


def phi_max_worm_pi(l):
    return pi / (4 * l)


def S_worm_pi(l, phi):
    phi_min = phi_min_worm_pi(l)
    phi_max = phi_max_worm_pi(l)

    if max(0, phi_min) < phi <= max(0, phi_max):
        return 2 * sqrt(pi * phi * l * l * l)

    return BIG_NUM


def rho_worm_pi(z, l, phi):
    mask = abs(z) <= sqrt(phi * l / pi)
    z_cut = z * mask
    rho = 2 * sqrt(phi * l / pi - z_cut**2) / l

    return rho * mask


""" Roll """


# Formulas in case of piº degree contact angle
def phi_min_roll_pi(l):
    return 0.25 * pi / l


def phi_max_roll_pi(l):
    return 1 + 0.25 * (pi - 4) / l


def S_roll_pi(l, phi):
    phi_min = phi_min_roll_pi(l)
    phi_max = phi_max_roll_pi(l)

    if max(0, phi_min) < phi <= max(0, phi_max):
        return 0.5 * pi * l + 2 * phi * l * l

    return BIG_NUM


def rho_roll_pi(z, l, phi):
    mask = abs(z) < 0.5
    z_cut = z * mask
    rho = 0.25 * (4 * phi * l + 4 * sqrt(1 - 4 * z_cut**2) - pi) / l

    return rho * mask


""" Perforation """


# Formulas in case of piº degree contact angle
def phi_min_perforation_pi(l):
    return 1 - 0.25 * pi + pi * pi / (8 * l) - pi / (6 * l * l)


def phi_max_perforation_pi(l):
    return 1 - pi * (5 / 12 - pi / 8) / (l * l)


def S_perforation_pi(l, phi):
    phi_min = phi_min_perforation_pi(l)
    phi_max = phi_max_perforation_pi(l)

    if max(0, phi_min) < phi <= max(0, phi_max):
        a = pi**3 / 16 - 2 * pi / 3 + 2 * phi * l * l
        b = sqrt(3 * pi**3 * (3 * pi**3 - 32 * pi + 192 * (1 - phi) * l * l))

        return a + b / 48

    return BIG_NUM


def rho_perforation_pi(z, l, phi):
    mask = abs(z) < 0.5
    z_cut = z * mask
    a = sqrt(3 * pi**3 + 192 * (1 - phi) * l * l - 32 * pi)
    b = sqrt(3 * pi) * (pi - 4 * sqrt(1 - 4 * z_cut**2))
    rho = 1 - (a + b) ** 2 / (192 * l * l)

    return rho * mask


""" Layer """


# Formulas in case of piº degree contact angle
def phi_min_layer_pi(l):
    return 0


def phi_max_layer_pi(l):
    return 1


def S_layer_pi(l, phi):
    phi_min = phi_min_layer_pi(l)
    phi_max = phi_max_layer_pi(l)

    if max(0, phi_min) < phi <= max(0, phi_max):
        return 2 * l * l

    return BIG_NUM


def rho_layer_pi(z, l, phi):
    return np.ones_like(z) * (abs(z) < 0.5 * phi)
