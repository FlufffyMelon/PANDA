import numpy as np
from numpy import pi, cos, sin, tan, sqrt, cbrt, abs

BIG_NUM = np.inf

""" Droplet """


# Formulas with non zero wetting layer and 180º contact angle
def r_droplet_delta(l, phi, delta):
    return cbrt(0.25 * 3 * phi * l**2 / pi)


def phi_min_droplet_delta(l, delta):
    return 0


def phi_max_droplet_delta(l, delta):
    return pi * (1 - 2 * delta) ** 3 / (6 * l**2)


def S_droplet_delta(l, phi, delta):
    phi_min = phi_min_droplet_delta(l, delta)
    phi_max = phi_max_droplet_delta(l, delta)

    if max(0, phi_min) < phi <= max(0, phi_max):
        r = r_droplet_delta(l, phi, delta)

        return 4 * pi * r**2

    return BIG_NUM


def y_droplet_delta(z, l, phi, delta, center=True):
    r = r_droplet_delta(l, phi, delta)

    mask = np.abs(z) < r
    z_cut = z * mask
    y = sqrt(r * r - z_cut**2)

    return y * mask


def rho_droplet_delta(z, l, phi, delta, center=True):
    return pi * y_droplet_delta(z, l, phi, delta, center=center) ** 2 / l**2


""" Doughnut """


# Formulas with non zero wetting layer and 180º contact angle
def d_doughnut_delta(l, phi, delta):
    a = -0.25 * pi * (1 - 2 * delta)
    b = (0.25 * pi**2 - 8 / 3) * (0.5 - delta) ** 2 + 4 * phi * l**2 / (
        pi * (1 - 2 * delta)
    )

    return a + sqrt(b)


def phi_min_doughnut_delta(l, delta):
    return pi * (1 - 2 * delta) ** 3 / (6 * l**2)


def phi_max_doughnut_delta(l, delta):
    a = 0.25 * pi * (1 - 2 * delta)
    b = 0.5 * pi * (1 - 0.25 * pi) * (1 - 2 * delta) ** 2 / l
    c = 0.25 * pi * (5 / 3 - 0.5 * pi) * (1 - 2 * delta) ** 3 / l**2

    return a - b + c


def S_doughnut_delta(l, phi, delta):
    phi_min = phi_min_doughnut_delta(l, delta)
    phi_max = phi_max_doughnut_delta(l, delta)

    if max(0, phi_min) < phi <= max(0, phi_max):
        d = d_doughnut_delta(l, phi, delta)

        return (
            0.5 * pi * d**2
            + 0.5 * pi**2 * d * (1 - 2 * delta)
            + pi * (1 - 2 * delta) ** 2
        )

    return BIG_NUM


def y_doughnut_delta(z, l, phi, delta):
    d = d_doughnut_delta(l, phi, delta)
    mask = abs(z) < 0.5 - delta
    z_cut = z * mask
    y = 0.5 * (sqrt((1 - 2 * delta) ** 2 - 4 * z_cut**2) + d)

    return y * mask


def rho_doughnut_delta(z, l, phi, delta):
    return pi * y_doughnut_delta(z, l, phi, delta) ** 2 / l**2


""" Worm """


# Formulas with non zero wetting layer and 180º contact angle
def r_worm_delta(l, phi, delta):
    return sqrt(phi * l / pi)


def phi_min_worm_delta(l, delta):
    return 0


def phi_max_worm_delta(l, delta):
    return 0.25 * pi * (1 - 2 * delta) ** 2 / l


def S_worm_delta(l, phi, delta):
    phi_min = phi_min_worm_delta(l, delta)
    phi_max = phi_max_worm_delta(l, delta)

    if max(0, phi_min) < phi <= max(0, phi_max):
        r = r_worm_delta(l, phi, delta)

        return 2 * pi * l * r

    return BIG_NUM


def y_worm_delta(z, l, phi, delta, center=True):
    r = r_worm_delta(l, phi, delta)

    mask = np.abs(z) < r
    z_cut = z * mask
    y = sqrt(r**2 - z_cut**2)

    return y * mask


def rho_worm_delta(z, l, phi, delta, center=True):
    return 2 * y_worm_delta(z, l, phi, delta, center=center) / l


""" Roll """


# Formulas with non zero wetting layer and 180º contact angle
def d_roll_delta(l, phi, delta):
    return phi * l / (1 - 2 * delta) - 0.25 * pi * (1 - 2 * delta)


def phi_min_roll_delta(l, delta):
    return 0.25 * pi * (1 - 2 * delta) ** 2 / l


def phi_max_roll_delta(l, delta):
    return (1 - 2 * delta) - (1 - 0.25 * pi) * (1 - 2 * delta) ** 2 / l


def S_roll_delta(l, phi, delta):
    phi_min = phi_min_roll_delta(l, delta)
    phi_max = phi_max_roll_delta(l, delta)

    if max(0, phi_min) < phi <= max(0, phi_max):
        d = d_roll_delta(l, phi, delta)

        return 2 * l * d + pi * l * (1 - 2 * delta)

    return BIG_NUM


def y_roll_delta(z, l, phi, delta):
    d = d_roll_delta(l, phi, delta)
    mask = abs(z) < 0.5 - delta
    z_cut = z * mask
    y = 0.5 * (sqrt((1 - 2 * delta) ** 2 - 4 * z_cut**2) + d)

    return y * mask


def rho_roll_delta(z, l, phi, delta):
    return 2 * y_roll_delta(z, l, phi, delta) / l


""" Perforation """


# Formulas with non zero wetting layer and 180º contact angle
def d_perforation_delta(l, phi, delta):
    a = 0.25 * pi * (1 - 2 * delta)
    b = (0.25 * pi**2 - 8 / 3) * (0.5 - delta) ** 2 + 4 * (
        1 - phi / (1 - 2 * delta)
    ) * l**2 / pi

    return a + sqrt(b)


def phi_min_perforation_delta(l, delta):
    a = (1 - 0.25 * pi) * (1 - 2 * delta)
    b = pi**2 * (1 - 2 * delta) ** 2 / (8 * l)
    c = pi * (1 - 2 * delta) ** 3 / (6 * l**2)

    return a + b - c


def phi_max_perforation_delta(l, delta):
    return 1 - 2 * delta - 0.25 * pi * (5 / 3 - pi / 2) * (1 - 2 * delta) ** 3 / l**2


def S_perforation_delta(l, phi, delta):
    phi_min = phi_min_perforation_delta(l, delta)
    phi_max = phi_max_perforation_delta(l, delta)

    if max(0, phi_min) < phi <= max(0, phi_max):
        d = d_perforation_delta(l, phi, delta)

        return (
            2 * l**2
            - 0.5 * pi * d**2
            - pi * (1 - 2 * delta) ** 2
            + 0.5 * pi**2 * d * (1 - 2 * delta)
        )

    return BIG_NUM


def y_perforation_delta(z, l, phi, delta):
    d = d_perforation_delta(l, phi, delta)

    mask = abs(z) < 0.5 - delta
    z_cut = z * mask
    y = 0.5 * (-sqrt((1 - 2 * delta) ** 2 - 4 * z_cut**2) + d)

    return y * mask


def rho_perforation_delta(z, l, phi, delta):
    rho = 1 - pi * y_perforation_delta(z, l, phi, delta) ** 2 / (l * l)

    return rho * (abs(z) < 0.5 - delta)
    # return 1 - pi * y_perforation_delta(z, l, phi, delta)**2 / (l * l)


""" Layer """


# Formulas with non zero wetting layer and 180º contact angle
def phi_min_layer_delta(l, delta):
    return 0


def phi_max_layer_delta(l, delta):
    return 1 - 2 * delta


def S_layer_delta(l, phi, delta):
    phi_min = phi_min_layer_delta(l, delta)
    phi_max = phi_max_layer_delta(l, delta)

    if max(0, phi_min) < phi <= max(0, phi_max):
        return 2 * l**2

    return BIG_NUM


def y_layer_delta(z, l, phi, delta):
    mask = abs(z) < 0.5 * phi

    return np.ones_like(z) * l * l * mask


def rho_layer_delta(z, l, phi, delta):
    return np.ones_like(z) * (abs(z) < 0.5 * phi)
