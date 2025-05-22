import numpy as np
from numpy import pi, cos, sin, tan, sqrt, cbrt, abs

BIG_NUM = np.inf

""" Droplet """


# General formulas in case of arbitrary theta
def r_droplet_theta(l, phi, th):
    return cbrt(3 * phi * l ** 2 / (4 * pi * (2 + cos(th)) * sin(0.5 * th) ** 4))


def phi_min_droplet_theta(l, th):
    return 0


def phi_max_droplet_theta(l, th):
    # phi_H = pi * (2 + cos(th)) / (3 * l * l * (1 - cos(th)))

    # if th >= 0.5 * pi:
    #     phi_l = pi * l * (2 + cos(th)) * sin(0.5 * th)**4 / 6
    # else:
    #     phi_l = pi * l * (2 + cos(th)) * sin(0.5 * th)**4 / sin(th)**3 / 6

    return pi * (2 + cos(th)) / (3 * l ** 2 * (1 - cos(th)))
    # return min(phi_H, phi_l)


def S_droplet_theta(l, phi, th):
    phi_min = phi_min_droplet_theta(l, th)
    phi_max = phi_max_droplet_theta(l, th)

    if max(0, phi_min) < phi <= max(0, phi_max):
        r = r_droplet_theta(l, phi, th)

        return pi * r * r * (2 + sin(th) - 2 * cos(th))

    return BIG_NUM


def y_droplet_theta(z, l, phi, th, center=False):
    r = r_droplet_theta(l, phi, th)
    z_dc = 3 * cos(0.5 * th) ** 4 / (2 + cos(th)) * r if center else 0
    mask = np.logical_and(cos(th) < (z + z_dc) / r, (z + z_dc) / r < 1)
    z_cut = z * mask
    y = sqrt(r ** 2 - (z_cut + z_dc) ** 2)

    return y * mask


def rho_droplet_theta(z, l, phi, th, center=True):
    return pi * y_droplet_theta(z, l, phi, th, center=center) ** 2 / (l * l)


""" Doughnut """


# General formulas in case of arbitrary theta
def d_doughnut_theta(l, phi, th):
    a = (pi - 2 * th) ** 2 / cos(th) ** 4
    b = (4 + 4 * (pi - 2 * th) * tan(th)) / cos(th) ** 2
    c = 8 * (cos(2 * th) - 5) / (3 * cos(th) ** 2)
    d = 64 * phi * l * l / pi - 4
    e = 0.5 * tan(th)
    f = 0.25 * (pi - 2 * th) / cos(th) ** 2

    return 0.25 * sqrt(a + b + c + d) - e + f


def phi_min_doughnut_theta(l, th):
    return (
        -pi
        * (-5 + cos(2 * th) + 3 * (pi - 2 * th) * tan(th))
        / (24 * l * l * cos(th) ** 2)
    )


def phi_max_doughnut_theta(l, th):
    a = 9 * cos(th) - cos(3 * th)
    b = l - pi + 2 * th + 2 * cos(th) + l * cos(2 * th) - sin(2 * th)

    return pi * (a + 6 * (1 + l * cos(th)) * b) / (48 * l * l * cos(th) ** 3)


def S_doughnut_theta(l, phi, th):
    phi_min = phi_min_doughnut_theta(l, th)
    phi_max = phi_max_doughnut_theta(l, th)

    if max(0, phi_min) < phi <= max(0, phi_max):
        d = d_doughnut_theta(l, phi, th)

        return (
            pi * d * d / 2
            - pi / cos(th)
            + pi * (pi - 2 * th) * (d + tan(th)) / (2 * cos(th))
        )

    return BIG_NUM


def y_doughnut_theta(z, l, phi, th):
    d = d_doughnut_theta(l, phi, th)
    mask = abs(z) < 0.5
    z_cut = z * mask
    y = 0.5 * (d + sqrt(tan(th) ** 2 + 1 - 4 * z_cut**2) + tan(th))

    return y * mask


def rho_doughnut_theta(z, l, phi, th):
    return pi * y_doughnut_theta(z, l, phi, th) ** 2 / (l * l)


""" Worm """


# General formulas in case of arbitrary theta
def r_worm_theta(l, phi, th):
    return sqrt(2 * phi * l / (2 * th - sin(2 * th)))


def phi_min_worm_theta(l, th):
    return 0


def phi_max_worm_theta(l, th):
    # phi_H = (2 * th - sin(2 * th)) / (2 * l * (1 - cos(th))**2)

    # if th > 0.5 * pi:
    #     phi_l = l * (2 * th - sin(2 * th)) / 8
    # else:
    #     phi_l = l * (2 * th - sin(2 * th)) / sin(th)**2 / 8

    return (2 * th - sin(2 * th)) / (2 * l * (1 - cos(th)) ** 2)
    # return min(phi_H, phi_l)


def S_worm_theta(l, phi, th):
    phi_min = phi_min_worm_theta(l, th)
    phi_max = phi_max_worm_theta(l, th)

    if max(0, phi_min) < phi <= max(0, phi_max):
        r = r_worm_theta(l, phi, th)

        return 2 * l * r * (th + sin(th))

    return BIG_NUM


def y_worm_theta(z, l, phi, th, center=True):
    r = r_worm_theta(l, phi, th)
    z_wc = 4 * sin(th) ** 3 * r / (3 * (2 * th - sin(2 * th))) if center else 0
    mask = np.logical_and(cos(th) < (z + z_wc) / r, (z + z_wc) / r < 1)
    z_cut = z * mask
    y = sqrt(r * r - (z_cut + z_wc) ** 2)

    return y * mask


def rho_worm_theta(z, l, phi, th, center=True):
    return 2 * y_worm_theta(z, l, phi, th, center=center) / l


""" Roll """


# General formulas in case of arbitrary theta
def d_roll_theta(l, phi, th):
    return 0.25 * (
        4 * phi * l + pi - 2 * th - 2 * tan(th) + (pi - 2 * th) * tan(th) ** 2
    )


def phi_min_roll_theta(l, th):
    return -0.25 * (pi - 2 * th) / (l * cos(th) ** 2) + 0.5 * tan(th) / l


def phi_max_roll_theta(l, th):
    return (
        1
        + 1 / (l * cos(th))
        - 0.25 * (pi - 2 * th) / (l * cos(th) ** 2)
        - 0.5 * tan(th) / l
    )


def S_roll_theta(l, phi, th):
    phi_min = phi_min_roll_theta(l, th)
    phi_max = phi_max_roll_theta(l, th)

    if max(0, phi_min) < phi <= max(0, phi_max):
        d = d_roll_theta(l, phi, th)

        return 2 * l * d + l * (pi - 2 * th) / cos(th)

    return BIG_NUM


def y_roll_theta(z, l, phi, th):
    d = d_roll_theta(l, phi, th)
    mask = abs(z) < 0.5
    z_cut = z * mask
    y = 0.5 * (d + sqrt(tan(th) ** 2 + 1 - 4 * z_cut**2) + tan(th))

    return y * mask


def rho_roll_theta(z, l, phi, th):
    return 2 * y_roll_theta(z, l, phi, th) / l


""" Perforation """


# General formulas in case of arbitrary theta
def d_perforation_theta(l, phi, th):
    a = (pi - 2 * th) ** 2 / cos(th) ** 4
    b = (12 - 4 * (pi - 2 * th) * tan(th)) / cos(th) ** 2
    c = 64 * (1 - phi) * l * l / pi + 4 / 3
    d = 0.5 * tan(th) - 0.25 * (pi - 2 * th) / cos(th) ** 2

    return 0.25 * sqrt(a - b + c) + d


def phi_min_perforation_theta(l, th):
    a = 1 - 0.25 * pi + pi / (12 * l * l)
    b = 2 + (pi - 2 * th) * (l - tan(th))

    return a + 0.25 * pi * tan(th) / l - pi * b / (8 * l * l * cos(th) ** 2)


def phi_max_perforation_theta(l, th):
    c = 1 + pi / (48 * l * l)
    e = cos(3 * th) - 29 * cos(th) + 8 * (pi - 2 * th + sin(2 * th))

    return c + pi * e / (64 * l * l * cos(th) ** 3)


def S_perforation_theta(l, phi, th):
    phi_min = phi_min_perforation_theta(l, th)
    phi_max = phi_max_perforation_theta(l, th)

    if max(0, phi_min) < phi <= max(0, phi_max):
        d = d_perforation_theta(l, phi, th)

        return (
            2 * l * l
            - 0.5 * pi * d * d
            + pi / cos(th)
            + 0.5 * pi * (pi - 2 * th) * (d - tan(th)) / cos(th)
        )

    return BIG_NUM


def y_perforation_theta(z, l, phi, th):
    d = d_perforation_theta(l, phi, th)
    mask = abs(z) < 0.5
    z_cut = z * mask
    y = 0.5 * (d - sqrt(tan(th) ** 2 + 1 - 4 * z_cut**2) - tan(th))

    return y * mask


def rho_perforation_theta(z, l, phi, th):
    rho = 1 - pi * y_perforation_theta(z, l, phi, th) ** 2 / (l * l)

    return rho * (abs(z) < 0.5)
    # return 1 - pi * y_perforation_theta(z, l, phi, th)**2 / (l * l)


""" Layer """


# General formulas in case of arbitrary theta
def phi_min_layer_theta(l, th):
    return 0


def phi_max_layer_theta(l, th):
    return 1


def S_layer_theta(l, phi, th):
    phi_min = phi_min_layer_theta(l, th)
    phi_max = phi_max_layer_theta(l, th)

    if max(0, phi_min) < phi <= max(0, phi_max):
        return 2 * l * l

    return BIG_NUM


def y_layer_theta(z, l, phi, th):
    mask = abs(z) < 0.5 * phi

    return np.ones_like(z) * l * l * mask


def rho_layer_theta(z, l, phi, th):
    return np.ones_like(z) * (abs(z) < 0.5 * phi)
