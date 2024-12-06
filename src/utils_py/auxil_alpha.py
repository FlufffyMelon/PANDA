import numpy as np
from numpy import pi, cos, sin, tan, sqrt, cbrt, abs, arcsin

BIG_NUM = np.Inf


""" Roll """


# General formulas in case of arbitrary theta
def d_roll_alpha(l, phi, th, delta):
    a = phi * l / (1 - 2 * delta)
    b = 0.25 * (1 - 2 * delta) * (pi - 2 * th) / cos(th) ** 2
    c = 0.5 * (1 - 2 * delta) * tan(th)

    return a + b - c


def phi_min_roll_alpha(l, th, delta):
    pass


def phi_max_roll_alpha(l, th, delta):
    pass


def S_roll_alpha(l, phi, th, delta):
    pass


def y_roll_alpha(z, l, phi, th, delta):
    d = d_roll_alpha(l, phi, th, delta)
    mask = abs(z) < 0.5 - delta
    z_cut = z * mask
    y = 0.5 * (
        d
        + sqrt((1 - 2 * delta) ** 2 / cos(th) ** 2 - 4 * z_cut**2)
        + (1 - 2 * delta) * tan(th)
    )

    return y * mask


def rho_roll_alpha(z, l, phi, th, delta):
    return 2 * y_roll_alpha(z, l, phi, th, delta) / l
