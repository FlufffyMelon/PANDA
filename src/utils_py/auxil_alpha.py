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

def grad_d_roll_alpha(l, phi, th, delta):
    # d_theta = 0.5 * (1 - 2 * delta) * (pi - 2 * th) * tan(th) / cos(th) ** 2
    d_theta = (1 - 2 * delta) * (0.5 * (pi - 2 * th) * tan(th) - 1) / cos(th) ** 2

    a = 2 * phi * l / (1 - 2 * delta) ** 2
    b = 0.5 * (pi - 2 * th) / cos(th) ** 2
    c = tan(th)
    d_delta = a - b + c

    return np.array([
        [d_theta],
        [d_delta]
        ])


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

def grad_y_roll_alpha(z, l, phi, th, delta):
    mask = abs(z) < 0.5 - delta
    z_cut = z * mask

    grad_d = grad_d_roll_alpha(l, phi, th, delta)

    sq = sqrt((1 - 2 * delta) ** 2 / cos(th) ** 2 - 4 * z_cut**2) * cos(th) ** 2
    grad_sqrt = np.array([
        [(1 - 2 * delta) ** 2 * tan(th)],
        [-2 * (1 - 2 * delta)]
    ]) / sq

    grad_rem = np.array([
        [(1 - 2 * delta) / cos(th) ** 2],
        [-2 * tan(th)]
    ])

    grad_y = 0.5 * (grad_d + grad_sqrt + grad_rem)

    return grad_y * mask

def rho_roll_alpha(z, l, phi, th, delta):
    return 2 * y_roll_alpha(z, l, phi, th, delta) / l

def grad_rho_roll_alpha(z, l, phi, th, delta):
    return 2 * grad_y_roll_alpha(z, l, phi, th, delta) / l
