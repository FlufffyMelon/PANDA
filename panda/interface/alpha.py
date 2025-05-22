import numpy as np
from numpy import pi, cos, sin, tan, sqrt, cbrt, abs, arcsin

BIG_NUM = np.Inf


""" Droplet """


# General formulas in case of arbitrary theta and delta
def r_droplet_alpha(l, phi, th, delta):
    return cbrt(3 * phi * l**2 / (4 * pi * (2 + cos(th)) * sin(0.5 * th) ** 4))


def grad_r_droplet_alpha(l, phi, th, delta):
    return NotImplementedError


def phi_min_droplet_alpha(l, th, delta):
    return 0


def phi_max_droplet_alpha(l, th, delta):
    phi_H = pi * (1 - 2 * delta) ** 3 * (2 + cos(th)) / (3 * l**2 * (1 - cos(th)))

    if th >= 0.5 * pi:
        phi_l = pi * l * (2 + cos(th)) * sin(0.5 * th) ** 4 / 6
    else:
        phi_l = pi * l * (2 + cos(th)) * sin(0.5 * th) ** 4 / sin(th) ** 3 / 6

    return min(phi_l, phi_H)


def S_droplet_alpha(l, phi, th, delta):
    return NotADirectoryError


def y_droplet_alpha(z, l, phi, th, delta, center=False):
    r = r_droplet_alpha(l, phi, th, delta)
    z_dc = 3 * cos(0.5 * th) ** 4 / (2 + cos(th)) * r if center else 0
    mask = np.logical_and(cos(th) < (z + z_dc) / r, (z + z_dc) / r < 1)
    z_cut = z * mask
    y = sqrt(r**2 - (z_cut + z_dc) ** 2)

    return y * mask


def grad_y_droplet_alpha(z, l, phi, th, delta, center=False):
    return NotImplementedError


def rho_droplet_alpha(z, l, phi, th, delta, center=True):
    return pi * y_droplet_alpha(z, l, phi, th, delta, center=center) ** 2 / l**2


def grad_rho_droplet_alpha(z, l, phi, th, delta, center=True):
    # return pi * grad_y_droplet_alpha(z, l, phi, th, delta, center=center) ** 2 / l**2
    return NotImplementedError


""" Doughnut """

# General formulas in case of arbitrary theta and delta
"""
d_{dou} = \frac{1}{4} \sqrt{\frac{(\pi-2\theta)^2}{\cos ^4 \theta }-\frac{12 - 4(\pi -2 \theta) \tan (\theta)}{\cos ^2 \theta} + \frac{64 \phi l^2}{\pi} + \frac{4}{3}} -\frac{\tan \theta}{2}+\frac{1}{4} \frac{\pi -2 \theta} {\cos ^2 \theta },
"""


def d_doughnut_alpha(l, phi, th, delta):
    a = (1 - 2 * delta) ** 2 * (pi - 2 * th) ** 2 / cos(th) ** 4
    b = (12 - 4 * (pi - 2 * th) * tan(th)) * (1 - 2 * delta) ** 2 / cos(th) ** 2
    c = 64 * phi * l * l / pi / (1 - 2 * delta)
    d = 4 * (1 - 2 * delta) ** 2 / 3
    e = 0.5 * (1 - 2 * delta) * tan(th)
    f = 0.25 * (1 - 2 * delta) * (pi - 2 * th) / cos(th) ** 2

    return 0.25 * sqrt(a - b + c + d) - e + f


def grad_d_doughnut_alpha(l, phi, th, delta):
    return NotImplementedError


def phi_min_doughnut_alpha(l, th, delta):
    a = pi * (1 - 2 * delta) ** 3 / l**2
    b = tan(th) ** 2 / 4 + 1 / 6
    c = (pi - 2 * th) * tan(th) / cos(th) ** 2 / 8

    return a * (b - c)


"""
\phi_{max}=\frac{\pi\left(1-2\delta\right)}{16l^{2}}\left[(1-2\delta)^{2}\tan^{2}(\theta)-\frac{2(1-2\delta)^{2}(\pi-2\theta)}{\cos(\theta)^{3}}-\frac{\left(1-2\delta\right)^{2}}{3}+4l^{2}-4(1-2\delta)l\tan(\theta)+ \frac{(1-2\delta)(7(1-2\delta)-2(\pi-2\theta)l)}{\cos(\theta)^{2}}-\frac{4(1-2\delta)((1-2\delta)\tan(\theta)-2l)}{\cos(\theta)}\right]
"""


def phi_max_doughnut_alpha(l, th, delta):
    a = pi * (1 - 2 * delta) / 16 / l**2
    b = (1 - 2 * delta) ** 2 * tan(th) ** 2
    c = 2 * (1 - 2 * delta) ** 2 * (pi - 2 * th) / cos(th) ** 3
    d = (1 - 2 * delta) ** 2 / 3 - 4 * l**2
    e = 4 * (1 - 2 * delta) * l * tan(th)
    f = (1 - 2 * delta) * (7 * (1 - 2 * delta) - 2 * (pi - 2 * th) * l) / cos(th) ** 2
    g = 4 * (1 - 2 * delta) * ((1 - 2 * delta) * tan(th) - 2 * l) / cos(th)

    return a * (b - c - d - e + f - g)


def S_doughnut_alpha(l, phi, th, delta):
    return NotImplementedError


def y_doughnut_alpha(z, l, phi, th, delta):
    d = d_doughnut_alpha(l, phi, th, delta)
    mask = abs(z) < 0.5 - delta
    z_cut = z * mask
    y = 0.5 * (
        d
        + sqrt((1 - 2 * delta) ** 2 / cos(th) ** 2 - 4 * z_cut**2)
        + (1 - 2 * delta) * tan(th)
    )

    return y * mask


def grad_y_doughnut_alpha(z, l, phi, th, delta):
    return NotImplementedError


def rho_doughnut_alpha(z, l, phi, th, delta):
    return pi * y_doughnut_alpha(z, l, phi, th, delta) ** 2 / l**2


def grad_rho_doughnut_alpha(z, l, phi, th, delta):
    # return pi * grad_y_doughnut_alpha(z, l, phi, th, delta) ** 2 / l**2
    return NotImplementedError


""" Worm """


# General formulas in case of arbitrary theta and delta
def r_worm_alpha(l, phi, th, delta):
    return sqrt(2 * phi * l / (2 * th - sin(2 * th)))


def grad_r_worm_alpha(l, phi, th, delta):
    return NotImplementedError


def phi_min_worm_alpha(l, th, delta):
    return 0


def phi_max_worm_alpha(l, th, delta):
    phi_H = (1 - 2 * delta) ** 2 * (2 * th - sin(2 * th)) / (2 * l * (1 - cos(th)) ** 2)

    if th > 0.5 * pi:
        phi_l = l * (2 * th - sin(2 * th)) / 8
    else:
        phi_l = l * (2 * th - sin(2 * th)) / sin(th) ** 2 / 8

    return min(phi_l, phi_H)


def S_worm_alpha(l, phi, th, delta):
    return NotImplementedError


def y_worm_alpha(z, l, phi, th, delta, center=True):
    r = r_worm_alpha(l, phi, th, delta)
    z_wc = 4 * sin(th) ** 3 * r / (3 * (2 * th - sin(2 * th))) if center else 0
    mask = np.logical_and(cos(th) < (z + z_wc) / r, (z + z_wc) / r < 1)
    z_cut = z * mask
    y = sqrt(r * r - (z_cut + z_wc) ** 2)

    return y * mask


def grad_y_worm_alpha(z, l, phi, th, delta, center=True):
    return NotImplementedError


def rho_worm_alpha(z, l, phi, th, delta, center=True):
    return 2 * y_worm_alpha(z, l, phi, th, delta, center=center) / l


def grad_rho_worm_theta(z, l, phi, th, delta, center=True):
    # return 2 * grad_y_worm_alpha(z, l, phi, th, delta, center=center) / l
    return NotImplementedError


""" Roll """


# General formulas in case of arbitrary theta
def d_roll_alpha(l, phi, th, delta):
    a = phi * l / (1 - 2 * delta)
    b = 0.25 * (1 - 2 * delta) * (pi - 2 * th) / cos(th) ** 2
    c = 0.5 * (1 - 2 * delta) * tan(th)

    return a + b - c


def grad_d_roll_alpha(l, phi, th, delta):
    d_theta = (1 - 2 * delta) * (0.5 * (pi - 2 * th) * tan(th) - 1) / cos(th) ** 2

    a = 2 * phi * l / (1 - 2 * delta) ** 2
    b = 0.5 * (pi - 2 * th) / cos(th) ** 2
    c = tan(th)
    d_delta = a - b + c

    return np.array([[d_theta], [d_delta]])


def phi_min_roll_alpha(l, th, delta):
    return (
        0.5 * (1 - 2 * delta) ** 2 * (tan(th) - 0.5 * (pi - 2 * th) / cos(th) ** 2) / l
    )


def phi_max_roll_alpha(l, th, delta):
    a = 1 / cos(th) - 0.25 * (pi - 2 * th) / cos(th) ** 2 - 0.5 * tan(th)

    return (1 - 2 * delta) + (1 - 2 * delta) ** 2 * a / l


def S_roll_alpha(l, phi, th, delta):
    return NotImplementedError


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
    grad_sqrt = (
        np.array([[(1 - 2 * delta) ** 2 * tan(th)], [-2 * (1 - 2 * delta)]]) / sq
    )

    grad_rem = np.array([[(1 - 2 * delta) / cos(th) ** 2], [-2 * tan(th)]])

    grad_y = 0.5 * (grad_d + grad_sqrt + grad_rem)

    return grad_y * mask


def rho_roll_alpha(z, l, phi, th, delta):
    return 2 * y_roll_alpha(z, l, phi, th, delta) / l


def grad_rho_roll_alpha(z, l, phi, th, delta):
    return 2 * grad_y_roll_alpha(z, l, phi, th, delta) / l


""" Perforation """


# General formulas in case of arbitrary theta and delta
def d_perforation_alpha(l, phi, th, delta):
    a = (1 - 2 * delta) ** 2 * (pi - 2 * th) ** 2 / cos(th) ** 4
    b = (12 - 4 * (pi - 2 * th) * tan(th)) * (1 - 2 * delta) ** 2 / cos(th) ** 2
    c = 64 * (1 - phi / (1 - 2 * delta)) * l**2 / pi
    d = 4 * (1 - 2 * delta) ** 2 / 3
    e = 0.5 * (1 - 2 * delta) * tan(th)
    f = 0.25 * (1 - 2 * delta) * (pi - 2 * th) / cos(th) ** 2

    return 0.25 * sqrt(a - b + c + d) + e - f


def phi_min_perforation_alpha(l, th, delta):
    a = (1 - 2 * delta) * (1 - pi / 4)
    b = pi * (1 - 2 * delta) ** 2 / 8 / l**2
    c = 2 * (1 - 2 * delta) / 3 + 2 * l * tan(th)
    d = (
        2 * (1 - 2 * delta)
        + (pi - 2 * th) * l
        - (1 - 2 * delta) * (pi - 2 * th) * tan(th)
    ) / cos(th) ** 2

    return a + b * (c - d)


def phi_max_perforation_alpha(l, th, delta):
    a = pi * (1 - 2 * delta) ** 3 / (16 * l**2 * cos(th) ** 3)
    b = 7 * cos(th) - cos(3 * th) / 3 - 2 * (pi - 2 * th + sin(2 * th))

    return 1 - 2 * delta - a * b


def S_perforation_alpha(l, phi, th, delta):
    return NotImplementedError


def y_perforation_alpha(z, l, phi, th, delta):
    d = d_perforation_alpha(l, phi, th, delta)
    mask = abs(z) < 0.5 - delta
    z_cut = z * mask
    y = 0.5 * (
        d
        - sqrt((1 - 2 * delta) ** 2 / cos(th) ** 2 - 4 * z_cut**2)
        - (1 - 2 * delta) * tan(th)
    )

    return y * mask


def grad_y_perforation_alpha(z, l, phi, th, delta):
    return NotImplementedError


def rho_perforation_alpha(z, l, phi, th, delta):
    mask = abs(z) < 0.5 - delta
    rho = 1 - pi * y_perforation_alpha(z, l, phi, th, delta) ** 2 / l**2

    return rho * mask


def grad_rho_perforation_theta(z, l, phi, th, delta):
    # mask = abs(z) < 0.5 - delta
    # rho = 1 - pi * grad_y_perforation_alpha(z, l, phi, th, delta) ** 2 / (l * l)

    # return rho * mask
    return NotImplementedError


""" Layer """


# General formulas in case of arbitrary theta and delta
def phi_min_layer_alpha(l, th, delta):
    return 0


def phi_max_layer_alpha(l, th, delta):
    return 1 - 2 * delta


def S_layer_alpha(l, phi, th, delta):
    return NotImplementedError


def y_layer_alpha(z, l, phi, th, delta):
    mask = abs(z) < 0.5 * phi

    return np.ones_like(z) * 0.5 * l * mask


def rho_layer_alpha(z, l, phi, th, delta):
    mask = abs(z) < 0.5 * phi

    return np.ones_like(z) * mask
