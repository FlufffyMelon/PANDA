import numpy as np
from numpy import pi, cos, sin, tan, sqrt, cbrt, abs

BIG_NUM = np.Inf

''' Droplet '''
# General formulas
def r_droplet(l, phi, th):
    return cbrt(3 * phi * l * l / (4 * pi * (2 + cos(th)) * sin(0.5 * th)**4))

def phi_min_droplet(l, th):
    return 0

def phi_max_droplet(l, th):
    return pi * (2 + cos(th)) / (3 * l * l * (1 - cos(th)))

def S_droplet(l, phi, th):
    phi_min = phi_min_droplet(l, th)
    phi_max = phi_max_droplet(l, th)

    if max(0, phi_min) < phi <= max(0, phi_max):
        r = r_droplet(l, phi, th)

        return pi * r * r * (2 + sin(th) - 2 * cos(th))

    return BIG_NUM

def y_droplet(z, l, phi, th, center=True):
    r = r_droplet(l, phi, th)
    z_dc = 3 * cos(0.5 * th)**4 / (2 + cos(th)) * r if center else 0
    mask = np.logical_and(cos(th) < (z + z_dc) / r, (z + z_dc) / r < 1)
    z_cut = z * mask
    y = sqrt(r * r - (z_cut + z_dc)**2)

    return y * mask

def rho_droplet(z, l, phi, th, center=True):
    return pi * y_droplet(z, l, phi, th, center=center)**2 / (l * l)

# Formulas in case of 180º degree contact angle
def phi_min_droplet_180(l):
    return 0

def phi_max_droplet_180(l):
    return pi / (6 * l * l)

def S_droplet_180(l, phi):
    phi_min = phi_min_droplet_180(l)
    phi_max = phi_max_droplet_180(l)

    if max(0, phi_min) < phi <= max(0, phi_max):
        return cbrt(36 * pi * phi * phi * l * l * l * l)

    return BIG_NUM

def rho_droplet_180(z, l, phi):
    mask = abs(z) <= cbrt(6 * pi * pi * phi * l * l) / (2 * pi)
    z_cut = z * mask
    # Equation from article (incorrect)
    # rho = (cbrt(36 * pi * phi * phi * l * l * l * l) - pi * z_cut**2) / (l * l)
    rho = (cbrt(36 * pi * phi * phi * l * l * l * l) / 4 - pi * z_cut**2) / (l * l)

    return rho * mask



''' Doughnut '''
# General formulas
def d_doughnut(l, phi, th):
    a = (pi - 2 * th)**2 / cos(th)**4
    b = (4 + 4 * (pi - 2 * th) * tan(th)) / cos(th)**2
    c = 8 * (cos(2 * th) - 5) / (3 * cos(th)**2)
    d = 64 * phi * l * l / pi - 4
    e = 0.5 * tan(th)
    f = 0.25 * (np.pi - 2 * th) / cos(th)**2

    return 0.25 * sqrt(a + b + c + d) - e + f

def phi_min_doughnut(l, th):
    return - pi * (-5 + cos(2 * th) + 3 * (pi - 2 * th) * tan(th)) / (24 * l * l * cos(th)**2)

def phi_max_doughnut(l, th):
    a = 9 * cos(th) - cos(3 * th)
    b = l - pi + 2 * th + 2 * cos(th) + l * cos(2 * th) - sin(2 * th)

    return pi * (a + 6 * (1 + l * cos(th)) * b) / (48 * l * l * np.cos(th)**3)

def S_doughnut(l, phi, th):
    phi_min = phi_min_doughnut(l, th)
    phi_max = phi_max_doughnut(l, th)

    if max(0, phi_min) < phi <= max(0, phi_max):
        d = d_doughnut(l, phi, th)

        return pi * d * d / 2 - pi / cos(th) + pi * (pi - 2 * th) * (d + tan(th)) / (2 * cos(th))

    return BIG_NUM

def y_doughnut(z, l, phi, th):
    d = d_doughnut(l, phi, th)
    mask = abs(z) < 0.5
    z_cut = z * mask
    y = 0.5 * (d + sqrt(tan(th)**2 + 1 - 4 * z_cut**2) + tan(th))

    return y * mask

def rho_doughnut(z, l, phi, th):
    return pi * y_doughnut(z, l, phi, th)**2 / (l * l)

# Formulas in case of 180º degree contact angle
def phi_min_doughnut_180(l):
    return pi / (6 * l * l)

def phi_max_doughnut_180(l):
    return pi * (0.25 + (pi - 4) / (8 * l) + (10 - 3 * pi) / (24 * l * l))

def S_doughnut_180(l, phi):
    phi_min = phi_min_doughnut_180(l)
    phi_max = phi_max_doughnut_180(l)

    if max(0, phi_min) < phi <= max(0, phi_max):
        a = -pi**3 / 16 + 2 * pi / 3 + 2 * phi * l * l
        b = sqrt(3 * pi**3 * (3 * pi**3 - 32 * pi + 192 * phi * l * l))

        return a + b / 48

    return BIG_NUM

def rho_doughnut_180(z, l, phi):
    # Equations from article (incorrect)
    # z_cut = z * (abs(z) < 0.5)
    # a = sqrt(pi * (3 * pi**3 + 192 * phi * l * l - 32 * pi) / 3) - pi * pi
    # b = 4 * pi * sqrt(1 - 4 * z_cut**2)

    # return (a + b) / (8 * l * l) * (abs(z) < 0.5)

    z_cut = z * (abs(z) < 0.5)
    a = (sqrt(pi * (3 * pi**3 + 192 * phi * l * l - 32 * pi) / 3) - pi * pi) / (4 * pi)
    b = sqrt(1 - 4 * z_cut**2)

    return 0.25 * pi * (a + b)**2 / (l * l) * (abs(z) < 0.5)



''' Worm '''
# General formulas
def r_worm(l, phi, th):
    return sqrt(2 * phi * l / (2 * th - sin(2 * th)))

def phi_min_worm(l, th):
    return 0

def phi_max_worm(l, th):
    return (2 * th - sin(2 * th)) / (2 * l * (1 - cos(th))**2)

def S_worm(l, phi, th):
    phi_min = phi_min_worm(l, th)
    phi_max = phi_max_worm(l, th)

    if max(0, phi_min) < phi <= max(0, phi_max):
        r = r_worm(l, phi, th)

        return 2 * l * r * (th + sin(th))

    return BIG_NUM

def y_worm(z, l, phi, th, center=True):
    r = r_worm(l, phi, th)
    z_wc = 4 * sin(th)**3 * r / (3 * (2 * th - sin(2 * th))) if center else 0
    mask = np.logical_and(cos(th) < (z + z_wc) / r, (z + z_wc) / r < 1)
    z_cut = z * mask
    y = sqrt(r * r - (z_cut + z_wc)**2)

    return y * mask

def rho_worm(z, l, phi, th, center=True):
    return 2 * y_worm(z, l, phi, th, center=center) / l

# Formulas in case of 180º degree contact angle
def phi_min_worm_180(l):
    return 0

def phi_max_worm_180(l):
    return pi / (4 * l)

def S_worm_180(l, phi):
    phi_min = phi_min_worm_180(l)
    phi_max = phi_max_worm_180(l)

    if max(0, phi_min) < phi <= max(0, phi_max):
        return 2 * sqrt(pi * phi * l * l * l)

    return BIG_NUM

def rho_worm_180(z, l, phi):
    mask = abs(z) <= sqrt(phi * l / pi)
    z_cut = z *  mask
    rho = 2 * sqrt(phi * l / pi - z_cut**2) / l

    return rho * mask



''' Roll '''
# General formulas
def d_roll(l, phi, th):
    return 0.25 * (4 * phi * l + pi - 2 * th - 2 * tan(th) + (pi - 2 * th) * tan(th)**2)

def phi_min_roll(l, th):
    return -0.25 * (pi - 2 * th) / (l * cos(th)**2) + 0.5 * tan(th) / l

def phi_max_roll(l, th):
    return 1 + 1 / (l * cos(th)) - 0.25 * (pi - 2 * th) / (l * cos(th)**2) - 0.5 * np.tan(th) / l

def S_roll(l, phi, th):
    phi_min = phi_min_roll(l, th)
    phi_max = phi_max_roll(l, th)

    if max(0, phi_min) < phi <= max(0, phi_max):
        d = d_roll(l, phi, th)

        return 2 * l * d + l * (pi - 2 * th) / cos(th)

    return BIG_NUM

def y_roll(z, l, phi, th):
    d = d_roll(l, phi, th)
    mask = abs(z) < 0.5
    z_cut = z * mask
    y = 0.5 * (d + sqrt(tan(th)**2 + 1 - 4 * z_cut**2) + tan(th))

    return y * mask

def rho_roll(z, l, phi, th):
    return 2 * y_roll(z, l, phi, th) / l

# Formulas in case of 180º degree contact angle
def phi_min_roll_180(l):
    return 0.25 * pi / l

def phi_max_roll_180(l):
    return 1 + 0.25 * (pi - 4) / l

def S_roll_180(l, phi):
    phi_min = phi_min_roll_180(l)
    phi_max = phi_max_roll_180(l)

    if max(0, phi_min) < phi <= max(0, phi_max):
        return 0.5 * pi * l + 2 * phi * l * l

    return BIG_NUM

def rho_roll_180(z, l, phi):
    mask = abs(z) < 0.5
    z_cut = z * mask
    rho = 0.25 * (4 * phi * l + 4 * sqrt(1 - 4 * z_cut**2) - pi) / l

    return rho * mask



''' Perforation '''
# General formulas
def d_perforation(l, phi, th):
    a = (pi - 2 * th)**2 / cos(th)**4
    b = (12 - 4 * (pi - 2 * th) * tan(th)) / cos(th)**2
    c = 64 * (1 - phi) * l * l / pi + 4 / 3
    d = 0.5 * tan(th) - 0.25 * (pi - 2 * th) / cos(th)**2

    return 0.25 * sqrt(a - b + c) + d

def phi_min_perforation(l, th):
    a = 1 - 0.25 * pi + pi / (12 * l * l)
    b = 2 + (pi - 2 * th) * (l - tan(th))

    return a + 0.25 * pi * tan(th) / l - pi * b / (8 * l * l * cos(th)**2)

def phi_max_perforation(l, th):
    c = 1 + pi / (48 * l * l)
    e = cos(3 * th) - 29 * cos(th) + 8 * (pi - 2 * th + sin(2 * th))

    return c + pi * e / (64 * l * l * cos(th)**3)

def S_perforation(l, phi, th):
    phi_min = phi_min_perforation(l, th)
    phi_max = phi_max_perforation(l, th)

    if max(0, phi_min) < phi <= max(0, phi_max):
        d = d_perforation(l, phi, th)

        return 2 * l * l - 0.5 * pi * d * d + pi / cos(th) + 0.5 * pi * (pi - 2 * th) * (d - tan(th)) / cos(th)

    return BIG_NUM

def y_perforation(z, l, phi, th):
    d = d_perforation(l, phi, th)
    mask = abs(z) < 0.5
    z_cut = z * mask
    y = 0.5 * (d - sqrt(tan(th)**2 + 1 - 4 * z_cut**2) - tan(th))

    return y * mask

def rho_perforation(z, l, phi, th):
    rho = 1 - pi * y_perforation(z, l, phi, th)**2 / (l * l)

    return rho * (abs(z) < 0.5)

# Formulas in case of 180º degree contact angle
def phi_min_perforation_180(l):
    return 1 - 0.25 * pi + pi * pi / (8 * l) - pi / (6 * l * l)

def phi_max_perforation_180(l):
    return 1 - pi * (5 / 12 - pi / 8) / (l * l)

def S_perforation_180(l, phi):
    phi_min = phi_min_perforation_180(l)
    phi_max = phi_max_perforation_180(l)

    if max(0, phi_min) < phi <= max(0, phi_max):
        a = pi**3 / 16 - 2 * pi / 3 + 2 * phi * l * l
        b = sqrt(3 * pi**3 * (3 * pi**3 - 32 * pi + 192 * (1 - phi) * l * l))

        return a + b / 48

    return BIG_NUM

def rho_perforation_180(z, l, phi):
    mask = abs(z) < 0.5
    z_cut = z * mask
    a = sqrt(3 * pi**3 + 192 * (1 - phi) * l * l - 32 * pi)
    b = sqrt(3 * pi) * (pi - 4 * sqrt(1 - 4 * z_cut**2))
    rho = 1 - (a + b)**2 / (192 * l * l)

    return rho * mask



''' Layer '''
# General formulas
def phi_min_layer(l, th):
    return 0

def phi_max_layer(l, th):
    return 1

def S_layer(l, phi, th):
    phi_min = phi_min_layer(l, th)
    phi_max = phi_max_layer(l, th)

    if max(0, phi_min) < phi <= max(0, phi_max):
        return 2 * l * l

    return BIG_NUM

def y_layer(z, l, phi, th):
    mask = abs(z) < 0.5 * phi

    return np.ones_like(z) * l * l * mask

def rho_layer(z, l, phi, th):
    return np.ones_like(z) * (abs(z) < 0.5 * phi)

# Formulas in case of 180º degree contact angle
def phi_min_layer_180(l):
    return 0

def phi_max_layer_180(l):
    return 1

def S_layer_180(l, phi):
    phi_min = phi_min_layer_180(l)
    phi_max = phi_max_layer_180(l)

    if max(0, phi_min) < phi <= max(0, phi_max):
        return 2 * l * l

    return BIG_NUM

def rho_layer_180(z, l, phi):
    return np.ones_like(z) * (abs(z) < 0.5 * phi)


