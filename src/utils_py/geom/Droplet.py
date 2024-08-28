from abc import ABC, abstractmethod
from dataclasses import dataclass, field
import numpy as np
from ..auxil import S_droplet, r_droplet, y_droplet, phi_min_droplet, phi_max_droplet
from .Shape import Shape

@dataclass
class Droplet(Shape):
    borders: np.array = field(default_factory=np.array)
    phi: float = 0
    th: float = 0
    H: float = 1

    def get_volume(self) -> float:
        return self.phi * np.prod(self.borders)

    def get_surface(self) -> float:
        return 42

    def check_point(self, point) -> bool:
        return (point[0]**2 + point[1]**2) / self.H**2 <= y_droplet(point[2] / self.H, self.borders[0] / self.H, self.phi, self.th, center=False)**2

    def generate_point(self) -> np.array:
        if not (phi_min_droplet(self.borders[0] / self.H, self.th) <= self.phi <= phi_max_droplet(self.borders[0] / self.H, self.th)):
            print('!!!droplet shape doesn`t exist, try another one!!!')
            return

        r = r_droplet(self.borders[0] / self.H, self.phi, self.th) * self.H
        point = np.array([np.random.uniform(-r, r),
                          np.random.uniform(-r, r),
                          np.random.uniform(r * np.cos(self.th), r)])

        while not self.check_point(point):
            point = np.array([np.random.uniform(-r, r),
                              np.random.uniform(-r, r),
                              np.random.uniform(r * np.cos(self.th), r)])

        point[2] -= r * np.cos(self.th)
        return point + self.center
