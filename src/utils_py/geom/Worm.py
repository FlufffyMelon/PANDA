from abc import ABC, abstractmethod
from dataclasses import dataclass, field
import numpy as np
from ..auxil import S_worm, r_worm, y_worm, phi_min_worm, phi_max_worm
from .Shape import Shape

@dataclass
class Worm(Shape):
    borders: np.array = field(default_factory=np.array)
    phi: float = 0
    th: float = 0
    H: float = 1

    def get_volume(self) -> float:
        return self.phi * np.prod(self.borders)

    def get_surface(self) -> float:
        return 42

    def check_point(self, point) -> bool:
        return np.abs(point[0] / self.H) <= y_worm(point[2] / self.H, self.borders[0] / self.H, self.phi, self.th, center=False)

    def generate_point(self) -> np.array:
        if not (phi_min_worm(self.borders[0] / self.H, self.th) <= self.phi <= phi_max_worm(self.borders[0] / self.H, self.th)):
            print('!!!worm shape doesn`t exist, try another one!!!')
            return

        r = r_worm(self.borders[0] / self.H, self.phi, self.th) * self.H
        point = np.array([np.random.uniform(-r, r),
                          np.random.uniform(-self.borders[1]/2, self.borders[1]/2),
                          np.random.uniform(r * np.cos(self.th), r)])

        while not self.check_point(point):
           point = np.array([np.random.uniform(-r, r),
                          np.random.uniform(-self.borders[1]/2, self.borders[1]/2),
                          np.random.uniform(r * np.cos(self.th), r)])

        point[2] -= r * np.cos(self.th)
        return point + self.center
