from abc import ABC, abstractmethod
from dataclasses import dataclass, field
import numpy as np
from ..auxil_theta import S_perforation, y_perforation, phi_min_perforation, phi_max_perforation
from .Shape import Shape

@dataclass
class Perforation(Shape):
    borders: np.array = field(default_factory=np.array)
    phi: float = 0
    th: float = 0
    H: float = 1

    def get_volume(self) -> float:
        return self.phi * np.prod(self.borders)

    def get_surface(self) -> float:
        return 42

    def check_point(self, point) -> bool:
        return (point[0]**2 + point[1]**2) / self.H**2 >= y_perforation(point[2] / self.H, self.borders[0] / self.H, self.phi, self.th)**2

    def generate_point(self) -> np.array:
        if not (phi_min_perforation(self.borders[0] / self.H, self.th) <= self.phi <= phi_max_perforation(self.borders[0] / self.H, self.th)):
            print('!!!perforation shape doesn`t exist, try another one!!!')
            return

        # y_max = y_perforation(0, self.borders[0] / self.H, self.phi, self.th) * self.H
        point = np.array([np.random.uniform(-self.borders[0] / 2, self.borders[0] / 2),
                          np.random.uniform(-self.borders[1] / 2, self.borders[1] / 2),
                          np.random.uniform(-self.borders[2] / 2, self.borders[2] / 2)])

        while not self.check_point(point):
            point = np.array([np.random.uniform(-self.borders[0] / 2, self.borders[0] / 2),
                              np.random.uniform(-self.borders[1] / 2, self.borders[1] / 2),
                              np.random.uniform(-self.borders[2] / 2, self.borders[2] / 2)])

        return point + self.center
