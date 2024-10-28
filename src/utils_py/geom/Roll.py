from abc import ABC, abstractmethod
from dataclasses import dataclass, field
import numpy as np
from ..auxil_theta import S_roll, y_roll, phi_min_roll, phi_max_roll
from .Shape import Shape

@dataclass
class Roll(Shape):
    borders: np.array = field(default_factory=np.array)
    phi: float = 0
    th: float = 0
    H: float = 1

    def get_volume(self) -> float:
        return self.phi * np.prod(self.borders)

    def get_surface(self) -> float:
        return 42

    def check_point(self, point) -> bool:
        return np.abs(point[0] / self.H) <= y_roll(point[2] / self.H, self.borders[0] / self.H, self.phi, self.th)

    def generate_point(self) -> np.array:
        if not (phi_min_roll(self.borders[0] / self.H, self.th) <= self.phi <= phi_max_roll(self.borders[0] / self.H, self.th)):
            print('!!!roll shape doesn`t exist, try another one!!!')
            return

        # Y_max = (self.y(0) if (self.theta > np.pi/2) else self.y(0.5)) * self.borders[2]
        y_max = y_roll(0, self.borders[0] / self.H, self.phi, self.th) * self.H
        point = np.array([np.random.uniform(-y_max, y_max),
                          np.random.uniform(-self.borders[1]/2, self.borders[1]/2),
                          np.random.uniform(-self.borders[2]/2, self.borders[2]/2)])

        while not self.check_point(point):
           point = np.array([np.random.uniform(-y_max, y_max),
                          np.random.uniform(-self.borders[1]/2, self.borders[1]/2),
                          np.random.uniform(-self.borders[2]/2, self.borders[2]/2)])

        return point + self.center
