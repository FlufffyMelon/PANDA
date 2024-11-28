from typing import Union

import numpy as np

from .Shape import Shape


class Box(Shape):
    def __init__(self, center: Union[list, np.array], borders: Union[list, np.array]):
        super().__init__()

        if isinstance(center, list):
            self.center = np.array(center)
        elif isinstance(center, np.ndarray):
            self.center = center

        if isinstance(borders, list):
            self.borders = np.array(borders)
        elif isinstance(borders, np.ndarray):
            self.borders = borders

    def get_volume(self) -> float:
        return np.prod(self.borders)

    def get_surface(self) -> float:
        return 2 * (
            self.borders[0] * self.borders[1]
            + self.borders[1] * self.borders[2]
            + self.borders[0] * self.borders[2]
        )

    def check_point(self, point) -> bool:
        return np.prod(np.abs(point - self.center) < self.borders / 2)

    def generate_point(self) -> np.array:
        return self.borders * np.random.uniform(-0.5, 0.5, 3) + self.center
