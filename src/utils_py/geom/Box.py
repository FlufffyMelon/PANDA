from typing import Union

import numpy as np

from ..utils import validate_list_and_array
from .Shape import Shape


class Box(Shape):
    def __init__(self, center: Union[list, np.array], borders: Union[list, np.array]):
        super().__init__()

        self.center = validate_list_and_array(center)
        self.borders = validate_list_and_array(borders)

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
