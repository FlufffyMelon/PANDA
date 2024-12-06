from typing import Union

import numpy as np

from ..utils import validate_list_and_array
from .Shape import Shape


class AntiBox(Shape):
    def __init__(
        self,
        box_center: Union[list, np.array],
        box_borders: Union[list, np.array],
        center: Union[list, np.array],
        borders: Union[list, np.array],
    ):
        super().__init__()

        self.box_center = validate_list_and_array(box_center)
        self.box_borders = validate_list_and_array(box_borders)
        self.center = validate_list_and_array(center)
        self.borders = validate_list_and_array(borders)

    def get_volume(self) -> float:
        return np.prod(self.borders) - np.prod(self.box_borders)

    def get_surface(self) -> float:
        pass

    def check_point(self, point) -> bool:
        return np.prod(np.abs(point - self.center) < self.borders / 2) and np.prod(
            np.abs(point - self.box_center) > self.box_borders / 2
        )

    def generate_point(self) -> np.array:
        relative_box_center = self.box_center - self.center
        z = np.random.uniform(-0.5, 0.5) * self.borders[2]
        if np.abs(z - relative_box_center[2]) > self.box_borders[2] / 2:
            x, y = np.random.uniform(-0.5, 0.5, 2) * self.borders[:2]
        else:
            y = np.random.uniform(-0.5, 0.5) * self.borders[1]
            if np.abs(y - relative_box_center[1]) > self.box_borders[1] / 2:
                x = np.random.uniform(-0.5, 0.5) * self.borders[0]
            else:
                x = np.random.choice(
                    [
                        np.random.uniform(
                            -self.borders[0] / 2,
                            relative_box_center[0] - self.box_borders[0] / 2,
                        ),
                        np.random.uniform(
                            relative_box_center[0] + self.box_borders[0] / 2,
                            self.borders[0] / 2,
                        ),
                    ]
                )

        return np.array([x, y, z]) + self.center
