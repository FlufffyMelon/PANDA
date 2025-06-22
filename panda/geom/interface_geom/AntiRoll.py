from typing import Union

import numpy as np

from . import Roll


class Antiroll(Roll):
    def __init__(
        self,
        center: Union[list, np.array],
        borders: Union[list, np.array],
        l: float,
        phi: float,
        theta: float = np.pi,
        delta: float = 0,
        extention: str = None,
    ):
        super().__init__(center, borders, l, phi, theta, delta, extention)

    def get_volume(self) -> float:
        return (1 - self.phi) * np.prod(self.borders)

    def check_point(self, point) -> bool:
        return not super().check_point(point)

    def generate_point(self) -> np.array:
        count = 0
        while count < 1000:
            point = np.array(
                [
                    np.random.uniform(-self.l / 2, self.l / 2),
                    np.random.uniform(-self.l / 2, self.l / 2),
                    np.random.uniform(-0.5, 0.5),
                ]
            )

            if self.check_point(point):
                break
            count += 1

        assert count < 1000, "Can't generate point"

        point *= self.H
        point += self.center

        return point
