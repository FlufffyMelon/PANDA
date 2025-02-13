from typing import Union

import numpy as np

from .Interface import Interface


class Perforation(Interface):
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
        super().__init__(center, borders, l, phi, theta, delta, "perforation", extention)

    def check_point(self, point) -> bool:
        y_min = self.y(point[2], self.l, self.phi, *self.attr)

        return (point[0] ** 2 + point[1] ** 2) >= y_min**2

    def generate_point(self) -> np.array:
        d = self.delta if self.extention in ["delta", "alpha"] else 0

        count = 0
        while count < 1000:
            point = np.array(
                [
                    np.random.uniform(
                        -0.5 * self.borders[0] / self.H, 0.5 * self.borders[0] / self.H
                    ),
                    np.random.uniform(
                        -0.5 * self.borders[1] / self.H, 0.5 * self.borders[1] / self.H
                    ),
                    np.random.uniform(-0.5 + d, 0.5 - d),
                ]
            )

            if self.check_point(point):
                break
            count += 1

        assert count < 1000, "Can't generate point"

        point *= self.H
        point += self.center

        return point
