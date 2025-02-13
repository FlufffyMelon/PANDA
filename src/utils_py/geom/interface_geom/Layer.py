from typing import Union

import numpy as np

from .Interface import Interface


class Layer(Interface):
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
        super().__init__(center, borders, l, phi, theta, delta, "layer", extention)

    def check_point(self, point) -> bool:
        y_max = self.y(point[2], self.l, self.phi, *self.attr)

        return (np.abs(point[0]) <= y_max) and (np.abs(point[1]) <= y_max)

    def generate_point(self) -> np.array:
        y_max = self.y(0, self.l, self.phi, *self.attr)
        d = self.delta if self.extention in ["delta", "alpha"] else 0

        count = 0
        while count < 1000:
            point = np.array(
                [
                    np.random.uniform(-y_max, y_max),
                    np.random.uniform(-y_max, y_max),
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
