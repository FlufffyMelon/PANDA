from typing import Union

import numpy as np

from .Interface import Interface


class Roll(Interface):
    def __init__(
        self,
        center: Union[list, np.array],
        borders: Union[list, np.array],
        l: float,
        phi: float,
        th: float = np.pi,
        delta: float = 0,
        extention: str = None,
    ):
        super().__init__(center, borders, l, phi, th, delta, "roll", extention)

        # print("l:", self.l)
        # print("phi:", self.phi)
        # print("attr:", self.attr)

    def check_point(self, point) -> bool:
        y_max = self.y(point[2], self.l, self.phi, *self.attr)

        return np.abs(point[0]) <= y_max

    def generate_point(self) -> np.array:
        y_max = self.y(0, self.l, self.phi, *self.attr)
        d = self.delta if self.extention == "delta" else 0
        # print("y_max:", y_max)

        count = 0
        while count < 1000:
            point = np.array(
                [
                    np.random.uniform(-y_max, y_max),
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
