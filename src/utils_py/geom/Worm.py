from typing import Union

import numpy as np

from .Interface import Interface


class Worm(Interface):
    def __init__(
        self,
        center: Union[list, np.array],
        borders: Union[list, np.array],
        l: float,
        phi: float,
        th: float = np.pi,
        delta: float = 0,
        extention: str = None,
        centering: bool = False,
    ):
        super().__init__(center, borders, l, phi, th, delta, "droplet", extention)

        self.centering = centering

    def check_point(self, point) -> bool:
        y_max = self.y(point[2], self.l, self.phi, *self.attr, center=self.centering)

        return np.abs(point[0]) <= y_max

    def generate_point(self) -> np.array:
        angle = self.th if self.extention == "theta" else np.pi

        count = 0
        while count < 1000:
            point = np.array(
                [
                    np.random.uniform(-self.rd, self.rd),
                    np.random.uniform(
                        -0.5 * self.borders[1] / self.H, 0.5 * self.borders[1] / self.H
                    ),
                    np.random.uniform(self.rd * np.cos(angle), self.rd),
                ]
            )

            if self.check_point(point):
                break
            count += 1

        assert count < 1000, "Can't generate point"

        point *= self.H
        point += self.center

        if not self.centering:
            point[2] -= self.center[2]
            point[2] += self.rd * self.H * np.abs(np.cos(angle))
            if self.extention == "delta":
                point[2] += self.delta * self.H

        return point
