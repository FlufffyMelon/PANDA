import numpy as np

from ..utils import validate_list_and_array
from .Shape import Shape


class UnionShape(Shape):
    def __init__(self, shapeA: Shape, shapeB: Shape):
        super().__init__()

        self.shapeA = shapeA
        self.shapeB = shapeB

    def get_volume(self) -> float:
        return self.shapeA.get_volume() + self.shapeB.get_volume()

    def get_surface(self) -> float:
        return self.shapeA.get_surface() + self.shapeB.get_surface()

    def check_point(self, point) -> bool:
        return self.shapeA.check_point(point) | self.shapeB.check_point(point)

    def generate_point(self) -> np.array:
        volumeA = self.shapeA.get_volume()
        volumeB = self.shapeB.get_volume()
        volume = volumeA + volumeB

        points = [self.shapeA.generate_point(), self.shapeB.generate_point()]
        return points[
            np.random.choice(
                2,
                p=[volumeA / volume, volumeB / volume],
            )
        ]
