from abc import ABC, abstractmethod
import numpy as np


class Shape(ABC):
    @abstractmethod
    def get_volume(self) -> float:
        pass

    @abstractmethod
    def get_surface(self) -> float:
        pass

    @abstractmethod
    def check_point(self) -> bool:
        pass

    @abstractmethod
    def generate_point(self) -> np.array:
        pass
