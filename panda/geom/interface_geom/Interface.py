from typing import Union

import numpy as np

from ... import interface
from ...utils import validate_list_and_array
from ..Shape import Shape


class Interface(Shape):
    def __init__(
        self,
        center: Union[list, np.array],
        borders: Union[list, np.array],
        l: float,
        phi: float,
        theta: float = np.pi,
        delta: float = 0,
        interface_type: str = None,
        extention: str = None,
        debug: bool = False
    ):
        super().__init__()

        self.center = validate_list_and_array(center)
        self.borders = validate_list_and_array(borders)

        self.phi = phi
        self.l = l
        self.H = borders[0] / l
        self.theta = theta
        self.delta = delta

        if extention is None:
            assert False, "Specify the type of the formulas that will be used"
        else:
            assert extention.lower() in [
                "pi",
                "theta",
                "delta",
                "alpha"
            ], "No such type of formulas"
            self.extention = extention.lower()

        if self.extention == "pi":
            self.attr = []
        elif self.extention == "theta":
            self.attr = [theta]
        elif self.extention == "delta":
            self.attr = [delta]
        elif self.extention == "alpha":
            self.attr = [theta, delta]

        if type is None:
            assert False, "Specify the type of the interface"
        else:
            assert interface_type.lower() in [
                "droplet",
                "doughnut",
                "worm",
                "roll",
                "perforation",
                "layer",
            ], "No such interface type"
            interface_type = interface_type.lower()

            self.y = getattr(interface, f"y_{interface_type}_{extention}")
            self.phi_min = getattr(interface, f"phi_min_{interface_type}_{extention}")(
                self.l, *self.attr
            )
            self.phi_max = getattr(interface, f"phi_max_{interface_type}_{extention}")(
                self.l, *self.attr
            )
            self.S = getattr(interface, f"S_{interface_type}_{extention}")(
                self.l, self.phi, *self.attr
            )

            if interface_type in ["droplet", "worm"]:
                self.rd = getattr(interface, f"r_{interface_type}_{extention}")(
                    self.l, self.phi, *self.attr
                )
            elif interface_type in ["doughnut", "roll", "perforation"]:
                self.rd = getattr(interface, f"d_{interface_type}_{extention}")(
                    self.l, self.phi, *self.attr
                )
            else:
                self.rd = None

        if not self.check_existence() and debug:
            print(f"Warning: This type of surface cannot exist under these conditions (l {self.l:.2f}; phi {self.phi:.2f}; theta {self.theta:.2f}; delta {self.delta:.2f}) ! (phi_min {self.phi_min:.2f}; phi_max {self.phi_max:.2f})")

    def get_volume(self) -> float:
        return self.phi * np.prod(self.borders)

    def get_surface(self) -> float:
        return self.S

    def check_existence(self) -> bool:
        return self.phi_min <= self.phi <= self.phi_max
