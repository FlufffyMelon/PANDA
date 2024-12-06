from typing import Union

import numpy as np

from .. import auxil
from ..utils import validate_list_and_array
from .Shape import Shape


class Interface(Shape):
    def __init__(
        self,
        center: Union[list, np.array],
        borders: Union[list, np.array],
        l: float,
        phi: float,
        th: float = np.pi,
        delta: float = 0,
        interface_type: str = None,
        extention: str = None,
    ):
        super().__init__()

        self.center = validate_list_and_array(center)
        self.borders = validate_list_and_array(borders)

        self.phi = phi
        self.l = l
        self.H = borders[0] / l
        self.th = th
        self.delta = delta

        if extention is None:
            assert False, "Specify the type of the formulas that will be used"
        else:
            assert extention.lower() in [
                "pi",
                "theta",
                "delta",
            ], "No such type of formulas"
            self.extention = extention.lower()

        if self.extention == "pi":
            self.attr = []
        elif self.extention == "theta":
            self.attr = [th]
        elif self.extention == "delta":
            self.attr = [delta]

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

            self.y = getattr(auxil, f"y_{interface_type}_{extention}")
            self.phi_min = getattr(auxil, f"phi_min_{interface_type}_{extention}")(
                self.l, *self.attr
            )
            self.phi_max = getattr(auxil, f"phi_max_{interface_type}_{extention}")(
                self.l, *self.attr
            )
            self.S = getattr(auxil, f"S_{interface_type}_{extention}")(
                self.l, self.phi, *self.attr
            )

            if interface_type in ["droplet", "worm"]:
                self.rd = getattr(auxil, f"r_{interface_type}_{extention}")(
                    self.l, self.phi, *self.attr
                )
            elif interface_type in ["doughnut", "roll", "perforation"]:
                self.rd = getattr(auxil, f"d_{interface_type}_{extention}")(
                    self.l, self.phi, *self.attr
                )
            else:
                self.rd = None

        assert (
            self.check_existence()
        ), "This type of surface cannot exist under these conditions"

    def get_volume(self) -> float:
        return self.phi * np.prod(self.borders)

    def get_surface(self) -> float:
        return self.S

    def check_existence(self) -> bool:
        return self.phi_min <= self.phi <= self.phi_max
