from typing import Union
import numpy as np
import math


def get_center_pbc(positions: np.array, box: np.array):
    positions, box = validate_positions_and_box(positions, box)

    theta = positions / box * 2 * np.pi
    center = np.zeros((positions.shape[0], 3))

    for i in range(3):
        phi = np.cos(theta[:, :, i])
        psi = np.sin(theta[:, :, i])

        phi_mean = np.mean(phi, axis=1)
        psi_mean = np.mean(psi, axis=1)

        theta_mean = np.arctan2(-psi_mean, -phi_mean) + np.pi
        center[:, i] = np.mean(box, axis=1)[:, i] * theta_mean / 2 / np.pi

    return center[:, np.newaxis, :]


def apply_pbc(positions: np.array, box: np.array):
    positions, box = validate_positions_and_box(positions, box)
    half_box_size = box / 2

    ids = abs(positions - half_box_size) >= half_box_size
    positions -= np.sign(positions) * box * ids

    return positions


def validate_positions(positions: np.array):
    assert (len(positions.shape) == 2) or (len(positions.shape) == 3), (
        f"The array of positions must be either two or three dimensional. Now {len(positions.shape)}"
    )

    if len(positions.shape) == 2:
        positions = positions[np.newaxis, :, :]

    return positions


def validate_box(box: np.array):
    assert (len(box.shape) == 1) or (len(box.shape) == 2) or (len(box.shape) == 3), (
        f"The box array must be either one or two dimensional. Now {len(box.shape)}"
    )

    if len(box.shape) == 1:
        box = box[np.newaxis, np.newaxis, :]
    elif len(box.shape) == 2:
        box = box[:, np.newaxis, :]

    assert box.shape[2] == 3, "Box should be orthoganal"

    return box


def validate_positions_and_box(positions: np.array, box: np.array):
    positions = validate_positions(positions)
    box = validate_box(box)

    assert box.shape[0] == positions.shape[0], (
        f"Lenghts of box and positions arrays are not the same. Box {box.shape[0]}, when Positions {positions.shape[0]}"
    )

    return positions, box


def validate_list_and_array(array: Union[list, np.array]):
    if isinstance(array, list):
        return np.array(array)
    elif isinstance(array, np.ndarray):
        return array
    elif array is None:
        return None
    else:
        return np.array(array)


def str2bool(s: str) -> bool:
    """Helper function to support boolean command line arguments."""
    if s.lower() in {"true", "t", "1"}:
        return True
    elif s.lower() in {"false", "f", "0"}:
        return False
    else:
        raise ValueError(f"Invalid boolean value: {s}")


def format_bytes(size):
    if size == 0:
        return "0 B"
    unit = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
    i = int(math.floor(math.log(size, 1024)))
    p = math.pow(1024, i)
    s = round(size / p, 2)
    return f"{s} {unit[i]}"


def serialize_component_config(components):
    """
    Convert a list of Component objects to a list of JSON-serializable dicts
    with region and insertion_region as dicts.
    """

    def region_to_dict(region):
        if region is None:
            return None
        d = {"type": type(region).__name__}
        if hasattr(region, "center"):
            center = getattr(region, "center")
            if hasattr(center, "tolist"):
                center = center.tolist()
            d["center"] = center
        if hasattr(region, "borders"):
            borders = getattr(region, "borders")
            if hasattr(borders, "tolist"):
                borders = borders.tolist()
            d["borders"] = borders
        return d

    result = []
    for comp in components:
        comp_dict = {
            "name": comp.name,
            "region": region_to_dict(getattr(comp, "region", None)),
            "insertion_region": region_to_dict(getattr(comp, "insertion_region", None)),
            "density": getattr(comp, "density", None),
            "numbers": getattr(comp, "numbers", None),
        }
        result.append(comp_dict)
    return result
