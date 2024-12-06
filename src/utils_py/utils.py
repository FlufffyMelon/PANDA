from typing import Union
import numpy as np


def get_numerical_density_profile(
    positions: np.array, box: np.array, sl: int, axis: str = "Z", center: bool = True
):
    # TODO: improve it
    assert axis == "Z", "Axis can only be 'X', 'Y' and 'Z'"
    positions, box = validate_positions_and_box(positions, box)

    axis, dz = np.linspace(
        0, np.mean(box, axis=0)[0, 2], sl + 1, endpoint=True, retstep=True
    )

    indexes = np.floor(positions[:, :, 2] / dz).astype(int)
    assert np.max(indexes) < sl, ""
    density_profiles = np.array([np.bincount(i, minlength=sl) for i in indexes])
    density_profiles = np.mean(density_profiles, axis=0)
    density_profiles /= np.mean(box, axis=0)[0, 0] * np.mean(box, axis=0)[0, 1] * dz

    if center:
        axis -= axis[-1] / 2

    return axis[:-1], density_profiles


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
    assert (
        (len(positions.shape) == 2) or (len(positions.shape) == 3)
    ), f"The array of positions must be either two or three dimensional. Now {len(positions.shape)}"

    if len(positions.shape) == 2:
        positions = positions[np.newaxis, :, :]

    return positions


def validate_box(box: np.array):
    assert (
        (len(box.shape) == 1) or (len(box.shape) == 2) or (len(box.shape) == 3)
    ), f"The box array must be either one or two dimensional. Now {len(box.shape)}"

    if len(box.shape) == 1:
        box = box[np.newaxis, np.newaxis, :]
    elif len(box.shape) == 2:
        box = box[:, np.newaxis, :]

    assert box.shape[2] == 3, "Box should be orthoganal"

    return box


def validate_positions_and_box(positions: np.array, box: np.array):
    positions = validate_positions(positions)
    box = validate_box(box)

    assert (
        box.shape[0] == positions.shape[0]
    ), f"Lenghts of box and positions arrays are not the same. Box {box.shape[0]}, when Positions {positions.shape[0]}"

    return positions, box

def validate_list_and_array(array: Union[list, np.array]):
    if isinstance(array, list):
        return np.array(array)
    elif isinstance(array, np.ndarray):
        return array
