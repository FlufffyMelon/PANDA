import numpy as np
from panda.utils import validate_positions_and_box


def delta_pbc(point: np.array, points: np.array, box: np.array) -> np.array:
    delta = points - point
    idx_pbc = np.abs(delta) * 2 >= box
    delta -= np.sign(delta) * box * idx_pbc
    return delta


def distance_pbc2(point: np.array, points: np.array, box: np.array) -> np.array:
    assert len(point.shape) == 1, f"point must be 1D array, got shape {point.shape}"
    assert len(points.shape) == 2, f"points must be 2D array, got shape {points.shape}"
    assert len(box.shape) == 1, f"box must be 1D array, got shape {box.shape}"
    assert point.shape[0] == 3, f"point must have 3 coordinates, got {point.shape[0]}"
    assert points.shape[1] == 3, (
        f"points must have 3 coordinates per point, got {points.shape[1]}"
    )
    assert box.shape[0] == 3, f"box must have 3 dimensions, got {box.shape[0]}"

    return np.sum(delta_pbc(point, points, box) ** 2, axis=1)


def distance_pbc(point1: np.array, points: np.array, box: np.array) -> np.array:
    return np.sqrt(distance_pbc2(point1, points, box))


# def apply_pbc_to_points(points: np.array, box_size: np.array, mask: np.array):
#     assert len(point.shape) == 1, f"point must be 1D array, got shape {point.shape}"
#     assert len(points.shape) == 2, f"points must be 2D array, got shape {points.shape}"
#     assert len(box.shape) == 1, f"box must be 1D array, got shape {box.shape}"
#     assert point.shape[0] == 3, f"point must have 3 coordinates, got {point.shape[0]}"
#     assert points.shape[1] == 3, (
#         f"points must have 3 coordinates per point, got {points.shape[1]}"
#     )
#     assert box.shape[0] == 3, f"box must have 3 dimensions, got {box.shape[0]}"

#     points_pbc = points.copy()
#     half_box_size = box_size[0:3] / 2

#     idx_pbc = abs(points_pbc - half_box_size) >= half_box_size
#     points_pbc -= np.sign(points_pbc) * box_size[0:3] * idx_pbc * mask

#     return points_pbc
