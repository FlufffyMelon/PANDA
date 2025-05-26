import numpy as np


def delta_pbc(point: np.array, points: np.array, box: np.array) -> np.array:
    assert len(point.shape) == 1, f"point must be 1D array, got shape {point.shape}"
    assert len(points.shape) == 2, f"points must be 2D array, got shape {points.shape}"
    assert len(box.shape) == 1, f"box must be 1D array, got shape {box.shape}"
    assert point.shape[0] == 3, f"point must have 3 coordinates, got {point.shape[0]}"
    assert points.shape[1] == 3, (
        f"points must have 3 coordinates per point, got {points.shape[1]}"
    )
    assert box.shape[0] == 3, f"box must have 3 dimensions, got {box.shape[0]}"

    delta = points - point
    idx_pbc = np.abs(delta) * 2 >= box
    delta -= np.sign(delta) * box * idx_pbc
    return delta


def distance_pbc2(point: np.array, points: np.array, box: np.array) -> np.array:
    return np.sum(delta_pbc(point, points, box) ** 2, axis=1)


def distance_pbc(point: np.array, points: np.array, box: np.array) -> np.array:
    return np.sqrt(distance_pbc2(point, points, box))
