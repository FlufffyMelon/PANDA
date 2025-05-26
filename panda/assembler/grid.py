import numpy as np
import itertools
from .utils import distance_pbc2
from tqdm import tqdm


class Grid:
    def __init__(self, box, max_diameter, buffer=0.0, label=None):
        """
        box: simulation box (3,)
        max_diameter: float, the largest molecule diameter
        buffer: float, optional, extra space to add to cell size
        """
        self.box = np.asarray(box)
        self.cell_size = max_diameter + buffer
        self.n_cells = np.floor(self.box / self.cell_size).astype(int)
        label_str = f" ({label})" if label else ""
        print(
            f"[Grid] {self.n_cells[0]}x{self.n_cells[1]}x{self.n_cells[2]}{label_str}"
        )
        self.n_cells[self.n_cells < 1] = 1
        self.cell_size = self.box / self.n_cells  # recalc to fit box exactly
        self.grid = np.empty(self.n_cells, dtype=object)
        self.neighbor_cell_indices = np.empty(self.n_cells, dtype=object)

        total_cells = np.prod(self.n_cells)  # Calculate total number of cells
        # Iterate over all cells and initialize grid and neighbor indices
        for ix, iy, iz in tqdm(
            itertools.product(
                range(self.n_cells[0]), range(self.n_cells[1]), range(self.n_cells[2])
            ),
            total=total_cells,
            desc="Cells",
            leave=False,
        ):
            self.grid[ix, iy, iz] = []

            self.neighbor_cell_indices[ix, iy, iz] = []
            for dx in [-1, 0, 1]:
                for dy in [-1, 0, 1]:
                    for dz in [-1, 0, 1]:
                        nidx = (
                            (ix + dx) % self.n_cells[0],
                            (iy + dy) % self.n_cells[1],
                            (iz + dz) % self.n_cells[2],
                        )
                        self.neighbor_cell_indices[ix, iy, iz].append(nidx)

    def get_cell_index(self, pos):
        idx = np.floor(pos / self.cell_size).astype(int) % self.n_cells
        return tuple(idx)

    def add(self, pos, obj_id):
        idx = self.get_cell_index(pos)
        self.grid[idx].append(obj_id)

    def check_collision(self, pos, all_positions, min_dist2):
        idx = self.get_cell_index(pos)
        for nidx in self.neighbor_cell_indices[idx]:
            ids = self.grid[nidx]

            if ids:
                if np.any(distance_pbc2(pos, all_positions[ids], self.box) < min_dist2):
                    return True
        return False

    def remove(self, pos, obj_id):
        idx = self.get_cell_index(pos)
        if obj_id in self.grid[idx]:
            self.grid[idx].remove(obj_id)

    def move(self, pos_old, pos_new, obj_id):
        idx_old = self.get_cell_index(pos_old)
        idx_new = self.get_cell_index(pos_new)
        if idx_old != idx_new:
            if obj_id in self.grid[idx_old]:
                self.grid[idx_old].remove(obj_id)
            self.grid[idx_new].append(obj_id)
