import numpy as np
from .utils import distance_pbc2
from panda.utils import validate_box, validate_positions_and_box


class Grid:
    def __init__(self, box, max_diameter, buffer=0.0):
        """
        box: simulation box (3,)
        max_diameter: float, the largest molecule diameter
        buffer: float, optional, extra space to add to cell size
        """
        self.box = np.asarray(box)
        self.cell_size = max_diameter + buffer
        self.n_cells = np.floor(self.box / self.cell_size).astype(int)
        print(f"Constructing grid with dimensions {'x'.join(map(str, self.n_cells))}")
        self.n_cells[self.n_cells < 1] = 1
        self.cell_size = self.box / self.n_cells  # recalc to fit box exactly
        self.grid = np.empty(self.n_cells, dtype=object)
        self.neighbor_cell_indices = np.empty(self.n_cells, dtype=object)
        for ix in range(self.n_cells[0]):
            for iy in range(self.n_cells[1]):
                for iz in range(self.n_cells[2]):
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
                # neigh_pos = np.array([all_positions[i] for i in ids])
                if np.any(distance_pbc2(pos, all_positions[ids], self.box) < min_dist2):
                    return True
        return False
