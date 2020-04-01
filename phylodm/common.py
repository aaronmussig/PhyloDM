import numpy as np


def create_mat(n: int, default) -> np.array:
    return np.full(n * (n + 1) // 2, default)


def row_idx_from_mat_coords(n: int, i: int, j: int) -> int:
    if i <= j:
        return int(i * n - (i - 1) * i / 2 + j - i)
    else:
        return int(j * n - (j - 1) * j / 2 + i - j)


def mat_shape_from_row_shape(n: int) -> int:
    return int(0.5 * (np.sqrt(8 * n + 1) - 1))
