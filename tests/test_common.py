import unittest

from phylodm.common import *


class TestCommon(unittest.TestCase):

    def test_row_idx_from_mat_coords(self):
        n = 1000
        expected = np.triu_indices(n)
        for true_idx, (i, j) in enumerate(zip(expected[0], expected[1])):
            test_idx_a = row_idx_from_mat_coords(n, i, j)
            test_idx_b = row_idx_from_mat_coords(n, j, i)
            self.assertTrue(true_idx == test_idx_a == test_idx_b)

    def test_mat_shape_from_row_shape(self):
        for n in range(1000):
            row_len = len(np.triu_indices(n)[0])
            test = mat_shape_from_row_shape(row_len)
            self.assertTrue(n == test)

