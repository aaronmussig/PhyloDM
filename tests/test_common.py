###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import itertools
import unittest

import numpy as np

from phylodm.common import create_mat_vector, compact_int_mat
from phylodm.pdm_c import row_idx_from_mat_coords, n_cartesian

N_TESTS = 100


class TestCommon(unittest.TestCase):

    def test_row_idx_from_mat_coords(self):
        expected = np.triu_indices(N_TESTS)
        for true_idx, (i, j) in enumerate(zip(expected[0], expected[1])):
            test_idx_a = row_idx_from_mat_coords(N_TESTS, i, j)
            test_idx_b = row_idx_from_mat_coords(N_TESTS, j, i)
            self.assertTrue(true_idx == test_idx_a == test_idx_b)

    def test_compact_int_mat(self):
        cases = (('int8', -128),
                 ('int16', -32768),
                 ('int32', -2147483648),
                 ('int64', -922337203685477),
                 ('uint8', 255),
                 ('uint16', 65535),
                 ('uint32', 4294967295),
                 ('uint64', 1844674407370955))
        for case, num in cases:
            true_arr = np.array([[1, 2], [3, num]], dtype=np.dtype('int64'))
            comp_arr = compact_int_mat(true_arr)
            self.assertEqual(case, comp_arr.dtype.name)

    def test_n_cartesian(self):
        for n_groups in range(1, N_TESTS):
            true_list = list()
            test_list = [0]

            # Generate a random number members for each group.
            count = 0
            for cur_group in range(n_groups):
                n_members = np.random.randint(1, 10)
                true_list.append(np.random.randint(1, 100, n_members))
                count += n_members
                test_list.append(count)

            test = n_cartesian(np.array(test_list, dtype=np.uint32))
            if n_groups == 1:
                true = len(true_list[0])
            else:
                true = 0
                for cur_item, next_item in itertools.combinations(true_list, 2):
                    true += np.transpose([np.tile(cur_item, len(next_item)),
                                          np.repeat(next_item, len(cur_item))]).shape[0]
            self.assertEqual(true, test)
