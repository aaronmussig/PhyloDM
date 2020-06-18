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

from typing import Union

import numpy as np


def create_mat_vector(n: int, default: Union[int, float]) -> np.array:
    """Create a vector which represents the upper/lower triangle of a matrix."""
    return np.full(n * (n + 1) // 2, default)


def compact_int_mat(mat: np.array) -> np.array:
    """For a matrix of integers, return a copy with the smallest data type."""
    signed = (('int8', -128, 127),
              ('int16', -32768, 32767),
              ('int32', -2147483648, 2147483647),
              ('int64', -9223372036854775808, 9223372036854775807))
    unsigned = (('uint8', 0, 255),
                ('uint16', 0, 65535),
                ('uint32', 0, 4294967295),
                ('uint64', 0, 18446744073709551615))
    m_min, m_max = mat.min(), mat.max()
    if m_min < 0:
        options = signed
    else:
        options = unsigned
    for d_type, d_min, d_max in options:
        if m_min >= d_min and m_max <= d_max:
            return mat.astype(np.dtype(d_type))
    return mat
