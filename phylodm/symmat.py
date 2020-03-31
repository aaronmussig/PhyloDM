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

import h5py
import numpy as np

from phylodm.indices import Indices


class SymMat(object):

    def __init__(self, n_indices: int, d_type, arr_default=0):
        self._n_indices = n_indices
        self._d_type = d_type
        self._arr_default = arr_default

        self._indices = Indices()
        self._data = np.full(n_indices * (n_indices + 1) // 2, arr_default)

    def __eq__(self, other):
        if isinstance(other, SymMat):
            return self._n_indices == other._n_indices and \
                   self._d_type == other._d_type and \
                   self._indices == other._indices and \
                   np.array_equal(self._data, other._data)
        return False

    def _data_idx(self, key_i: str, key_j: str) -> int:
        i = self._indices.get_key_idx(key_i)
        j = self._indices.get_key_idx(key_j)
        if i <= j:
            return int(i * self._n_indices - (i - 1) * i / 2 + j - i)
        else:
            return int(j * self._n_indices - (j - 1) * j / 2 + i - j)

    def get_value(self, key_i: str, key_j: str):
        data_idx = self._data_idx(key_i, key_j)
        return self._data[data_idx]

    def set_value(self, key_i: str, key_j: str, value):
        if not self._indices.contains(key_i):
            self._indices.add_key(key_i)
        if not self._indices.contains(key_j):
            self._indices.add_key(key_j)
        self._data[self._data_idx(key_i, key_j)] = value

    def as_matrix(self):
        out = np.full((self._n_indices, self._n_indices), 0)
        out[np.triu_indices(self._n_indices)] = self._data
        diag = out[np.diag_indices(self._n_indices)]
        out = out + out.T
        out[np.diag_indices(self._n_indices)] = out[np.diag_indices(self._n_indices)] - diag
        return out


    def write(self, path: str):
        with h5py.File(path, 'w') as f:
            f.create_dataset('indices',
                             data=[t.encode('ascii') for t in self._indices.get_keys()],
                             dtype=h5py.string_dtype(encoding='ascii'))
            f.create_dataset('data', data=self._data, chunks=True)

    def read(self, path: str):
        self._indices = Indices()
        with h5py.File(path, 'r') as hf:
            self._data = hf['data'][()]

            for idx, key in enumerate(hf['indices'][()]):
                self._indices.add_key(key.decode('utf-8'))
