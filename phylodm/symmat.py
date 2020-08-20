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

from typing import Collection, Tuple, Optional, Union

import h5py
import numpy as np
from phylodm.pdm_c import row_idx_from_mat_coords

from phylodm.common import create_mat_vector
from phylodm.indices import Indices


class SymMat(object):

    def __init__(self):
        """Variables are set when called by any of the static methods."""
        self._d_type: Optional[np.dtype] = None
        self._arr_default: Optional[Union[int, float]] = None
        self._indices: Optional[Indices] = None
        self._data: Optional[np.array] = None

    def __eq__(self, other) -> bool:
        """Two SymMats are equal if the data, defaults, and indices are equal."""
        if isinstance(other, SymMat):
            return self._d_type == other._d_type and \
                   self._arr_default == other._arr_default and \
                   self._indices == other._indices and \
                   np.array_equal(self._data, other._data)
        return False

    @staticmethod
    def get_from_shape(n_indices: int, d_type: np.dtype,
                       arr_default: Union[int, float] = 0) -> 'SymMat':
        """Create a blank SymMat given the specifications."""
        return SymMat()._get_from_shape(n_indices, d_type, arr_default)

    @staticmethod
    def get_from_indices(indices: Collection[str], d_type: np.dtype,
                         arr_default: Union[int, float] = 0) -> 'SymMat':
        """Create a blank SymMat given the indices."""
        return SymMat()._get_from_indices(indices, d_type, arr_default)

    @staticmethod
    def get_from_path(path: str) -> 'SymMat':
        """Load the SymMat from a cache."""
        return SymMat()._get_from_path(path)

    def _idx_from_key(self, key_i: str, key_j: str) -> int:
        """Determine the row vector index given the indices names."""
        i = self._indices.get_key_idx(key_i)
        j = self._indices.get_key_idx(key_j)
        return row_idx_from_mat_coords(len(self._indices), i, j)

    def _get_from_shape(self, n_indices: int, d_type: np.dtype,
                        arr_default: Union[int, float] = 0) -> 'SymMat':
        """Create a blank SymMat given the specifications."""
        self._d_type = d_type
        self._arr_default = arr_default
        self._indices = Indices()
        self._data = create_mat_vector(n_indices, arr_default)
        return self

    def _get_from_indices(self, indices: Collection[str], d_type: np.dtype,
                          arr_default: Union[int, float] = 0) -> 'SymMat':
        """Create a blank SymMat given the indices."""
        self._d_type = d_type
        self._arr_default = arr_default
        self._indices = Indices()
        self._indices.add_keys(indices)
        self._data = create_mat_vector(len(indices), arr_default)
        return self

    def _get_from_path(self, path: str) -> 'SymMat':
        """Load the SymMat from a cache."""
        self._indices = Indices()
        with h5py.File(path, 'r') as hf:
            self._arr_default = hf['arr_default'][()]
            self._data = hf['data'][()]
            for idx, key in enumerate(hf['indices'][()]):
                self._indices.add_key(key.decode('utf-8'))
        self._d_type = self._data.dtype
        return self

    def save_to_path(self, path: str):
        """Save the SymMat to a cache."""
        with h5py.File(path, 'w') as f:
            f.create_dataset('indices',
                             data=[t.encode('ascii') for t in self._indices.get_keys()],
                             dtype=h5py.string_dtype(encoding='ascii'))
            f.create_dataset('data', data=self._data, dtype=self._d_type)
            f.create_dataset('arr_default', data=self._arr_default)

    def get_value(self, key_i: str, key_j: str) -> Union[int, float]:
        """Get the value for the specified keys."""
        data_idx = self._idx_from_key(key_i, key_j)
        return self._data[data_idx]

    def set_value(self, key_i: str, key_j: str, value: Union[float, int]):
        """Set a specific value given the keys."""
        self._data[self._idx_from_key(key_i, key_j)] = value

    def remove_keys(self, keys: Collection[str]):
        """Remove the key and associated values from the matrix."""
        key_idx = self._indices.get_key_indices()
        keep_keys = sorted(key_idx.keys() - set(keys))
        new_indices = Indices()
        new_indices.add_keys(keep_keys)

        # Create the new matrix and import all of the keys across.
        new_mat = SymMat.get_from_indices(keep_keys, d_type=self._d_type,
                                          arr_default=self._arr_default)

        # Determine the mapping for importing the data across.
        for i in range(len(keep_keys)):
            for j in range(i + 1):
                new_mat.set_value(keep_keys[i], keep_keys[j],
                                  self.get_value(keep_keys[i], keep_keys[j]))
        self._data = new_mat._data
        self._indices = new_mat._indices

    def as_matrix(self) -> Tuple[Tuple[str], np.array]:
        """Return a symmetric numpy matrix given the SymMat."""
        n_indices = len(self._indices)
        mat = np.full((n_indices, n_indices), 0, dtype=self._d_type)
        mat[np.triu_indices_from(mat)] = self._data
        diag = mat[np.diag_indices_from(mat)]
        mat = mat + mat.T
        mat[np.diag_indices_from(mat)] = mat[np.diag_indices_from(mat)] - diag
        return self._indices.get_keys(), mat
