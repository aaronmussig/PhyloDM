import cython
cimport numpy as np
import numpy as np

ctypedef unsigned long long ULLong
ctypedef unsigned int UInt


cpdef ULLong row_idx_from_mat_coords(ULLong n, ULLong i, ULLong j):
    if i <= j:
        return <ULLong> (i * n - (i - 1) * i / 2 + j - i)
    else:
        return <ULLong> (j * n - (j - 1) * j / 2 + i - j)


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef UInt n_cartesian(UInt[:] groupings):
    cdef UInt total, i, i_start, i_end, i_size, j, n_groups

    n_groups = groupings.shape[0] - 1
    if n_groups == 1:
        return groupings[1]

    total = 0
    for i in range(n_groups):
        i_start = groupings[i]
        i_end = groupings[i + 1]
        i_size = i_end - i_start
        for j in range(i + 1, n_groups):
            total += i_size * (groupings[j+1] - groupings[j])
    return total


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef cartesian_sum(ULLong n_indices, UInt[:] groupings, UInt[:] group_idxs, double[:] group_vals, double[:] data):

    cdef UInt n_comb = n_cartesian(groupings)
    cdef np.ndarray[np.uint64_t, ndim=1] idx_arr = np.empty(n_comb, dtype=np.uint64)
    cdef np.ndarray[np.float64_t, ndim=1] idx_val = np.empty(n_comb, dtype=np.float64)

    cdef UInt i, i_start_idx, i_end_idx, j, j_start_idx, j_end_idx, cur_i, tax_i, cur_j, n_iter, n_groups
    cdef double dist_i, dist_j

    n_iter = 0
    for i in range(groupings.shape[0] - 1):
        i_start_idx = groupings[i]
        i_end_idx = groupings[i+1]

        for j in range(i):
            j_start_idx = groupings[j]
            j_end_idx = groupings[j + 1]

            for cur_i in range(i_start_idx, i_end_idx):
                tax_i = group_idxs[cur_i]
                dist_i = group_vals[cur_i]
                for cur_j in range(j_start_idx, j_end_idx):
                    data[row_idx_from_mat_coords(n_indices, tax_i, group_idxs[cur_j])] = dist_i + group_vals[cur_j]
                    n_iter += 1
    return n_iter
