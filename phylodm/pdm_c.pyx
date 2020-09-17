#!python
#cython: language_level=3

import cython


ctypedef unsigned char uchar
ctypedef unsigned short ushort
ctypedef unsigned int uint
ctypedef unsigned long ulong

ctypedef fused fused_dtypes:
    char
    unsigned char
    short
    unsigned short
    int
    unsigned int
    long
    unsigned long
    float
    double


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef Py_ssize_t row_idx_from_mat_coords(Py_ssize_t n, Py_ssize_t i, Py_ssize_t j) nogil:
    """Get the row vector index from symmetric matrix coordinates."""
    if i <= j:
        return <Py_ssize_t> (i * n - (i - 1) * i / 2 + j - i)
    else:
        return <Py_ssize_t> (j * n - (j - 1) * j / 2 + i - j)


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef uint n_cartesian(uint[:] groupings) nogil:
    cdef Py_ssize_t total, i, i_start, i_end, i_size, j, n_groups

    with nogil:
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
cpdef Py_ssize_t cartesian_sum(Py_ssize_t n_indices, uint[:] groupings, uint[:] group_idxs, double[:] group_vals, double[:] data) nogil:

    cdef Py_ssize_t i, i_start_idx, i_end_idx, j, j_start_idx, j_end_idx, cur_i, tax_i, cur_j, n_iter, n_groups
    cdef double dist_i, dist_j

    with nogil:
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




@cython.boundscheck(False)
@cython.wraparound(False)
cdef void row_vec_to_symmat_fused(fused_dtypes[:] vec, fused_dtypes[:, :] mat) nogil:
    """Convert the row vector to a symmetric matrix."""

    cdef Py_ssize_t i, j, n, cur_idx
    with gil:
        n = mat.shape[0]

    # Set each of the elements.
    for i in range(n):
        for j in range(i):
            cur_idx = row_idx_from_mat_coords(n, i, j)
            mat[i, j] = vec[cur_idx]
            mat[j, i] = vec[cur_idx]

        # Set the diagonal.
        mat[i, i] = 0


cpdef row_vec_to_symmat(vec, mat):
    print(vec.dtype.name)
    if vec.dtype.name == 'int8':
        row_vec_to_symmat_fused[char](vec, mat)
    elif vec.dtype.name == 'uint8':
        row_vec_to_symmat_fused[uchar](vec, mat)
    elif vec.dtype.name == 'int16':
        row_vec_to_symmat_fused[short](vec, mat)
    elif vec.dtype.name == 'uint16':
        row_vec_to_symmat_fused[ushort](vec, mat)
    elif vec.dtype.name == 'int32':
        row_vec_to_symmat_fused[int](vec, mat)
    elif vec.dtype.name == 'uint32':
        row_vec_to_symmat_fused[uint](vec, mat)
    elif vec.dtype.name == 'int64':
        row_vec_to_symmat_fused[long](vec, mat)
    elif vec.dtype.name == 'uint64':
        row_vec_to_symmat_fused[ulong](vec, mat)
    elif vec.dtype.name == 'float32':
        row_vec_to_symmat_fused[float](vec, mat)
    elif vec.dtype.name == 'float64':
        row_vec_to_symmat_fused[double](vec, mat)
    else:
        raise TypeError(f'Unknown type: {vec.dtype}')


