use numpy::ndarray::Array2;

pub fn row_idx_from_mat_coords(n: i64, i: i64, j: i64) -> usize {
    /*
    Return the row vector index corresponding to the symmetric matrix coordinates (i, j).
     */
    return if i <= j {
        i * n - (i - 1) * i / 2 + j - i
    } else {
        j * n - (j - 1) * j / 2 + i - j
    } as usize;
}

pub fn row_vec_to_symmat(row_vec: &Vec<f64>, array: &mut Array2<f64>) {
    /*
    Convert a row vector to a symmetric matrix.
     */
    let n_leaves = array.shape()[0];
    for i in 0..n_leaves {
        for j in 0..n_leaves {
            array[[i, j]] = row_vec[row_idx_from_mat_coords(n_leaves as i64, i as i64, j as i64)];
        }
    }
}

