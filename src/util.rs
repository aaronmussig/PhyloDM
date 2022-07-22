
use numpy::ndarray::Array2;

/// Return the row vector index corresponding to the symmetric matrix coordinates (i, j).
pub fn row_idx_from_mat_coords(n: usize, i: usize, j: usize) -> usize {
    let n = n as isize;
    let i = i as isize;
    let j = j as isize;
    if i <= j {
        (i * n - (i - 1) * i / 2 + j - i) as usize
    } else {
        (j * n - (j - 1) * j / 2 + i - j) as usize
    }
}

/// Calculate the size of a row vector given the number of leaf nodes.
fn row_vec_size_for_num_leaf_nodes(n_leaves: usize) -> usize {
    (n_leaves * (n_leaves + 1)) / 2
}

/// Create a row vector given the number leaf nodes in a symmetric disance matrix.
#[must_use]
pub fn create_row_vec_from_mat_dims(n_leaves: usize) -> Vec<f64> {
    return vec![0.0; row_vec_size_for_num_leaf_nodes(n_leaves)];
}

/// Calculate the number of leaf nodes given a row vector size.
pub fn num_leaves_from_row_vec_size(size: usize) -> usize {
    /*
    Return the number of leaves in a row vector.
     */
    (
        0.5 * ((((size * 8 + 1) as f64).sqrt()) - 1.0)
    ) as usize
}

/// Convert a row vector into a symmetric square distance matrix.
pub fn row_vec_to_symmat(row_vec: &[f64]) -> Array2<f64> {
    let num_leaf = num_leaves_from_row_vec_size(row_vec.len());
    let mut array = Array2::<f64>::zeros([num_leaf, num_leaf]);
    for i in 0..num_leaf {
        for j in 0..num_leaf {
            array[[i, j]] = row_vec[row_idx_from_mat_coords(num_leaf, i, j)];
        }
    }
    array
}


#[test]
fn test_num_leaves_from_row_vec_size() {
    assert_eq!(num_leaves_from_row_vec_size(1), 1);
    assert_eq!(num_leaves_from_row_vec_size(3), 2);
    assert_eq!(num_leaves_from_row_vec_size(6), 3);
    assert_eq!(num_leaves_from_row_vec_size(10), 4);
    assert_eq!(num_leaves_from_row_vec_size(15), 5);
    assert_eq!(num_leaves_from_row_vec_size(21), 6);
    assert_eq!(num_leaves_from_row_vec_size(28), 7);
    assert_eq!(num_leaves_from_row_vec_size(36), 8);
    assert_eq!(num_leaves_from_row_vec_size(45), 9);
    assert_eq!(num_leaves_from_row_vec_size(55), 10);
}


#[test]
fn test_row_idx_from_mat_coords() {
    assert_eq!(row_idx_from_mat_coords(3, 0, 0), 0);
    assert_eq!(row_idx_from_mat_coords(3, 0, 1), 1);
    assert_eq!(row_idx_from_mat_coords(3, 0, 2), 2);
    assert_eq!(row_idx_from_mat_coords(3, 1, 0), 1);
    assert_eq!(row_idx_from_mat_coords(3, 1, 1), 3);
    assert_eq!(row_idx_from_mat_coords(3, 1, 2), 4);
    assert_eq!(row_idx_from_mat_coords(3, 2, 0), 2);
    assert_eq!(row_idx_from_mat_coords(3, 2, 1), 4);
    assert_eq!(row_idx_from_mat_coords(3, 2, 2), 5);
}

#[test]
fn test_row_vec_to_symmat() {
    let row_vec = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0];
    let array = row_vec_to_symmat(&row_vec);
    assert_eq!(array[[0, 0]], 0.0);
    assert_eq!(array[[0, 1]], 1.0);
    assert_eq!(array[[0, 2]], 2.0);
    assert_eq!(array[[1, 0]], 1.0);
    assert_eq!(array[[1, 1]], 3.0);
    assert_eq!(array[[1, 2]], 4.0);
    assert_eq!(array[[2, 0]], 2.0);
    assert_eq!(array[[2, 1]], 4.0);
    assert_eq!(array[[2, 2]], 5.0);
}

#[test]
fn test_integration() {
    for n_leaf in 1..=10 {
        let mut array = Array2::<f64>::zeros([n_leaf, n_leaf]);
        let mut row_vec = create_row_vec_from_mat_dims(n_leaf);
        let mut value = 0.0;
        for i in 0..n_leaf {
            for j in 0..=i {
                array[[i, j]] = value;
                array[[j, i]] = value;
                row_vec[row_idx_from_mat_coords(n_leaf, i, j)] = value;
                value += 1.0;
            }
        }
        let mat = row_vec_to_symmat(&row_vec);
        assert_eq!(mat, array);
    }
}

