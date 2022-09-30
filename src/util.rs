use ndarray::Array2;

/// Return the row vector index corresponding to the symmetric matrix coordinates (i, j).
///
/// # Arguments
///
/// * `n`: - The number of rows/columns in the matrix.
/// * `i`: - The row index.
/// * `j`: - The column index.
///
/// # Examples
///
/// ```
/// use phylodm::util::row_idx_from_mat_coords;
/// assert_eq!(row_idx_from_mat_coords(3, 0, 0), 0);
/// ```
#[must_use]
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

/// Calculate the size of a row vector given the number of rows/columns in the matrix.
///
/// # Arguments
///
/// * `size`: - The number of rows/columns in the matrix.
///
/// # Examples
///
/// ```
/// use phylodm::util::row_vec_size_from_mat_size;
/// assert_eq!(row_vec_size_from_mat_size(3), 6);
/// ```
#[must_use]
pub fn row_vec_size_from_mat_size(size: usize) -> usize {
    (size * (size + 1)) / 2
}

#[test]
fn test_row_vec_size_from_mat_size() {
    assert_eq!(row_vec_size_from_mat_size(1), 1);
    assert_eq!(row_vec_size_from_mat_size(2), 3);
    assert_eq!(row_vec_size_from_mat_size(3), 6);
    assert_eq!(row_vec_size_from_mat_size(4), 10);
    assert_eq!(row_vec_size_from_mat_size(5), 15);
    assert_eq!(row_vec_size_from_mat_size(6), 21);
}

/// Create a row vector given the number rows/columns in a symmetric distance matrix.
///
/// # Arguments
///
/// * `size`: - The number of rows/columns in the distance matrix.
///
/// # Examples
///
/// ```
/// use phylodm::util::create_row_vec_from_mat_dims;
/// assert_eq!(create_row_vec_from_mat_dims(3), vec![0.0; 6]);
/// ```
#[must_use]
pub fn create_row_vec_from_mat_dims(size: usize) -> Vec<f64> {
    vec![0.0; row_vec_size_from_mat_size(size)]
}

#[test]
fn test_create_row_vec_from_mat_dims() {
    assert_eq!(create_row_vec_from_mat_dims(1), vec![0.0; 1]);
    assert_eq!(create_row_vec_from_mat_dims(2), vec![0.0; 3]);
    assert_eq!(create_row_vec_from_mat_dims(3), vec![0.0; 6]);
    assert_eq!(create_row_vec_from_mat_dims(4), vec![0.0; 10]);
    assert_eq!(create_row_vec_from_mat_dims(5), vec![0.0; 15]);
    assert_eq!(create_row_vec_from_mat_dims(6), vec![0.0; 21]);
}

/// Calculate the number of rows/columns in the matrix given a row vector size.
///
/// # Arguments
///
/// * `size`: - The number of elements in a row vector.
///
/// # Examples
///
/// ```
/// use phylodm::util::mat_size_from_row_vec_size;
/// assert_eq!(mat_size_from_row_vec_size(6), 3);
/// ```
#[must_use]
pub fn mat_size_from_row_vec_size(size: usize) -> usize {
    (0.5 * ((((size * 8 + 1) as f64).sqrt()) - 1.0)) as usize
}

#[test]
fn test_mat_size_from_row_vec_size() {
    assert_eq!(mat_size_from_row_vec_size(1), 1);
    assert_eq!(mat_size_from_row_vec_size(3), 2);
    assert_eq!(mat_size_from_row_vec_size(6), 3);
    assert_eq!(mat_size_from_row_vec_size(10), 4);
    assert_eq!(mat_size_from_row_vec_size(15), 5);
    assert_eq!(mat_size_from_row_vec_size(21), 6);
}

/// Convert a row vector into a symmetric distance matrix.
///
/// # Arguments
///
/// * `row_vec`: - The row vector to convert.
///
/// ```
/// use phylodm::util::row_vec_to_symmat;
/// row_vec_to_symmat(&vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0]);
/// ```
#[must_use]
pub fn row_vec_to_symmat(row_vec: &[f64]) -> Array2<f64> {
    let num_leaf = mat_size_from_row_vec_size(row_vec.len());
    let mut array = Array2::<f64>::zeros([num_leaf, num_leaf]);
    for i in 0..num_leaf {
        for j in 0..num_leaf {
            array[[i, j]] = row_vec[row_idx_from_mat_coords(num_leaf, i, j)];
        }
    }
    array
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
