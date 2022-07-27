use ndarray::Array2;
use phylodm::util::{create_row_vec_from_mat_dims, row_idx_from_mat_coords, row_vec_to_symmat};

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
