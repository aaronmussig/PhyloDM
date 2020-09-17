from timeit import Timer

import numpy as np

from phylodm.symmat import SymMat

N_TESTS = 2500


def generate_test_matrix(n):
    expected = np.full((n, n), -1, dtype=np.int64)
    expected[np.triu_indices(n, 1)] = np.random.randint(0, 1e6, (n * (n - 1)) // 2)
    expected = expected + expected.T
    expected[np.diag_indices(n)] = np.random.randint(0, 1e6, n)
    assert (np.all(np.abs(expected - expected.T) < 1e-8))
    return expected


def generate_test_symmat(n):
    # Create a test matrix and populate it with expected data.
    expected = generate_test_matrix(n)
    mat = SymMat.get_from_indices(list(map(str, range(n))), np.dtype('int64'), -1)
    for i in range(n):
        for j in range(i + 1):
            mat.set_value(str(i), str(j), expected[i][j])
    return mat


def test_as_matrix(in_mat):
    in_mat.as_matrix()


if __name__ == '__main__':
    test_symmat = generate_test_symmat(N_TESTS)

    t = Timer(lambda: test_as_matrix(test_symmat))
    print(t.timeit(number=100))

    pass
