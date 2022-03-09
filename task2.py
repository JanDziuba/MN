from math import hypot

import numpy as np
import scipy.linalg


def _get_givens_c_s(v_i, v_j):
    """Compute matrix entries for Givens rotation."""
    r = hypot(v_i, v_j)
    c = v_i / r
    s = -v_j / r
    return c, s


def _givens_rotate_vector(v, c, s, i, j):
    row_i_value = v[i]
    row_j_value = v[j]

    v[i] = c * row_i_value - s * row_j_value
    v[j] = s * row_i_value + c * row_j_value


def _givens_rotate_matrix(A, c, s, i, j):
    for column in range(A.shape[1]):
        row_i_value = A[i][column]
        row_j_value = A[j][column]

        A[i][column] = c * row_i_value - s * row_j_value
        A[j][column] = s * row_i_value + c * row_j_value


def _givens_rotate_matrix_last_row(A, last_row, c, s, i):
    for column in range(A.shape[1]):
        row_i_value = A[i][column]
        last_row_value = last_row[column]

        A[i][column] = c * row_i_value - s * last_row_value
        last_row[column] = s * row_i_value + c * last_row_value


class LZNK:
    def __init__(self, A):
        self.N = A.shape[0]
        self.p = A.shape[1]

        # Initialize orthogonal matrix Q and upper triangular matrix R.
        QT = []
        R = np.copy(A)

        # Iterate over lower triangular matrix.
        (rows, cols) = np.tril_indices(self.N, -1, self.p)
        for (row, col) in zip(rows, cols):

            # Compute Givens rotation matrix and
            # zero-out lower triangular matrix entries.
            if R[row][col] != 0:
                (c, s) = _get_givens_c_s(R[col][col], R[row][col])

                _givens_rotate_matrix(R, c, s, col, row)
                QT.append((c, s, col, row))

        R = np.delete(R, np.s_[R.shape[1]:], 0)
        self.QT = QT
        self.R = R

    def QT_multiply(self, b):
        b = b.copy()
        for givens in self.QT:
            _givens_rotate_vector(b, *givens)
        return b

    def solve(self, b):
        return scipy.linalg.solve_triangular(self.R, np.delete(self.QT_multiply(b), np.s_[self.R.shape[1]:], 0))

    def addrow(self, r):
        r = r.copy()
        self.N += 1

        for col in range(len(r)):
            (c, s) = _get_givens_c_s(self.R[col][col], r[col])
            _givens_rotate_matrix_last_row(self.R, r, c, s, col)
            self.QT.append((c, s, col, self.N - 1))

        return self

