from math import sqrt, cos, sin, pi
import numpy as np

def rotation(mtrx, deg, ax):
    """Returns new coordinates after rotation on
     n degrees around ax axis

     Input:
     mtrx - matrix of coordinates
     deg - degree of angle
     ax - the axis around which the figure is rotated"""

    def matrix(n, axis):
        def cos_deg(x):
            return round(cos(pi * (x / 180)), 5)

        def sin_deg(x):
            return round(sin(pi * (x / 180)), 5)

        if axis == 'x':
            M = [[1, 0, 0],
                 [0, cos_deg(n), -sin_deg(n)],
                 [0, sin_deg(n), cos_deg(n)]]
            return M
        elif axis == 'y':
            M = [[cos_deg(n), 0, sin_deg(n)],
                 [0, 1, 0],
                 [-sin_deg(n), 0, cos_deg(n)]]
            return M
        elif axis == 'z':
            M = [[cos_deg(n), -sin_deg(n), 0],
                 [sin_deg(n), cos_deg(n), 0],
                 [0, 0, 1]]
            return M

    return np.array([np.dot(matrix(deg, ax), vect) for vect in mtrx])

def reflection(mtrx, ax):
    """Function returns coordinates after reflection about the axis

    Input:
    mtrx - matrix of coordinates
    ax - the axis around which figure is reflected
    """
    if ax == 'x':
        M = np.diag([-1, 1, 1])
    if ax == 'y':
        M = np.diag([1, -1, 1])
    if ax == 'z':
        M = np.diag([1, 1, -1])
    return np.array([np.dot(M, vect) for vect in mtrx])

def translation(mtrx, add):
    """Function returns coordinates after translation about the vector

    Input:
    mtrx - matrix of coordinates
    add - vector that adder to every node of pharmacophore
    """
    return np.array([vect + add for vect in mtrx])

