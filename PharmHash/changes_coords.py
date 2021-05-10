from math import sqrt, cos, sin, pi
import numpy as np


def axis_validator(f):
    def inner(*args):
        if not args[1] in ('x', 'y', 'z'):
            raise ValueError('Axis must be x, y or z')
        return f(*args)
    return inner


@axis_validator
def rotation(coords, ax, deg):
    """Returns new coordinates after rotation on n degrees around ax axis
     
    :param mtrx: - matrix of coordinates
    :param deg: - degree of angle
    :param ax: - the axis around which the figure is rotated
    """

    def matrix(n, axis):
        def cos_deg(x):
            return round(cos(pi * (x/180)), 5)
    
        def sin_deg(x):
            return round(sin(pi * (x/180)), 5)
        
        if axis == 'x':
            M = [[1, 0, 0],
                 [0, cos_deg(n), -sin_deg(n)],
                 [0, sin_deg(n), cos_deg(n)]]

        elif axis == 'y':
            M = [[cos_deg(n), 0, sin_deg(n)],
                 [0, 1, 0],
                 [-sin_deg(n), 0, cos_deg(n)]]

        else:
            M = [[cos_deg(n), -sin_deg(n), 0],
                 [sin_deg(n), cos_deg(n), 0],
                 [0, 0, 1]]

        return M

    return np.array([np.dot(matrix(deg, ax), vect) for vect in coords])


@axis_validator
def reflection(coords, ax):
    """Function returns coordinates after reflection about the axis

    :param coords: - matrix of coordinates
    :param ax: - the axis around which figure is reflected
    """
    if ax == 'x':
        M = np.diag([-1, 1, 1])
    elif ax == 'y':
        M = np.diag([1, -1, 1])
    else:
        M = np.diag([1, 1, -1])

    return np.array([np.dot(M, vect) for vect in coords])

def translation(coords, add):
    """Function returns coordinates after translation about the vector

    :param coords: - matrix of coordinates
    :param add: - vector that adder to every node of pharmacophore
    """
    return np.array([vect + add for vect in coords])

