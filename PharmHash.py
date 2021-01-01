import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import sqrt
import random as rnd


def get_coords(pharmacophore):
    """
    Get coordinates from pharmacophore
    :param pharmacophore:
    :return: list of nodes and array of coords
    """
    nodes = [i[0] for i in pharmacophore]
    coords = [list(i[1]) for i in pharmacophore]
    return nodes, np.array(coords)


def center_mass(coords, ns):
    """Return center of mass considering assigned weights

    Input(coords, ns):
    'coords' - matrix of cartesian coordinates of centers.
    'ns' - nodes
    """
    # Assign 'mass' to every type of active centers: Negative, Positive, Donor, Acceptor, Hydrophobic and aRomatic
    w = {'N': sqrt(2), 'P': sqrt(3), 'D': sqrt(5), 'A': sqrt(7), 'H': sqrt(11), 'a': sqrt(13)}
    weights_sum = sum([w[n] for n in ns])

    # Make matrix of weighted coordinates by multiplication (weight/weights_sum) of nodes on their coordinates
    coords_w = [(w[n] / weights_sum) * coords[i] for i, n in enumerate(ns)]

    # And return sum by columns
    return np.array(coords_w).sum(axis=0)


def new_coordinates(coords, center):
    """Return new coordinates after moving center in origin

    Input(coords, cent):
    'coords' - matrix of coordinates of pharmacophore centers
    'cent' - vector of coordinates of center of mass
    """
    # from every node subtract coordinates of center of mass
    return np.array([c - center for c in coords])


def new_center():
    """Returns new center that is in origin"""
    return [0, 0, 0]


def inertia_tensor(coords, ns):
    """Returns tensor of inertia of molecule

    Input(coords, ns):
    'coords' - coordinates of pharmacophore centers after centering
    'ns' - nodes
    """
    # Assign 'mass' to every type of active centers: Negative, Positive, Donor, Acceptor, Hydrophobic and aRomatic
    weights = {'N': sqrt(2), 'P': sqrt(3), 'D': sqrt(5), 'A': sqrt(7), 'H': sqrt(11), 'a': sqrt(13)}

    # Calculate principal axes of inertia tensor
    Ixx = sum([weights[ns[i]] * (coords[i][1] ** 2 + coords[i][2] ** 2) for i in range(len(coords))])
    Iyy = sum([weights[ns[i]] * (coords[i][0] ** 2 + coords[i][2] ** 2) for i in range(len(coords))])
    Izz = sum([weights[ns[i]] * (coords[i][0] ** 2 + coords[i][1] ** 2) for i in range(len(coords))])

    # Calculate other values
    Ixy = -sum([weights[ns[i]] * coords[i][0] * coords[i][1] for i in range(len(coords))])
    Ixz = -sum([weights[ns[i]] * coords[i][0] * coords[i][2] for i in range(len(coords))])
    Iyz = -sum([weights[ns[i]] * coords[i][1] * coords[i][2] for i in range(len(coords))])

    return np.array([[Ixx, Ixy, Ixz],
                     [Ixy, Iyy, Iyz],
                     [Ixz, Iyz, Izz]])


def eigen(mat):
    """Function returns eigenvectors and eigenvalues of matrix.
    After that, function sort this vectors and values descending


    :param: (mat)
    'mat' - matrix (inertia tensor)
    """
    w, v = np.linalg.eig(mat)
    a = sorted(zip(w, v.T), reverse=True)
    return np.array([el[1] for el in a]).T


def only_right_vect(mat):
    """
    Returns only right-handed coordinate systems
    by triple product of vectors

    Input (mat)
    'mat' - matrix of eigenvectros
    """

    def generate_numbers(N, M: int, prefix=None):
        """
        Function returns all permutation with repetiton, where N is elements,
        M is number of positions
        """
        if M == 0:
            # write down every set of numbers in 'num' and return None
            nums.append(prefix[:])
            return
        prefix = prefix or []
        for digit in N:
            prefix.append(digit)
            generate_numbers(N, M - 1, prefix)
            prefix.pop()

    nums = []

    generate_numbers([1, -1], 3)

    Rs = []
    for i in nums:
        # multiply by numbers from 'num'
        b = np.dot(mat, np.diag(i))
        # and calculate triple products = e_x dot (e_y cross e_z)
        V = np.dot(b[:, 0], np.cross(b[:, 1], b[:, 2]))
        if V > 0:
            Rs.append(b)
    return Rs


def arounded(mat, m=1):
    """Enter vector or matrix and get rounded elements

    Input: (mat, m=1)
    'mat' - matrix of vector
    'm' - rounding accuracy
    """

    def rounded(n, m):
        """
        Returns rounded number 'n' accurate to 'm'.
        """
        a = n - n % m
        return a if n % m < m / 2 else a + m

    if type(mat[0]) in [list, np.ndarray]:
        s = [[rounded(x, m) for x in i] for i in mat]
    else:
        s = [rounded(x, m) for x in mat]
    return np.array(s)


def coordinates_2(coo, Rs, m=1):
    """Returns rounded coordinates after multiplication on eigenvector matrix

    Input(coo, Rs)
    'coo' - coordinates
    'Rs' - list of matrices of right handed eigenvectors
    """
    return [arounded(coo @ R, m) for R in Rs]


def final_coordinates(coordinate, ns):
    """
    Returns tuple of tuples with names of nodes
    :coordinate
    :ns
    """
    weights = {'N': sqrt(2), 'P': sqrt(3), 'D': sqrt(5), 'A': sqrt(7), 'H': sqrt(11), 'a': sqrt(13)}
    tuples = [[tuple(list(ns[i]) + [tuple(j)]) for i, j in enumerate(M)] for M in coordinate]
    tuples = [sorted(M, key=lambda x: weights[x[0]], reverse=True) for M in tuples]
    return sorted(list(map(tuple, tuples)))


def show_ph(coords, cent, nodes):
    """
    Demonstration of graph in 3D.

    input(coords, cent):
    'coords' - matrix of coordinates
    'cent' - vector of center mass
    'nodes' - nodes
    """
    # Assign special color to every type of centers
    colors = {"A": "red", "a": "orange", "D": "plum", "H": "cyan", "P": "olive", "N": "lime"}

    fig = plt.figure(figsize=(30, 24))
    ax = fig.add_subplot(111, projection='3d')

    # Form arrays of xs, ys, and zs from coordinates

    xs = coords[:, 0]
    ys = coords[:, 1]
    zs = coords[:, 2]

    for i in range(len(xs)):
        ax.scatter(xs[i], ys[i], zs[i], c=colors[nodes[i]], marker='o', s=300)

    # Make simple, bare axis lines through space:
    ax.plot((min(xs) - 1, max(xs) + 1), (0, 0), (0, 0), 'k')
    ax.plot((0, 0), (min(ys) - 1, max(ys) + 1), (0, 0), 'k')
    ax.plot((0, 0), (0, 0), (min(zs) - 1, max(zs) + 1), 'k')

    # Nodes on plot:
    # ax.scatter(xs, ys, zs, c='r', marker='o', s=1000)

    # Center of mass:
    ax.scatter(cent[0], cent[1], cent[2], c='k', marker='o', s=500)

    # Labels on centers:
    for i in range(len(xs)):
        ax.text(xs[i]+rnd.random()*0.1, ys[i]+rnd.random()*0.1, zs[i]+rnd.random()*0.1, nodes[i], fontsize=20)

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    plt.show()


def make_hash(pharmacophore, m=1):
    """Function returns tuple of coordinates and mark of every node in pharmacophore

    Input:
    coords - coordinates of nodes
    nodes - list of nodes
    m - accuracy of rounding
    """

    nodes, coords = get_coords(pharmacophore)
    cnt = center_mass(coords, nodes)
    new_c = new_coordinates(coords, cnt)
    I = inertia_tensor(new_c, nodes)
    eig = eigen(I)
    RS = only_right_vect(eig)
    C2 = coordinates_2(new_c, RS, m)
    CF = final_coordinates(C2, nodes)

    return CF

def view_ph(pharmacophore):
    """
    Shows pharmacophore
    :param pharmacophore:
    :return: nothing
    """
    nodes, coords = get_coords(pharmacophore)
    c = center_mass(coords, nodes)
    show_ph(coords, c, nodes)