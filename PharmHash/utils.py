from typing import List, NamedTuple, Tuple
import numpy as np
from constants import WEIGHTS, PERMUTATIONS


class PharmCenter(NamedTuple):
    node: str
    coords: np.ndarray

    @property
    def tpl(self) -> Tuple[str, Tuple]:
        return self.node, tuple(self.coords)

    @property
    def weight(self) -> float:
        return WEIGHTS[self.node]


def rounded(arr: np.ndarray, accuracy: float) -> np.ndarray:
    """
    Enter vector or matrix and get rounded elements accurate to 'm'.
    If less than half of rounding step, then round to little number
    else to bigger number
    """
    if  accuracy <= 0:
        raise ValueError('Accuracy must be positive')
    cond: np.ndarray[bool] = (arr % accuracy >= accuracy / 2)
    return arr - arr % accuracy + accuracy * (cond)


def get_inertia_tensor(nodes: List[str], coords: np.ndarray) -> np.ndarray:
    Ixx: float = sum([WEIGHTS[nodes[i]] * (coords[i][1] ** 2 + coords[i][2] ** 2) for i in range(len(coords))])
    Iyy: float = sum([WEIGHTS[nodes[i]] * (coords[i][0] ** 2 + coords[i][2] ** 2) for i in range(len(coords))])
    Izz: float = sum([WEIGHTS[nodes[i]] * (coords[i][0] ** 2 + coords[i][1] ** 2) for i in range(len(coords))])
    # Calculate other values
    Ixy: float = -sum([WEIGHTS[nodes[i]] * coords[i][0] * coords[i][1] for i in range(len(coords))])
    Ixz: float = -sum([WEIGHTS[nodes[i]] * coords[i][0] * coords[i][2] for i in range(len(coords))])
    Iyz: float = -sum([WEIGHTS[nodes[i]] * coords[i][1] * coords[i][2] for i in range(len(coords))])
    # Calculate inertia tensor
    inertia_tensor = np.array([
        [Ixx, Ixy, Ixz],
        [Ixy, Iyy, Iyz],
        [Ixz, Iyz, Izz],
    ])
    return inertia_tensor


def get_orientations(matrix: np.ndarray) -> List[np.ndarray]:
    """
    1. Calculate eigenvectors and eigenvalues of inertia tensor and sort it
       thus max eigenvalue and its eigenvector became first and min - last
       and then safe only eigenvectors.

    2. Multiply by diagonal matrix from every permutation
       and calculate determinant of matrix.
       only right-handed coordinate systems must be saved
    """
    w, v = np.linalg.eig(matrix)
    sorted_vectors = sorted(zip(w, v.T), reverse=True)
    eigenvectors: np.ndarray = np.array([el[1] for el in sorted_vectors]).T

    orientations: List[np.ndarray] = []
    for permutation in PERMUTATIONS:
        b: np.ndarray = eigenvectors @ np.diag(permutation)
        if np.linalg.det(b) > 0: 
            orientations.append(b)
    return orientations

def get_center_mass(nodes: List[str], coords: np.ndarray) -> np.ndarray:
    weights_sum = sum(WEIGHTS[node] for node in nodes)
    # Make matrix of weighted coordinates by multiplication (weight/weights_sum) of nodes on their coordinates
    weighted_coords = np.array([(WEIGHTS[node] / weights_sum) * coords[i] for i, node in enumerate(nodes)])
    center: np.ndarray = weighted_coords.sum(axis=0)        # and return sum by columns
    return center

def get_hash(accuracy: float, coords, orientations, nodes):
    # multiply every coordinates by position matrix to turn it
    # and round with given accuracy
    rounded_coordinates: List[np.ndarray] = [rounded(coords @ orientation, accuracy) for orientation in orientations]

    # make tuples that looks like (('N', (1, 0, 3)), ...)
    oriented_pharmacophores = [
        [PharmCenter(node, coord) for node, coord in zip(nodes, M)]
        for M in rounded_coordinates
    ]

    # sort tuples by weights of nodes
    tuples = [sorted(phash, key=lambda x: x.weight, reverse=True) for phash in oriented_pharmacophores]
    # return maxs
    m = [tuple(pc.tpl for pc in i) for i in tuples]
    
    maximal = max(m)
    return maximal




if __name__ == "__main__":
    a = PharmCenter('d', np.array([1,2,3]))
    b = PharmCenter('d', np.array([1,2,3]))
    c = PharmCenter('d', np.array([1,2,3]))
    d = PharmCenter('d', np.array([1,2,3]))
    tpls = a,b,c,d
    print([t.tpl for t in tpls])