import numpy as np
from math import sqrt



class cached_property(object):
    def __init__(self, func):
        self.func = func

    def __get__(self, instance, cls=None):
        result = instance.__dict__[self.func.__name__] = self.func(instance)
        return result


class Pharmacophore:
    """
    Class that represents pharmacophores
    You can calculate hash of pharmacophore with given accuracy
    To change accuracy you should
    """
    # Assign 'mass' to every type of active centers: Negative, Positive, Donor, Acceptor, Hydrophobic and aromatic
    w = {'N': sqrt(2), 'P': sqrt(3), 'D': sqrt(5), 'A': sqrt(7), 'H': sqrt(11), 'a': sqrt(13)}
    # All permutations with repetitions to find only right-hand vectors
    nums = ((-1, -1, -1), (-1, -1, 1), (-1, 1, -1), (-1, 1, 1), (1, -1, -1), (1, -1, 1), (1, 1, -1), (1, 1, 1))


    def __init__(self, name, ph_coords, rmsd):
        """
        :param name:
        :param ph_coords:
        :param rmsd:
        """
        self.name = name
        self.nodes = [i[0] for i in ph_coords]
        self.coords = np.array([list(i[1]) for i in ph_coords])
        self.rmsd = rmsd


    @cached_property
    def positions(self):
        """
        Return positions of pharmacophore without rounding
        """
        weights_sum = sum(Pharmacophore.w[node] for node in self.nodes)

        # Make matrix of weighted coordinates by multiplication (weight/weights_sum) of nodes on their coordinates
        coords_w = np.array([(Pharmacophore.w[node] / weights_sum) * self.coords[i] for i, node in enumerate(self.nodes)])

        center = coords_w.sum(axis=0)        # and return sum by columns

        # Move center at the origin, move all nodes 
        # From every node subtract coordinates of center of mass
        self.coords = np.array([coord - center for coord in self.coords])

        # Calculate principal axes of inertia tensor
        Ixx = sum([Pharmacophore.w[self.nodes[i]] * (self.coords[i][1] ** 2 + self.coords[i][2] ** 2) for i in range(len(self.coords))])
        Iyy = sum([Pharmacophore.w[self.nodes[i]] * (self.coords[i][0] ** 2 + self.coords[i][2] ** 2) for i in range(len(self.coords))])
        Izz = sum([Pharmacophore.w[self.nodes[i]] * (self.coords[i][0] ** 2 + self.coords[i][1] ** 2) for i in range(len(self.coords))])
        # Calculate other values
        Ixy = -sum([Pharmacophore.w[self.nodes[i]] * self.coords[i][0] * self.coords[i][1] for i in range(len(self.coords))])
        Ixz = -sum([Pharmacophore.w[self.nodes[i]] * self.coords[i][0] * self.coords[i][2] for i in range(len(self.coords))])
        Iyz = -sum([Pharmacophore.w[self.nodes[i]] * self.coords[i][1] * self.coords[i][2] for i in range(len(self.coords))])
        # Calculate inertia tensor
        inertia_tensor =  np.array([[Ixx, Ixy, Ixz],
                                    [Ixy, Iyy, Iyz],
                                    [Ixz, Iyz, Izz]])

        w, v = np.linalg.eig(inertia_tensor)        # Calculate eigenvectors and eigenvalues of inertia tensor and sort it
                                                    # thus max eigenvalue and its eigenvector became first and min - last
        eigenvectors = np.array([el[1] for el in sorted(zip(w, v.T), reverse=True)]).T  # and then safe only eigenvectors           

        positions = []
        for num in Pharmacophore.nums:
            b = np.dot(eigenvectors, np.diag(num))            # multiply by vectors from 'num'
            V = np.dot(b[:, 0], np.cross(b[:, 1], b[:, 2]))   # and calculate triple products = e_x dot (e_y cross e_z)
            if V > 0:                                         # only right-handed coordinate systems must be saved
                positions.append(b)

        return positions


    def rounded(n, m):
        """
        Enter vector or matrix and get rounded elements accurate to 'm'
        """
        if  m <= 0:
            raise ValueError
        cond = (n%m >= m/2)         # if less than half of rounding step, then round to little number
        return n - n%m + m*(cond)   # else to big number
        

    def hash(self, acc):
        """Function returns tuple of coordinates and mark of every node in pharmacophore
        :param acc: - accuracy of rounding
        """     
        # multiply every coordinates by position matrix to turn it
        # and round with 'm'-accuracy
        rounded_coordinates = [Pharmacophore.rounded(self.coords @ pos, acc) for pos in self.positions]

        # make tuples that looks like (('N', (1, 0, 3)), ...)
        tuples = [[tuple(list(self.nodes[i]) + [tuple(coord)]) for i, coord in enumerate(M)] for M in rounded_coordinates]

        # sort tuples by weights of nodes
        tuples = [sorted(M, key=lambda x: Pharmacophore.w[x[0]], reverse=True) for M in tuples]
        # return max
        return sorted(list(map(tuple, tuples)))[0]


    def is_match(self, phrm, acc):
        if not isinstance(phrm, Pharmacophore):
            raise TypeError(f'{phrm} must be Pharmacophore object')
        return phrm.hash(acc) == self.hash(acc)

