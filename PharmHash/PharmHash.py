from __future__ import annotations

from typing import Callable, List, Tuple

import numpy as np

from .utils import get_center_mass, get_hash, get_inertia_tensor, get_orientations


class cached_property(object):
    def __init__(self, func: Callable) -> None:
        self.func: Callable = func

    def __get__(self, instance, cls=None):
        result = instance.__dict__[self.func.__name__] = self.func(instance)
        return result


class Pharmacophore:
    """
    Class that represents pharmacophores
    You can calculate hash of pharmacophore with given accuracy
    To change accuracy you should
    """

    def __init__(self, name: str, ph_coords, rmsd: float):
        """
        :param name:
        :param ph_coords:
        :param rmsd:
        """
        self.name: str = name
        self.nodes: List[str] = [i[0] for i in ph_coords]
        self.coords: np.ndarray = np.array([list(i[1]) for i in ph_coords])
        self.rmsd = rmsd

    @cached_property
    def orientations(self) -> List[np.ndarray]:
        """
        Return all possible positions of pharmacophore without rounding.

        Move center at the origin, move all nodes 
        From every node subtract coordinates of center of mass
        """
        mass_center = get_center_mass(self.nodes, self.coords)
        self.coords: np.ndarray = np.array([coord - mass_center for coord in self.coords])
        inertia_tensor = get_inertia_tensor(self.nodes, self.coords)
        orientations = get_orientations(inertia_tensor)
        return orientations

    def hash(self, accuracy: float) -> Tuple:
        """
        Function returns tuple of coordinates and mark of every node in pharmacophore
        :param acc: - accuracy of rounding
        """
        return get_hash(accuracy, self.coords, self.orientations, self.nodes)

    def is_match(self, phrm: Pharmacophore, accuracy: float) -> bool:
        if not isinstance(phrm, Pharmacophore):
            raise TypeError(f'{phrm} must be Pharmacophore object')
        return phrm.hash(accuracy) == self.hash(accuracy)
    
