from math import sqrt
from typing import Dict, Tuple


WEIGHTS: Dict[str, float] = {
    'N': sqrt(2),
    'P': sqrt(3),
    'D': sqrt(5),
    'A': sqrt(7),
    'H': sqrt(11),
    'a': sqrt(13),
}

    # Assign 'mass' to every type of active centers: Negative, Positive, Donor, Acceptor, Hydrophobic and aromatic

    # All permutations with repetitions to find only right-hand vectors
PERMUTATIONS: Tuple[Tuple[int, int, int]] = (
    (-1, -1, -1), 
    (-1, -1, 1), 
    (-1, 1, -1), 
    (-1, 1, 1), 
    (1, -1, -1), 
    (1, -1, 1), 
    (1, 1, -1), 
    (1, 1, 1),
)
