from typing import Any, Dict, List, Tuple


def load_from_coords(fname: str) -> Dict[str, Any]:
    """
    Reads pharmacophore from  coords file.
    :param fname: coords-file name
    :return: dictionary of coords
    """
    dic = {}

    with open(fname, 'r') as f:
        feature_coords: List[Tuple[str, Tuple[float, float, float]]] = []
        for line in f:       
            line = line.strip().split()
            if not line:   # if line has no symbols -> add pharmacophore to our dict
                dic[name] = tuple(feature_coords)
                feature_coords = []
            elif len(line) == 1:  # if line has one str -> this is name
                name = line[0]  
            else:    
                label, *coords = line
                coords: Tuple[float, float, float] = tuple(map(float, coords))
                feature_coords.append((label, coords))
    return dic


def load_rmsd(fname: str) -> Dict[str, float]:
    """
    Reads RMSD from .rms file
    :param fname: .rms file
    :returns dictionary of RMSD's:
    """
    dic = {}
    with open (fname) as f:
        for line in f:
            line  = line.strip().split()
            dic[line[0]+'_'+line[1]] = float(line[2])
    return dic

