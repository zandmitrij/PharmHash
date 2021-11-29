import os
import re

from PharmHash.reading import load_from_coords, load_rmsd
from PharmHash.changes_coords import reflection
from PharmHash.PharmHash import Pharmacophore as ph

import logging
logging.basicConfig(level=logging.ERROR, filename="logs/my_file", filemode='w')

logger = logging.getLogger()


data = os.path.join(os.getcwd(), 'data')
ph_coords = {}
RMSD = {}
for file in os.listdir(data):
    if re.search('\.coords$', file):
        my_file = os.path.join(data, file)
        ph_coords.update(load_from_coords(my_file))
    if re.search('\.rms$', file):
        my_file = os.path.join(data, file)
        RMSD.update(load_rmsd(my_file))

accs = [0.1, 0.2, 0.5, 1, 1.5, 2]

with open('data\log_file3.csv', 'w') as f:
    a = ','.join(map(str, accs))
    f.write(f'name,{a},RMSD\n')

    for name, coord in ph_coords.items():
        str_to_write = f'{name},'

        ph1 = ph(name, coord, RMSD[name])
        ph2 = ph(name, coord, RMSD[name])
        ph2.coords = reflection(ph2.coords, 'x')

        for acc in accs:
            try:
                ss = f'{ph1.is_match(ph2, acc)},'
            except:
                logger.exception(f'\nError was happened at {name}, {acc}')
                ss = f'{False},'
            finally:
                str_to_write += ss
        str_to_write += f'{ph1.rmsd}\n'
        f.write(str_to_write)

