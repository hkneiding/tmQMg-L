import os

import pandas as pd
from tqdm import tqdm


optimization_directory = './opt+sp/'
if not os.path.exists(optimization_directory):
    os.mkdir(optimization_directory)

# load ligand fingerprints
ligands_fingerprints = pd.read_csv('../1-extraction/ligands_fingerprints.csv', sep=';')
# load stable ligands
stable_ligands = pd.read_csv('../3-stable/stable.csv', sep=';')

ligands_data = ligands_fingerprints.merge(stable_ligands, on='name', how='inner')

# load xyzs
with open('../1-extraction/ligands_xyzs.xyz', 'r') as fh:
    xyzs = {xyz.split('\n')[1]: xyz for xyz in fh.read().split('\n\n')}

# load template
with open('./gaussian-opt+sp-template.com', 'r') as fh:
    sp_template = fh.read()

for ligand in ligands_data.to_dict(orient="records"):

    subgraph_name = ligand['stable_occurrence']
    charge = ligand['charge']
    xyz = '\n'.join(xyzs[subgraph_name].split('\n')[2:])
    with open(optimization_directory + 'opt+sp-' +  subgraph_name + '.com', 'w') as fh:
        fh.write(sp_template.replace('<name>', 'opt-' + subgraph_name).replace('<charge>', str(charge)).replace('<xyz>', xyz))






