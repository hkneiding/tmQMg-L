import os

import pandas as pd
from tqdm import tqdm


single_point_directory = './sp/'
if not os.path.exists(single_point_directory):
    os.mkdir(single_point_directory)

# load ligand library
ligands_fingerprints = pd.read_csv('../1-extraction/ligands_fingerprints.csv', sep=';')
ligands_misc = pd.read_csv('../1-extraction/ligands_misc_info.csv', sep=';')[['name', 'occurrence', 'parent_metal_occurrences', 'metal_bond_node_idx_groups']]
with open('../1-extraction/ligands_xyzs.xyz', 'r') as fh:
    xyzs = {xyz.split('\n')[1]: xyz for xyz in fh.read().split('\n\n')}
# merge csvs
ligands_data = ligands_fingerprints.merge(ligands_misc, on='name', how='inner')

# load template
with open('./gaussian-sp-template.com', 'r') as fh:
    sp_template = fh.read()

for ligand in ligands_data.to_dict(orient="records"):

    charge = ligand['charge']
    for subgraph_name in eval(ligand['metal_bond_node_idx_groups']).keys():
        xyz = '\n'.join(xyzs[subgraph_name].split('\n')[2:])
        with open(single_point_directory + 'sp-' +  subgraph_name + '.com', 'w') as fh:
            fh.write(sp_template.replace('<name>', 'sp-' + subgraph_name).replace('<charge>', str(charge)).replace('<xyz>', xyz))

        # prepare only one input for single atom ligands
        if ligand['n_atoms'] == 1:
            break




