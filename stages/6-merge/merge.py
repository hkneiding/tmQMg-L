import os
import pandas as pd


# load dfs
ligands_misc = pd.read_csv('../1-extraction/ligands_misc_info.csv', sep=';')
ligands_fingerprints = pd.read_csv('../1-extraction/ligands_fingerprints.csv', sep=';')
ligands_stable = pd.read_csv('../3-stable/stable.csv', sep=';')
rdkit_descriptors = pd.read_csv('../5-descriptors/rdkit_descriptors.csv', sep=';')
steric_descriptors = pd.read_csv('../5-descriptors/steric_descriptors.csv', sep=';')
electronic_descriptors = pd.read_csv('../5-descriptors/electronic_descriptors.csv', sep=';')

# load xyzs
with open('../1-extraction/ligands_xyzs.xyz', 'r') as fh:
    xyzs = {xyz.split('\n')[1]: xyz for xyz in fh.read().split('\n\n')}

with open('../4-optimizations/ligands_stable_opt_xyzs.xyz', 'r') as fh:
    stable_opt_xyzs = {xyz.split('\n')[1]: xyz for xyz in fh.read().split('\n\n')}

# merge descriptor CSVs
full_descriptors = pd.merge(
                    rdkit_descriptors,
                    steric_descriptors,
                    on='name',
                    how='inner'
)
full_descriptors = pd.merge(
                    full_descriptors,
                    electronic_descriptors,
                    on='name',
                    how='inner'
)

# obtain selected ligand names
names = full_descriptors['name'].to_numpy().tolist()

# subselect ligands CSVs
ligands_misc = ligands_misc[ligands_misc['name'].isin(names)]
ligands_fingerprints = ligands_fingerprints[ligands_fingerprints['name'].isin(names)]

# add stable occurrence to ligands_misc
ligands_misc = pd.merge(
                ligands_misc,
                ligands_stable,
                on='name',
                how='inner'
)

# subselect xyzs
subselected_xyzs = []
subselected_stable_xyzs = []
subselected_stable_opt_xyzs = []

for ligand in ligands_misc.to_dict(orient='records'):
    occurrences = list(eval(ligand['metal_bond_node_idx_groups']).keys())
    for occurrence in occurrences:
        subselected_xyzs.append(xyzs[occurrence])

    stable_occurrence = ligand['stable_occurrence']
    subselected_stable_xyzs.append(xyzs[stable_occurrence])
    subselected_stable_opt_xyzs.append(stable_opt_xyzs[stable_occurrence])

# write data
target_directory = '../../'
ligands_misc.to_csv(target_directory + 'ligands_misc_info.csv', index=False, sep=';')
ligands_fingerprints.to_csv(target_directory + 'ligands_fingerprints.csv', index=False, sep=';')
full_descriptors.to_csv(target_directory + 'ligands_descriptors.csv', index=False, sep=';')

target_directory = '../../xyz/'
os.makedirs(target_directory, exist_ok=True)
with open(target_directory + 'ligand_xyzs.xyz', 'w') as fh:
    fh.write('\n\n'.join(subselected_xyzs))
with open(target_directory + 'ligand_stable_xyzs.xyz', 'w') as fh:
    fh.write('\n\n'.join(subselected_stable_xyzs))
with open(target_directory + 'ligand_stable_opt_xyzs.xyz', 'w') as fh:
    fh.write('\n\n'.join(subselected_stable_opt_xyzs))




