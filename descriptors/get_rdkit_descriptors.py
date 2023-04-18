import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


# The RDKit derived properties do not depend on the specific
# coordinates of different instances of the same ligand. All
# occurrences of the same ligand should have the same properties.
# Therefore, simply the corresponding SMILES string can be used
# for each ligand.

# set print options
np.set_printoptions(threshold=np.inf)

# read data
ligand_data = pd.read_csv('/home/hkneiding/Documents/UiO/Data/tmQMg-L/ligands_misc_info.csv', sep=';')
# load xyz
xyzs = {}
with open('/home/hkneiding/Documents/UiO/Data/tmQMg-L/ligands_xyz.xyz', 'r') as fh:
    for xyz in fh.read().split('\n\n'):
        xyzs[xyz.split('\n')[1]] = xyz

# xyzs to write
xyz_names = []

# result variable
rdkit_descriptor_list = []

# go through ligands
for i, ligand in enumerate(tqdm(ligand_data.to_dict(orient='records'))):

    print('')
    print(list(eval(ligand['parent_metal_occurrences']).values()))
    xyz_names.append(list(eval(ligand['parent_metal_occurrences']).values())[0][0])

    # read smiles
    smi = ligand['smiles']

    # read molecule to RDKit
    mol = Chem.MolFromSmiles(smi, sanitize=False)

    # custom sanitisation
    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(
        mol,
        Chem.SanitizeFlags.SANITIZE_SETAROMATICITY|
        Chem.SanitizeFlags.SANITIZE_SETCONJUGATION|
        Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION|
        Chem.SanitizeFlags.SANITIZE_CLEANUP
    )

    # morgan fingerprint
    morgan_fingerprint = np.array(rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius=3, nBits=4096))
    # logp
    logp = rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
    # n aliphatic rings
    n_aliphatic_rings = rdMolDescriptors.CalcNumAliphaticRings(mol)
    # n aromatic rings
    n_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    # n saturated rings
    n_saturated_rings = rdMolDescriptors.CalcNumSaturatedRings(mol)
    # n rotatable bonds
    n_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)

    rdkit_descriptor_list.append({
        'ligand_name': ligand['name'],
        'smiles': Chem.MolToSmiles(mol, allHsExplicit=True),
        'morgan_fingerprint': morgan_fingerprint,
        'logp': logp,
        'n_aliphatic_rings': n_aliphatic_rings,
        'n_aromatic_rings': n_aromatic_rings,
        'n_saturated_rings': n_saturated_rings,
        'n_rotatable_bonds': n_rotatable_bonds
    })

rdkit_descriptor_df = pd.DataFrame(rdkit_descriptor_list)
rdkit_descriptor_df.to_csv('./rdkit_descriptors.csv', sep=';', index=False)
