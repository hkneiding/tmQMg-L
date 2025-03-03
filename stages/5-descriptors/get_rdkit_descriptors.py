import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


# This script calculates a selection of RDKit descriptors for
# each ligand. The RDKit derived properties do not depend
# on the specific coordinates of different instances of the
# same ligand. All occurrences of the same ligand should have
# the same properties. Therefore, simply the corresponding
# SMILES string is used for each ligand.

# set print options
np.set_printoptions(threshold=np.inf)

# read data
ligand_data = pd.read_csv('./../1-extraction/ligands_misc_info.csv', sep=';')
# load xyz
xyzs = {}
with open('./../1-extraction/ligands_xyzs.xyz', 'r') as fh:
    for xyz in fh.read().split('\n\n'):
        xyzs[xyz.split('\n')[1]] = xyz

# result variables
morgan_fingerprint_list = []
rdkit_descriptor_list = []

# go through ligands
for i, ligand in enumerate(tqdm(ligand_data.to_dict(orient='records'))):

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

    morgan_fingerprint_list.append({
        'name': ligand['name'],
        'morgan_fingerprint': str(morgan_fingerprint).replace('[','').replace(']','')
    })

    rdkit_descriptor_list.append({
        'name': ligand['name'],
        'logp': logp,
        'n_aliphatic_rings': n_aliphatic_rings,
        'n_aromatic_rings': n_aromatic_rings,
        'n_saturated_rings': n_saturated_rings,
        'n_rotatable_bonds': n_rotatable_bonds
    })

morgan_fingerprint_df = pd.DataFrame(morgan_fingerprint_list)
# morgan_fingerprint_df.to_csv('morgan_fingerprint.csv', sep=';', index=False)

rdkit_descriptor_df = pd.DataFrame(rdkit_descriptor_list)
rdkit_descriptor_df.to_csv('./rdkit_descriptors.csv', sep=';', index=False)
