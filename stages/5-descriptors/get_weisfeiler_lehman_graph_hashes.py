import pandas as pd
from pysmiles import read_smiles
import networkx as nx
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


df = pd.read_csv('./../1-extraction/ligands_misc_info.csv', sep=';')

hash_dict_list = []
for ligand in tqdm(df.to_dict(orient="records")):

    # read molecule to RDKit
    mol = Chem.MolFromSmiles(ligand['smiles'])

    smi = Chem.MolToSmiles(mol)

    G = read_smiles(smi, explicit_hydrogen=True, strict=False)

    base_hash = nx.weisfeiler_lehman_graph_hash(G, iterations=3)
    atom_attribution_hash = nx.weisfeiler_lehman_graph_hash(G, node_attr='element', iterations=3)
    bond_attribution_hash = nx.weisfeiler_lehman_graph_hash(G, edge_attr='order', iterations=3)
    atom_bond_attribution_hash = nx.weisfeiler_lehman_graph_hash(G, node_attr='element', edge_attr='order', iterations=3)

    hash_dict = {
        'name': ligand['name'],
        'base_hash': base_hash,
        'atom_attribution_hash': atom_attribution_hash,
        'bond_attribution_hash': bond_attribution_hash,
        'atom_bond_attribution_hash': atom_bond_attribution_hash
    }

    hash_dict_list.append(hash_dict)

hash_df = pd.DataFrame(hash_dict_list)
hash_df.to_csv('./weisfeiler_lehman_graph_hashes.csv', sep=';', index=False)
