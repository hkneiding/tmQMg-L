import os
import time
import pickle
import pandas as pd
from tqdm import tqdm

def flatten_list(input_list):

    """Flattens a irregular list. Embeds any sublist as individual values in main list.

    Returns:
        list[]: The flattend list.
    """

    flattend_list = []
    for element in input_list:
        if isinstance(element, list):
            flattend_list.extend(flatten_list(element))
        else:
            flattend_list.append(element)

    return flattend_list

# get misc data
ligand_data = pd.read_csv('/home/hkneiding/Documents/UiO/Data/tmQMg-L/ligands_misc_info.csv', sep=';')

# get most stable ligand names
with open('./stable.txt', 'r') as fh:
    ligand_stable_names = {_.split(' | ')[0]: _.split(' | ')[1] for _ in fh.read().strip().split('\n')}

# get SP output result files
sp_output_path = './sp_outputs/'
result_files = [f for f in os.listdir(sp_output_path) if os.path.isfile(os.path.join(sp_output_path, f))]

# list to store extracted data
extracted_data = []
for result_file in tqdm(result_files):

    with open(sp_output_path + result_file, 'rb') as fh:
        sp_result = pickle.load(fh)

    for i in range(len(sp_result)):

        # get ligand and occurrence
        ligand_name = sp_result[i]['id'].split('/')[0]
        occurrence_name = sp_result[i]['id'].split('/')[1].replace('sp-', '').replace('.out', '')

        # get metal bond node indices
        metal_bond_node_idx_group_occurrence_dict = eval(ligand_data.loc[ligand_data['name'] == ligand_name]['metal_bond_node_idx_groups'].item())
        metal_bond_node_idx = flatten_list(metal_bond_node_idx_group_occurrence_dict[occurrence_name])


        print(metal_bond_node_idx)
        # get orbital information
        for _ in sp_result[i]['molecular_orbital_data']:
            print(_)
            exit()
        # append dict
        extracted_data.append({
            'name': ligand_name,
            'occurrence_name': occurrence_name,
            'dipole_moment': sp_result[i]['dipole_moment'],
        })

        print(extracted_data)
        exit()
