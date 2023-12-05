import os
import re
import sys
import time
import pickle
import numpy as np
import pandas as pd
from tqdm import tqdm


# This script extracts electronic descriptors based on
# DFT calculations of the respective ligands.

def flatten_list(input_list):

    """Flattens a irregular list. Embeds any sublist as individual values in main list.

    Returns:
        list[]: The flattened list.
    """

    flattened_list = []
    for element in input_list:
        if isinstance(element, list):
            flattened_list.extend(flatten_list(element))
        else:
            flattened_list.append(element)

    return flattened_list

def get_orbital_symmetry(occupations: dict):

    """Calculates the relative orbital symmetry character in terms of s, p, d contributions.

    Arguments:
        dict: A dictionary of molecular orbitals.

    Returns:
        float: The relative s orbital contribution.
        float: The relative p orbital contribution.
        float: The relative d orbital contribution.
    """

    s_occupation = 0
    p_occupation = 0
    d_occupation = 0
    f_occupation = 0

    for key in occupations.keys():

        if key.split('-')[-1] == 's':
            s_occupation += occupations[key]
        elif key.split('-')[-1] == 'p':
            p_occupation += occupations[key]
        elif key.split('-')[-1] == 'd':
            d_occupation += occupations[key]
        elif key.split('-')[-1] == 'f':
            f_occupation += occupations[key]

    occupation_sum = s_occupation + p_occupation + d_occupation + f_occupation

    return s_occupation / occupation_sum, p_occupation / occupation_sum, d_occupation / occupation_sum, f_occupation / occupation_sum

# get input argument
output_path = sys.argv[1]

# get misc data
ligand_data = pd.read_csv('./../ligands_misc_info_full.csv', sep=';')

# get most stable ligand names
with open('./../stable.csv', 'r') as fh:
    ligand_stable_names = {_.split(';')[0]: _.split(';')[1] for _ in fh.read().strip().split('\n')}

# get output result files
result_files = [f for f in os.listdir(output_path) if os.path.isfile(os.path.join(output_path, f))]

ligand_names = []

# list to store extracted data
failures = []
electronic_descriptor_list = []
for result_file in tqdm(result_files):

    with open(output_path + result_file, 'rb') as fh:
        sp_result = pickle.load(fh)

    for i in range(len(sp_result)):

        try:

            # get ligand and occurrence
            ligand_name = sp_result[i]['id'].split('/')[0]

            if ligand_name in ligand_names:
                continue

            ligand_names.append(ligand_name)

            # get occurrence name
            occurrence_name = sp_result[i]['id'].split('/')[1].replace('opt+', '').replace('sp-', '').replace('.out', '')

            # get metal bond node indices
            metal_bond_node_idx_group_occurrence_dict = eval(ligand_data.loc[ligand_data['name'] == ligand_name]['metal_bond_node_idx_groups'].item())
            metal_bond_node_idx = flatten_list(metal_bond_node_idx_group_occurrence_dict[occurrence_name])

            metal_bound_homo_idx = None
            metal_bound_lumo_idx = None
            homo_idx = 0

            # get orbital information
            for j, _ in enumerate(sp_result[i]['molecular_orbital_data']):

                orbital_atom_idx = [int(re.sub('[^0-9]', '', atom_id)) - 1 for atom_id in _['occupations'].keys()]

                # check if there is overlap between the metal bound atom idx and this orbitals atom idx
                # if not continue to next orbital
                if not bool(set(metal_bond_node_idx) & set(orbital_atom_idx)):
                    continue

                if _['type'] == 'occ':

                    # set id of HOMO
                    homo_idx = j

                    # if metal bound HOMO is not set, set at first occurrence
                    if metal_bound_homo_idx is None:
                        metal_bound_homo_idx = j
                    # otherwise check if energy is higher than current metal bound HOMO
                    else:
                        if sp_result[i]['molecular_orbital_data'][j]['energy'] > sp_result[i]['molecular_orbital_data'][metal_bound_homo_idx]['energy']:
                            metal_bound_homo_idx = j

                elif _['type'] == 'vir':

                    # if metal bound LUMO is not set, set at first occurrence
                    if metal_bound_lumo_idx is None:
                        metal_bound_lumo_idx = j
                    # otherwise check if energy is lower than current metal bound LUMO
                    else:
                        if sp_result[i]['molecular_orbital_data'][j]['energy'] < sp_result[i]['molecular_orbital_data'][metal_bound_lumo_idx]['energy']:
                            metal_bound_lumo_idx = j

                else:
                    print('Error: orbital type not recognized.')
                    exit()


            # calculate HOMO-LUMO
            homo_lumo_gap = sp_result[i]['molecular_orbital_data'][homo_idx + 1]['energy'] - sp_result[i]['molecular_orbital_data'][homo_idx]['energy']

            # get metal bound HOMO/LUMO symmetries
            metal_bound_homo_symmetries = get_orbital_symmetry(sp_result[i]['molecular_orbital_data'][metal_bound_homo_idx]['occupations'])
            metal_bound_lumo_symmetries = get_orbital_symmetry(sp_result[i]['molecular_orbital_data'][metal_bound_lumo_idx]['occupations'])

            # append dict
            property_dict = {
                'name': ligand_name,
                'occurrence_name': occurrence_name,
                'dipole_moment': sp_result[i]['dipole_moment'],
                'homo_lumo_gap': homo_lumo_gap,
                'metal_bound_homo_energy': sp_result[i]['molecular_orbital_data'][metal_bound_homo_idx]['energy'],
                'metal_bound_homo_s': metal_bound_homo_symmetries[0],
                'metal_bound_homo_p': metal_bound_homo_symmetries[1],
                'metal_bound_homo_d': metal_bound_homo_symmetries[2],
                'metal_bound_homo_f': metal_bound_homo_symmetries[3],
                'metal_bound_lumo_energy': sp_result[i]['molecular_orbital_data'][metal_bound_lumo_idx]['energy'],
                'metal_bound_lumo_s': metal_bound_lumo_symmetries[0],
                'metal_bound_lumo_p': metal_bound_lumo_symmetries[1],
                'metal_bound_lumo_d': metal_bound_lumo_symmetries[2],
                'metal_bound_lumo_f': metal_bound_lumo_symmetries[3]
            }

            # add values from optimisation
            if 'isotropic_polarisability' in list(sp_result[i].keys()):
                property_dict['polarisability'] = sp_result[i]['isotropic_polarisability']
            if 'frequencies' in list(sp_result[i].keys()):
                property_dict['largest_frequency'] = sp_result[i]['frequencies'][-1]
            if 'principal_moments' in list(sp_result[i].keys()):

                line = sp_result[i]['principal_moments'].split('--')[-1][2:]

                # add principal moment ratios if no error in Gaussian output
                if '*' not in line:

                    I1 = float(line[0:10])
                    I2 = float(line[10:20])
                    I3 = float(line[20:30])

                    property_dict['I1/I3'] = I1 / I3
                    property_dict['I2/I3'] = I2 / I3

            if 'molar_volume' in list(sp_result[i].keys()):
                property_dict['molar_volume'] = sp_result[i]['molar_volume']

            electronic_descriptor_list.append(property_dict)

        except Exception as e:
            print(e)
            print(sp_result[i]['id'])
            failures.append(sp_result[i]['id'])

print('Number of failures:', len(failures))
print('Failures:', failures)

electronic_descriptors_df = pd.DataFrame(electronic_descriptor_list)
electronic_descriptors_df.to_csv('./electronic_descriptors_' + output_path.split('_')[0] + '.csv', sep=';', index=False)
