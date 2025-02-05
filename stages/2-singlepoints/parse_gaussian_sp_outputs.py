import os
import re

import pandas as pd

from data_parser import DataParser


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

sp_directory = './sp/'

# load ligand library
ligands_misc = pd.read_csv('../1-extraction/ligands_misc_info.csv', sep=';')[['name', 'metal_bond_node_idx_groups']]

single_point_data = []
for ligand in ligands_misc.to_dict(orient="records"):

    for subgraph_name, metal_bond_node_idx_groups in eval(ligand['metal_bond_node_idx_groups']).items():

        file_name = 'sp-' + subgraph_name + '.out'

        if not os.path.exists(sp_directory + file_name):
            print(file_name + ': file not found.')
            continue
        else:
            print('Processing: ' + file_name)

        sp_result = DataParser(sp_directory + file_name).parse()

        if sp_result is None:
            print(file_name + ': calculation failed.')
            continue

        metal_bond_node_idx = flatten_list(metal_bond_node_idx_groups)

        metal_bound_homo_idx = None
        metal_bound_lumo_idx = None
        homo_idx = 0

        # get orbital information
        for j, _ in enumerate(sp_result['molecular_orbital_data']):

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
                    if sp_result['molecular_orbital_data'][j]['energy'] > sp_result['molecular_orbital_data'][metal_bound_homo_idx]['energy']:
                        metal_bound_homo_idx = j

            elif _['type'] == 'vir':

                # if metal bound LUMO is not set, set at first occurrence
                if metal_bound_lumo_idx is None:
                    metal_bound_lumo_idx = j
                # otherwise check if energy is lower than current metal bound LUMO
                else:
                    if sp_result['molecular_orbital_data'][j]['energy'] < sp_result['molecular_orbital_data'][metal_bound_lumo_idx]['energy']:
                        metal_bound_lumo_idx = j

            else:
                print('Error: orbital type not recognized.')
                exit()


        # calculate HOMO-LUMO
        homo_lumo_gap = sp_result['molecular_orbital_data'][homo_idx + 1]['energy'] - sp_result['molecular_orbital_data'][homo_idx]['energy']

        # get metal bound HOMO/LUMO symmetries
        metal_bound_homo_symmetries = get_orbital_symmetry(sp_result['molecular_orbital_data'][metal_bound_homo_idx]['occupations'])
        metal_bound_lumo_symmetries = get_orbital_symmetry(sp_result['molecular_orbital_data'][metal_bound_lumo_idx]['occupations'])

        # append dict
        property_dict = {
            'occurrence': subgraph_name,
            'energy': sp_result['electronic_energy'],
            'dipole_moment': sp_result['dipole_moment'],
            'homo_lumo_gap': homo_lumo_gap,
            'metal_bound_homo_energy': sp_result['molecular_orbital_data'][metal_bound_homo_idx]['energy'],
            'metal_bound_homo_s': metal_bound_homo_symmetries[0],
            'metal_bound_homo_p': metal_bound_homo_symmetries[1],
            'metal_bound_homo_d': metal_bound_homo_symmetries[2],
            'metal_bound_homo_f': metal_bound_homo_symmetries[3],
            'metal_bound_lumo_energy': sp_result['molecular_orbital_data'][metal_bound_lumo_idx]['energy'],
            'metal_bound_lumo_s': metal_bound_lumo_symmetries[0],
            'metal_bound_lumo_p': metal_bound_lumo_symmetries[1],
            'metal_bound_lumo_d': metal_bound_lumo_symmetries[2],
            'metal_bound_lumo_f': metal_bound_lumo_symmetries[3]
        }

        single_point_data.append(property_dict)

pd.DataFrame(single_point_data).to_csv('sp_summary.csv', index=False, sep=';')
