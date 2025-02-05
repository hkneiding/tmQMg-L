import os
import re

import pandas as pd

from data_parser import DataParser


element_identifiers = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
                       'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K',
                       'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',
                       'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
                       'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd',
                       'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs',
                       'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
                       'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta',
                       'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb',
                       'Bi', 'Po', 'At', 'Rn']

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

opt_sp_directory = './opt+sp/'
# load ligand library
ligands_misc = pd.read_csv('../1-extraction/ligands_misc_info.csv', sep=';')[['name', 'metal_bond_node_idx_groups']]

opt_sp_xyzs = []
opt_sp_data = []
for ligand in ligands_misc.to_dict(orient="records"):

    for subgraph_name, metal_bond_node_idx_groups in eval(ligand['metal_bond_node_idx_groups']).items():

        file_name = 'opt+sp-' + subgraph_name + '.out'

        print(file_name)

        if not os.path.exists(opt_sp_directory + file_name):
            print(file_name + ': file not found.')
            continue

        opt_sp_result = DataParser(opt_sp_directory + file_name).parse()

        if opt_sp_result is None:
            print(file_name + ': calculation failed.')
            continue


        metal_bond_node_idx = flatten_list(metal_bond_node_idx_groups)

        metal_bound_homo_idx = None
        metal_bound_lumo_idx = None
        homo_idx = 0

        # get orbital information
        for j, _ in enumerate(opt_sp_result['molecular_orbital_data']):

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
                    if opt_sp_result['molecular_orbital_data'][j]['energy'] > opt_sp_result['molecular_orbital_data'][metal_bound_homo_idx]['energy']:
                        metal_bound_homo_idx = j

            elif _['type'] == 'vir':

                # if metal bound LUMO is not set, set at first occurrence
                if metal_bound_lumo_idx is None:
                    metal_bound_lumo_idx = j
                # otherwise check if energy is lower than current metal bound LUMO
                else:
                    if opt_sp_result['molecular_orbital_data'][j]['energy'] < opt_sp_result['molecular_orbital_data'][metal_bound_lumo_idx]['energy']:
                        metal_bound_lumo_idx = j

            else:
                print('Error: orbital type not recognized.')
                exit()


        # calculate HOMO-LUMO
        homo_lumo_gap = opt_sp_result['molecular_orbital_data'][homo_idx + 1]['energy'] - opt_sp_result['molecular_orbital_data'][homo_idx]['energy']

        # get metal bound HOMO/LUMO symmetries
        metal_bound_homo_symmetries = get_orbital_symmetry(opt_sp_result['molecular_orbital_data'][metal_bound_homo_idx]['occupations'])
        metal_bound_lumo_symmetries = get_orbital_symmetry(opt_sp_result['molecular_orbital_data'][metal_bound_lumo_idx]['occupations'])

        # append dict
        property_dict = {
            'occurrence': subgraph_name,
            'energy': opt_sp_result['electronic_energy'],
            'dipole_moment': opt_sp_result['dipole_moment'],
            'homo_lumo_gap': homo_lumo_gap,
            'metal_bound_homo_energy': opt_sp_result['molecular_orbital_data'][metal_bound_homo_idx]['energy'],
            'metal_bound_homo_s': metal_bound_homo_symmetries[0],
            'metal_bound_homo_p': metal_bound_homo_symmetries[1],
            'metal_bound_homo_d': metal_bound_homo_symmetries[2],
            'metal_bound_homo_f': metal_bound_homo_symmetries[3],
            'metal_bound_lumo_energy': opt_sp_result['molecular_orbital_data'][metal_bound_lumo_idx]['energy'],
            'metal_bound_lumo_s': metal_bound_lumo_symmetries[0],
            'metal_bound_lumo_p': metal_bound_lumo_symmetries[1],
            'metal_bound_lumo_d': metal_bound_lumo_symmetries[2],
            'metal_bound_lumo_f': metal_bound_lumo_symmetries[3]
        }

        property_dict['polarisability'] = opt_sp_result['isotropic_polarisability']

        if 'frequencies' in opt_sp_result.keys():
            property_dict['largest_frequency'] = opt_sp_result['frequencies'][-1]

        # add principal moment ratios if no error in Gaussian output
        if 'principal_moments' in opt_sp_result.keys():
            line = opt_sp_result['principal_moments'].split('--')[-1][2:]
            if '*' not in line:

                I1 = float(line[0:10])
                I2 = float(line[10:20])
                I3 = float(line[20:30])

                property_dict['I1/I3'] = I1 / I3
                property_dict['I2/I3'] = I2 / I3

        property_dict['molar_volume'] = opt_sp_result['molar_volume']

        # build optimized xyz structure
        xyz_lines = [str(len(opt_sp_result['atomic_numbers'])), subgraph_name]
        for atomic_number, position in zip(opt_sp_result['atomic_numbers'], opt_sp_result['geometric_data']):
            
            xyz_line = []
            xyz_line.append(element_identifiers[atomic_number - 1])
            xyz_line.append('{:.4f}'.format(position[0]))
            xyz_line.append('{:.4f}'.format(position[1]))
            xyz_line.append('{:.4f}'.format(position[2]))

            xyz_lines.append(' '.join(xyz_line))
        
        xyz = '\n'.join(xyz_lines)

        # append results
        opt_sp_xyzs.append(xyz)
        opt_sp_data.append(property_dict)

# write xyzs
with open('ligands_stable_opt_xyzs.xyz', 'w') as fh:
    fh.write('\n\n'.join(opt_sp_xyzs))

# write electronic data
pd.DataFrame(opt_sp_data).to_csv('opt+sp_summary.csv', index=False, sep=';')
