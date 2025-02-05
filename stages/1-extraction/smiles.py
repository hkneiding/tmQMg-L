import copy
from math import isclose

import numpy as np


def get_smiles(xyz: str):

    """Gets the SMILES string of a given xyz structure.

    Arguments:
        xyz (str): The xyz structure.

    Returns:
        str: The SMILES string.
    """

    from openbabel import openbabel

    # setup converter for reading xyz
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("xyz", "mdl")

    # setup molecule
    mol = openbabel.OBMol()
    obConversion.ReadString(mol, xyz.split('\n')[0] + '\n\n' + '\n'.join(xyz.split('\n')[2:-1]))

    # set converter for getting smiles
    obConversion.SetInAndOutFormats("mdl", "smi")

    return obConversion.WriteString(mol).strip()

def get_smiles_and_xy_positions(xyz: str):

    """Gets the SMILES string of a given xyz structure and the xy atom positions in the order they appear in the SMILES string.

    Arguments:
        xyz (str): The xyz structure.

    Returns:
        str: The SMILES string.
        list[list[float]]: The list of xy atom positions.
    """

    from openbabel import openbabel

    # setup converter for reading xyz
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("xyz", "mdl")

    # setup molecule
    mol = openbabel.OBMol()
    obConversion.ReadString(mol, xyz.split('\n')[0] + '\n\n' + '\n'.join(xyz.split('\n')[2:]))

    # set converter for getting smiles
    obConversion.SetInAndOutFormats("mdl", "smi")
    obConversion.AddOption("x", obConversion.OUTOPTIONS)
    obConversion.AddOption("h", obConversion.OUTOPTIONS)

    smiles, xy_string = obConversion.WriteString(mol).strip().split()

    xy_string_split = xy_string.split(',')
    xy_positions = [
        [
            float(xy_string_split[i]),
            float(xy_string_split[i+1])
        ]
        for i in range(0, len(xy_string_split), 2)
    ]

    return smiles, xy_positions

def get_smiles_and_smiles_metal_connecting_index_group(xyz: str, xyz_metal_connecting_index_group: list[list[int]]):

    """Gets the SMILES string of a given xyz structure and the corresponding SMILES metal connecting indices.

    Arguments:
        xyz (str): The xyz structure.
        xyz_metal_connecting_index_group (list[list[int]]): The xyz metal connecting index group.

    Returns:
        str: The SMILES string.
        list[list[int]]: The SMILES metal connecting index group.
    """

    smiles, xy_positions = get_smiles_and_xy_positions(xyz)
    xyz_positions = [[float(_) for _ in line.split()[1:]] for line in xyz.split('\n')[2:]]

    smiles_metal_connecting_index_group = copy.deepcopy(xyz_metal_connecting_index_group)

    for i in range(len(xyz_metal_connecting_index_group)):
        for j in range(len(xyz_metal_connecting_index_group[i])):
            is_found = False
            for k in range(len(xy_positions)):

                if isclose(
                    xy_positions[k][0],
                    np.round(xyz_positions[xyz_metal_connecting_index_group[i][j]][0], decimals=4),
                    abs_tol=0.0002) \
                and isclose(
                    xy_positions[k][1],
                    np.round(xyz_positions[xyz_metal_connecting_index_group[i][j]][1],decimals=4),
                    abs_tol=0.0002):

                    smiles_metal_connecting_index_group[i][j] = k
                    is_found = True
                    break

            if not is_found:
                print('ERROR MATCHING XY AND XYZ POSITIONS')
                exit()

    return smiles, smiles_metal_connecting_index_group
