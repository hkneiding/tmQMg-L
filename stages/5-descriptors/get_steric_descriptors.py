import os
import numpy as np
import pandas as pd
from tqdm import tqdm

from morfeus import BuriedVolume, ConeAngle, SolidAngle, SASA, read_xyz


# This script calculates a selection of steric descriptors
# for each ligand using the MORFEUS package. Some descriptors
# are calculated in relation to the metal they are connected
# to (exact cone angle, buried volume, solid angle), others
# only basedon the organic part of the ligand (solvent
# accessible area and volume). For the latter the descriptors
# are calculated for both: the most stable occurrence as it
# appears in the TMC and the DFT optimised structure.

transition_metal_identifiers = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                                'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
                                'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
                                'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os',
                                'Ir', 'Pt', 'Au', 'Hg']
tmp_file_name = './temp.xyz'

# get tmc xyzs
with open('./../1-extraction/tmc_xyzs.xyz', 'r') as fh:
    tmc_xyzs = {xyz.split('\n')[1]: xyz for xyz in fh.read().split('\n\n')}

# get ligands xyzs
with open('./../1-extraction/ligands_xyzs.xyz', 'r') as fh:
    ligand_xyzs = {xyz.split('\n')[1]: xyz for xyz in fh.read().split('\n\n')}

# get optimized ligands xyzs
with open('./../4-optimizations/ligands_stable_opt_xyzs.xyz', 'r') as fh:
    ligand_xyzs_free = {xyz.split('\n')[1]: xyz for xyz in fh.read().split('\n\n')}

# get stable ligand
stable_ligands = pd.read_csv('../3-stable/stable.csv', sep=';')

# result variable
steric_descriptor_list = []

for name, stable_occurrence_name in tqdm(zip(stable_ligands['name'], stable_ligands['stable_occurrence'])):

    # retrieve ligand xyz stable
    ligand_xyz = ligand_xyzs[stable_occurrence_name].strip()

    # retrieve tmc xyz
    tmc_xyz = tmc_xyzs[stable_occurrence_name.split('-')[0]]

    # piece together ligand and tmc metal
    reduced_tmc_xyz_lines = []
    # get metal element
    for line in tmc_xyz.split('\n'):

        line_split = line.split()
        if line_split[0] in transition_metal_identifiers:
            reduced_tmc_xyz_lines.append(line)
            break
    # get ligand
    reduced_tmc_xyz_lines.extend(ligand_xyz.split('\n')[2:])

    # --- metal based descriptors --- #

    # write to temporary file
    with open(tmp_file_name, 'w') as fh:
        fh.write('\n'.join([str(len(reduced_tmc_xyz_lines)), ''] + reduced_tmc_xyz_lines))

    # read xyz into MORFEUS
    elements, coordinates = read_xyz(tmp_file_name)

    try:
        # calculate cone angle
        exact_cone_angle = ConeAngle(elements, coordinates, 1, method="internal").cone_angle

    except Exception as e:
        print('Failed exact cone angle:', e)
        print(str(len(reduced_tmc_xyz_lines)) + '\n\n' + '\n'.join(reduced_tmc_xyz_lines))
        print('\n')
        exact_cone_angle = None
    try:
        # calculate buried volume
        buried_volume = BuriedVolume(elements, coordinates, 1).fraction_buried_volume
    except Exception as e:
        print('Failed buried volume:', e)
        print(str(len(reduced_tmc_xyz_lines)) + '\n\n' + '\n'.join(reduced_tmc_xyz_lines))
        print('\n')
        buried_volume = None

    try:
        # calculate solid angle
        solid_angle_calculator = SolidAngle(elements, coordinates, 1)
        solid_angle = solid_angle_calculator.solid_angle
        solid_cone_angle = solid_angle_calculator.cone_angle
        G_parameter = solid_angle_calculator.G
    except Exception as e:
        print('Failed solid angle:', e)
        print(str(len(reduced_tmc_xyz_lines)) + '\n\n' + '\n'.join(reduced_tmc_xyz_lines))
        print('\n')
        solid_angle = None
        solid_cone_angle = None
        G_parameter = None

    # --- metal free descriptors stable --- #

    # write to temporary file
    with open(tmp_file_name, 'w') as fh:
        fh.write(ligand_xyz)

    # read xyz into MORFEUS
    elements, coordinates = read_xyz(tmp_file_name)

    # calculate solvent accessible area and volume
    sasa_calculator_stable = SASA(elements, coordinates)
    sasa_area_stable = sasa_calculator_stable.area
    sasa_volume_stable = sasa_calculator_stable.volume

    # --- metal free descriptors free --- #

    try:
        # retrieve ligand xyz free
        ligand_xyz_free = ligand_xyzs_free[stable_occurrence_name]

        # write to temporary file
        with open(tmp_file_name, 'w') as fh:
            fh.write(ligand_xyz_free)

        # read xyz into MORFEUS
        elements, coordinates = read_xyz(tmp_file_name)

        # calculate solvent accessible area and volume
        sasa_calculator_free = SASA(elements, coordinates)
        sasa_area_free = sasa_calculator_free.area
        sasa_volume_free = sasa_calculator_free.volume
    except KeyError as e:
        print('No optimized free ligand found')
        sasa_area_free = None
        sasa_volume_free = None

    steric_descriptor_list.append({
        'name': name,
        'exact_cone_angle': exact_cone_angle,
        'buried_volume': buried_volume,
        'solid_angle': solid_angle,
        'solid_cone_angle': solid_cone_angle,
        'G_parameter': G_parameter,
        'sasa_area_stable': sasa_area_stable,
        'sasa_volume_stable': sasa_volume_stable,
        'sasa_area_free': sasa_area_free,
        'sasa_volume_free': sasa_volume_free
    })

steric_descriptor_df = pd.DataFrame(steric_descriptor_list)
steric_descriptor_df.to_csv('./steric_descriptors.csv', sep=';', index=False)

# clean up tmp file
os.remove(tmp_file_name)
