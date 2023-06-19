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

# get xyzs of tmQMg TMCs
tmc_dir = '/home/hkneiding/Documents/UiO/Data/tmQMg/xyz/'
tmc_names = sorted([ f.name for f in os.scandir(tmc_dir) if not f.is_dir() ])
tmc_xyzs = {}
for tmc_name in tmc_names:

    with open(tmc_dir + tmc_name, 'r') as fh:
        tmc_xyzs[tmc_name.split('.')[0]] = fh.read()

# get xyzs of tmQMg ligands stable
with open('./../xyz/ligands_stable_xyzs.xyz', 'r') as fh:
    ligand_xyzs = {xyz.split('\n')[1]: xyz for xyz in fh.read().split('\n\n')}

# get xyzs of tmQMg ligands free
with open('./../xyz/ligands_stable_xyzs_opt.xyz', 'r') as fh:
    ligand_xyzs_free = {xyz.split('\n')[1]: xyz for xyz in fh.read().split('\n\n')}

# get most stable ligand
with open('./../stable.txt', 'r') as fh:
    ligand_stable_names = [[_.split(' | ')[0], _.split(' | ')[1]] for _ in fh.read().strip().split('\n')]

#print(ligand_stable_names)
#print(tmc_xyzs)
#print(ligand_xyzs)

# result variable
steric_descriptor_list = []

for ligand_stable_name in tqdm(ligand_stable_names):

    # retrieve ligand xyz stable
    ligand_xyz = ligand_xyzs[ligand_stable_name[1]].strip()
    # retrieve ligand xyz free
    ligand_xyz_free = ligand_xyzs_free[ligand_stable_name[1]]

    # retrieve tmc xyz
    tmc_xyz = tmc_xyzs[ligand_stable_name[1].split('-')[0]]

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

    # write to temporary file
    with open(tmp_file_name, 'w') as fh:
        fh.write(ligand_xyz_free)

    # read xyz into MORFEUS
    elements, coordinates = read_xyz(tmp_file_name)

    # calculate solvent accessible area and volume
    sasa_calculator_free = SASA(elements, coordinates)
    sasa_area_free = sasa_calculator_stable.area
    sasa_volume_free = sasa_calculator_stable.volume


    steric_descriptor_list.append({
        'name': ligand_stable_name[0],
        'stable_occurrence_name': ligand_stable_name[1],
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

    #print(steric_descriptor_list)
    #exit()

steric_descriptor_df = pd.DataFrame(steric_descriptor_list)
steric_descriptor_df.to_csv('./steric_descriptors.csv', sep=';', index=False)

# clean up tmp file
os.remove(tmp_file_name)
