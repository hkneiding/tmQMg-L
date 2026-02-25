import pandas as pd


# load ligand library
ligands_misc = pd.read_csv('../1-extraction/ligands_misc_info.csv', sep=';')[['name', 'metal_bond_node_idx_groups']]
# load singlepoint summary
single_point_summary = pd.read_csv('../2-singlepoints/sp_summary.csv', sep=';')

stable_data = []
for ligand in ligands_misc.to_dict(orient="records"):

    lowest_energy = None
    lowest_energy_subgraph_name = None
    for subgraph_name in eval(ligand['metal_bond_node_idx_groups']).keys():

        if not subgraph_name in single_point_summary['occurrence'].tolist():
            continue

        subgraph_energy = single_point_summary.loc[single_point_summary['occurrence'] == subgraph_name]['energy'].item()
        if lowest_energy is None or subgraph_energy < lowest_energy:
            lowest_energy = subgraph_energy
            lowest_energy_subgraph_name = subgraph_name

    if lowest_energy is None:
        continue

    stable_data.append({
        'name': ligand['name'],
        'stable_occurrence': lowest_energy_subgraph_name
    })

pd.DataFrame(stable_data).to_csv('stable.csv', index=False, sep=';')
