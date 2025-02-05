import os
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm

from HyDGL import GraphGenerator, GraphGeneratorSettings, QmData, Graph, ElementLookUpTable, FileHandler
from HyDGL.enums import NodeFeature, EdgeFeature, HydrogenMode, BondOrderType, EdgeType

from tools import *
from chem import *
from smiles import *
from neighbour import *

# noble gas atomic numbers
NOBLE_GAS_ATOMIC_NUMBERS = [0, 2, 10, 18, 36, 54, 86]


def remove_metal_center(graph: Graph):

    """Returns a graph with the metal centre removed.

    Arguments:
        graph (Graph): The molecular Graph of which the metal centre should be removed.

    Returns:
        Graph: The graph without metal centre.
    """

    nodes = []
    # build new node list without metal
    for i, node in enumerate(graph.nodes):
        # determine metal node index
        if node.features['atomic_number'] in ElementLookUpTable.transition_metal_atomic_numbers:
            metal_node_idx = i
        else:
            nodes.append(node)

    # check the number of metal atoms in molecular graph
    if len(nodes) == len(graph.nodes):
        raise RuntimeError('No metal atom in molecular graph.')
    elif len(graph.nodes) - len(nodes) > 1:
        raise RuntimeError('There are more than one (1) metal centres in molecular graph.')

    # variables to store node ids of atoms that bind to metal
    metal_bond_node_ids = []
    metal_bond_node_distances = []
    metal_bond_nbo_bond_node_ids = {}

    edges=[]
    # build new edge list without edge to metal
    for edge in graph.edges:
        # determine nodes with edges to metal
        if metal_node_idx in edge.node_indices:
            metal_bond_node_ids.append(graph.nodes[edge.node_indices[(edge.node_indices.index(metal_node_idx) + 1) % 2]].id)
            metal_bond_node_distances.append(edge.features['bond_distance'])

            if edge.label == 'NBO':
                metal_bond_nbo_bond_node_ids[metal_bond_node_ids[-1]] = edge.features['n_bn']

        else:
            for i in range(len(edge._node_indices)):
                if edge._node_indices[i] > metal_node_idx:
                    edge._node_indices[i] -= 1
            edges.append(edge)

    # build new graph
    return Graph(nodes, edges, {}, {}, meta_data={
        'id': graph.id,
        'metal_node_idx': metal_node_idx,
        'metal_bond_node_ids': metal_bond_node_ids,
        'metal_bond_node_distances': metal_bond_node_distances,
        'metal_bond_nbo_bond_node_ids': metal_bond_nbo_bond_node_ids
    })

def calculate_total_valency(graph: Graph):

    """Calculates the total classical valency based on atomic numbers.

    Arguments:
        graph (Graph): The molecular graph.

    Returns:
        int: The total valency.
    """

    total_valency = 0
    for node in graph.nodes:
        
        # accumulate total valence total
        total_valency +=  min([node.features['atomic_number'] - _ if node.features['atomic_number'] - _ > 0 else np.inf for _ in NOBLE_GAS_ATOMIC_NUMBERS])

        # d orbital correction for elements after d block
        if ElementLookUpTable.get_atomic_number(node.label) >= 31 and ElementLookUpTable.get_atomic_number(node.label) <= 36 or \
         ElementLookUpTable.get_atomic_number(node.label) >= 49 and ElementLookUpTable.get_atomic_number(node.label) <= 54 or \
         ElementLookUpTable.get_atomic_number(node.label) >= 81 and ElementLookUpTable.get_atomic_number(node.label) <= 86:
            total_valency -= 10

    return total_valency

def calculate_nbo_electron_occupation(graph: Graph, metal_bond_nbo_bond_node_ids: list):
    
    """Calculates the electron occupation from NBO LPs and BDs.

    Arguments:
        graph (Graph): The molecular graph.
        metal_bond_nbo_bond_node_ids (list[int]): The list of node ids of nodes involved in metal NBO bonds.

    Returns:
        int: The NBO electron occupation.
    """

    # calculate nbo electron occupation
    nbo_electron_occupation = 0
    
    # iterate through nodes
    for node in graph.nodes:

        # add 2 electrons per lone pair
        nbo_electron_occupation += 2 * node.features['n_lone_pairs']
        
        # add 2 electrons per NBO bond to metal center
        if node.id in metal_bond_nbo_bond_node_ids.keys():
            nbo_electron_occupation += 2 * metal_bond_nbo_bond_node_ids[node.id]

    # iterate through edges
    for edge in graph.edges:
        # add 2 electrons per NBO bond
        if edge.features['nbo_type'] == 'BD':
            nbo_electron_occupation += 2 * edge.features['n_bn']
        # return None if 3C interaction is detected
        elif edge.features['nbo_type'] == '3C':
            return None     

    return nbo_electron_occupation       

def calculate_predicted_charge(graph: Graph, metal_bond_nbo_bond_node_ids: list):

    """Calculates the predicted charge for a molecular graph.

    Arguments:
        graph (Graph): The molecular graph.
        metal_bond_nbo_bond_node_ids (list[int]): The list of node ids of nodes involved in metal NBO bonds.
    
    Returns:
        int: The predicted charge
    """

    # calculate predicted charge
    subgraph_total_valency = calculate_total_valency(graph)
    nbo_electron_occupation = calculate_nbo_electron_occupation(graph, metal_bond_nbo_bond_node_ids)

    if nbo_electron_occupation is None:
        return None

    return subgraph_total_valency - nbo_electron_occupation

def get_modified_charge(predicted_charge: int, xyz: str, metal_bond_node_idx_groups: list[list[int]]):

    """Determines the modified charge of a given ligand based on a set of heuristic rules.

    Arguments:
        predicted_charge (int): The predicted charge.
        fingerprint (list[int]): The fingerprint of the ligand.
    
    Returns:
        int: The predicted charge
    """

    # set C2-eta and C4-eta ligands with charge -2 to 0
    # caveat: only ligands where all metal connecting carbons have 3 substituents

    # check charge
    if predicted_charge == -2:
        
        neighbour_list = get_neighbour_list(xyz)
        
        for metal_bond_node_idx_group in metal_bond_node_idx_groups:

            # look for haptic carbon coordination modes
            has_only_carbon = True
            xyz_split = xyz.split('\n')
            for i, idx in enumerate(metal_bond_node_idx_group):
                if xyz_split[idx + 2].split(' ')[0] != 'C':
                    has_only_carbon = False
                    break

            # skip for non-carbon
            if not has_only_carbon:
                continue
            
            # check for eta2
            if len(metal_bond_node_idx_group) == 2:
                # check for substitution pattern
                if len(get_neighbours(metal_bond_node_idx_group[0], neighbour_list)) == 3 and \
                   len(get_neighbours(metal_bond_node_idx_group[1], neighbour_list)) == 3:
                    return 0

            # check for eta4
            if len(metal_bond_node_idx_group) == 4:
                # check for substitution pattern
                if len(get_neighbours(metal_bond_node_idx_group[0], neighbour_list)) == 3 and \
                   len(get_neighbours(metal_bond_node_idx_group[1], neighbour_list)) == 3 and \
                   len(get_neighbours(metal_bond_node_idx_group[2], neighbour_list)) == 3 and \
                   len(get_neighbours(metal_bond_node_idx_group[3], neighbour_list)) == 3:
                    return 0

    return predicted_charge

def get_metal_bond_idx(nodes, metal_bond_node_ids, metal_bond_node_distances):

    """Gets the indices of atoms bonded to the metal centre.

    Arguments:
        nodes (list[Node]): The list of nodes.
        metal_bond_node_ids (list[int]): The ids of atoms bonded to the metal centre.
        metal_bond_node_distances (list[list]): The distances of atoms bonded to the metal centre.

    Returns:
        list[int]: The indices of atoms bonded to the metal centre.
        list[int]: Their corresponding distances to the metal centre.
        list[int]: Their corresponding ids.
    """

    # get subgraph meta data
    subgraph_metal_bond_node_idx = []
    subgraph_metal_bond_node_ids = []
    subgraph_metal_bond_node_distances = []

    for i, node in enumerate(nodes):

        # get indices of nodes that bind to metal
        if node.id in metal_bond_node_ids:
            subgraph_metal_bond_node_idx.append(i)
            subgraph_metal_bond_node_ids.append(node.id)
            subgraph_metal_bond_node_distances.append(metal_bond_node_distances[metal_bond_node_ids.index(node.id)])

    return subgraph_metal_bond_node_idx, subgraph_metal_bond_node_distances, subgraph_metal_bond_node_ids

def get_metal_bond_idx_groups(graph, metal_bond_node_idx):

    """Gets the indices of atoms bound to the metal grouped in terms of contiguous nodes.

    Arguments:
        graph (Graph): The molecular graph.
        metal_bond_node_idx (list[int]): The list of node indices connecting to the metal.

    Returns:
        list[list[int]]: The list of contiguous indices sublists.
    """

    # get the adjacency matrix
    adj = graph.get_adjacency_matrix()

    node_id_groups = []
    for i in range(len(metal_bond_node_idx)):
        node_id_group = [metal_bond_node_idx[i]]
        for j in range(len(metal_bond_node_idx)):
            if adj[metal_bond_node_idx[i]][metal_bond_node_idx[j]] == 1:
                node_id_group.append(metal_bond_node_idx[j])
        if sorted(node_id_group) not in node_id_groups:
            node_id_groups.append(sorted(node_id_group))

    # merge sublists with common elements.
    merged_node_id_groups = merge_lists_with_common_elements(node_id_groups)

    return merged_node_id_groups

def node_match(n1_features, n2_features):
    
    """Compares two nodes based on their atomic number.

    Arguments:
        n1_features (dict): The feature dict of node 1.
        n2_features (dict): The feature dict of node 2.

    Returns:
        bool: Flag indicating whether nodes match or not.
    """

    if n1_features['feature_atomic_number'] == n2_features['feature_atomic_number']:
        return True

    return False

def get_fingerprint_dict(xyz, charge, metal_bond_node_idx_groups, is_alternative_charge=0):

    """Gets the fingerprint of a ligand as a dictionary.

    Arguments:
        xyz (str): The xyz data.
        charge (int): The assigned charge.
        metal_bond_node_idx_groups (list[int]): The list of bond node index groups.
        is_alternative_charge (int): The flag indicating whether or not the assigned charge is canonical or not.

    Returns:
        dict: The fingerprint as dictionary.
    """

    # elements to consider
    elements = ['H', 'B', 'C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'As', 'Se', 'Br', 'I']
    # atoms in molecule
    atoms = [xyz.strip().split('\n')[2:][i].split(' ')[0] for i in range(len(xyz.strip().split('\n')[2:]))]

    chem_fingerprint = {key: atoms.count(key) for key in elements}

    # get counts of elements that are metal-bound individually
    dentic_fingerprint = {'dentic_'+ key: 0 for key in elements}

    # counts of elements that are metal-bound contiguously
    haptic_fingerprint = {'haptic_'+ key: 0 for key in elements}
    
    # fill dentic and haptic fingerprint
    for metal_bond_node_idx_group in metal_bond_node_idx_groups:

        # dentic atoms
        if len(metal_bond_node_idx_group) == 1:
            dentic_fingerprint['dentic_' + atoms[metal_bond_node_idx_group[0]]] += 1
        # haptic atoms
        else:
            for haptic_atom_idx in metal_bond_node_idx_group:
                haptic_fingerprint['haptic_' + atoms[haptic_atom_idx]] += 1

    # build fingerprint    
    fingerprint = {}
    fingerprint['charge'] = charge
    fingerprint['n_atoms'] = len(atoms)
    fingerprint['n_metal_bound'] = sum(dentic_fingerprint.values()) + sum(haptic_fingerprint.values())
    fingerprint['n_dentic_bound'] = sum(dentic_fingerprint.values())
    fingerprint['n_haptic_bound'] = sum(haptic_fingerprint.values())
    fingerprint['is_alternative_charge'] = is_alternative_charge
    
    fingerprint = fingerprint | chem_fingerprint | dentic_fingerprint | haptic_fingerprint
    
    return fingerprint

def get_fingerprint(xyz, charge, metal_bond_node_idx_groups, is_alternative_charge=0):

    """Gets the fingerprint of a ligand as a list.

    Arguments:
        xyz (str): The xyz data.
        charge (int): The assigned charge.
        metal_bond_node_idx_groups (list[int]): The list of bond node index groups.
        is_alternative_charge (int): The flag indicating whether or not the assigned charge is canonical or not.

    Returns:
        list: The fingerprint as list.
    """

    return list(get_fingerprint_dict(xyz, charge, metal_bond_node_idx_groups, is_alternative_charge=is_alternative_charge).values())

def get_nbo_extracted_ligands(data_dir: str):

    """Extracts ligands based on their u-NatQG and NBO data for determining their charge.

    Arguments:
        data_dir (str): The path to the directory containing the json formatted raw data.

    Returns:
        pd.DataFrame: The extracted ligands and their (meta) information.
    """

    files = [x for x in os.listdir(data_dir)]

    # set up graph generator
    ggs = GraphGeneratorSettings(
        node_features=[NodeFeature.ATOMIC_NUMBER, NodeFeature.LONE_PAIR_AVERAGE],
        edge_types=[EdgeType.BOND_ORDER_METAL, EdgeType.NBO_BONDING_ORBITALS],
        bond_order_mode=BondOrderType.WIBERG,
        bond_threshold=0.5,
        bond_threshold_metal=0.05,
        hydrogen_mode=HydrogenMode.EXPLICIT,
        hydrogen_count_threshold=0.5,
        edge_features=[EdgeFeature.BOND_DISTANCE, EdgeFeature.NBO_TYPE, EdgeFeature.BOND_ORBITAL_AVERAGE],
        sopa_edge_features=[],
        graph_features=[],
        targets=[],
        sopa_contribution_threshold=None,
        sopa_interaction_threshold=None,
        sopa_resolution_mode=None,
        max_bond_distance=2.75
    )
    gg = GraphGenerator(ggs)

    ligands = {}

    for file in tqdm(files):

        # skip files that are not in JSON format
        if file.split('.')[-1] != 'json':
            continue   

        # TMCs for which to skip ligand extraction
        skip = [
            'SURLEX', 'DAKZEU', 'MCTPCR', 'AHIDEB', 'NANPYZ', 'PDTCNI', 'MQUPTI', 'PHGFAN',
            'SOXNOI', 'ACPXNI', 'HGTHUR', 'PHGBAN', 'COTNUU', 'CPMONC', 'GECHAX', 'NABGON',
            'JAPROH', 'TEFZOT', 'IPTCNI', 'DAMTIU', 'VIYKET', 'IPTCHG', 'DEFJUT', 'CMEEAM',
            'JEBPUB', 'SANLEY', 'CFMNTV', 'KORGED', 'MEBCUU', 'MECBHG', 'XEDREE', 'NEBQIU'
        ]
        if file.split('.')[0] in skip:
            continue

        qm_data = QmData.from_dict(FileHandler.read_dict_from_json_file(data_dir + file))
        # generate graph
        graph = gg.generate_graph(qm_data)

        # skip if disconnected
        if not graph.is_connected():
            print(graph.id)
            continue

        # get graph without metal center
        new_graph = remove_metal_center(graph)
        metal_node_idx = new_graph.meta_data['metal_node_idx']

        # get ligand subgraphs
        subgraphs = new_graph.get_disjoint_sub_graphs()

        # iterate through subgraphs and determine meta data
        for subgraph in subgraphs:

            # subgraph.visualise()
            
            # get metal bond ids and distances
            subgraph_metal_bond_node_idx, subgraph_metal_bond_node_distances, subgraph_metal_bond_node_ids = get_metal_bond_idx(subgraph.nodes, 
                                                                                    new_graph.meta_data['metal_bond_node_ids'], 
                                                                                    new_graph.meta_data['metal_bond_node_distances'])

            metal_bond_node_idx_groups = get_metal_bond_idx_groups(subgraph, subgraph_metal_bond_node_idx)

            # flag to determine whether to make a new entry for this ligand
            make_new_entry = True

            # calculate spectrum and hash it as a heuristic to rule out all ligands that are not cospectral
            spectrum_hash = hash(tuple(np.around(subgraph.get_spectrum(), decimals=5)))

            # get ids of cospectral graphs
            cospectral_graphs = [k for k,v, in ligands.items() if v['spectrum_hash'] == spectrum_hash]

            # obtain ligand xyz
            xyz = subgraph.get_xyz_data().strip()

            # obtain fingerprint for ligand
            fingerprint_sub = get_fingerprint_dict(
                xyz,
                charge=None,
                metal_bond_node_idx_groups=metal_bond_node_idx_groups
            )

            # get charge
            predicted_charge = calculate_predicted_charge(subgraph, new_graph.meta_data['metal_bond_nbo_bond_node_ids'])
            if predicted_charge is None:
                continue

            # modify charge based on heuristic rules
            predicted_charge = get_modified_charge(predicted_charge, xyz, metal_bond_node_idx_groups)

            # if spectrum_hash in [ligands[subgraph_id]['spectrum_hash'] for subgraph_id in ligands.keys()]:
            if len(cospectral_graphs) > 0:

                # loop through cospectral graphs to find possible identities 
                for key in cospectral_graphs:

                    # obtain fingerprint for cospectral ligand
                    fingerprint_lig = get_fingerprint_dict(
                        ligands[key]['xyzs'][key],
                        charge=None,
                        metal_bond_node_idx_groups=ligands[key]['metal_bond_node_idx_groups'][key]
                    )

                    # check if fingerprints match
                    if fingerprint_sub == fingerprint_lig:
                        # check for isomorphy via bruteforce VF2
                        if nx.is_isomorphic(subgraph.get_networkx_graph_object(), ligands[key]['networkx'], node_match=node_match):

                            if graph.nodes[metal_node_idx].label in ligands[key]['parent_metal_occurrences'].keys():
                                ligands[key]['parent_metal_occurrences'][graph.nodes[metal_node_idx].label].append(subgraph.id)
                            else:
                                ligands[key]['parent_metal_occurrences'][graph.nodes[metal_node_idx].label] = [subgraph.id]
                            
                            if predicted_charge in ligands[key]['charge_occurrences'].keys():
                                ligands[key]['charge_occurrences'][predicted_charge].append(subgraph.id)
                            else:
                                ligands[key]['charge_occurrences'][predicted_charge] = [subgraph.id]

                            ligands[key]['xyzs'][subgraph.id] = xyz
                            ligands[key]['metal_bond_node_idx_groups'][subgraph.id] = metal_bond_node_idx_groups

                            # set flag to false
                            make_new_entry = False

            # if it does not exist make new entry
            if make_new_entry:

                smiles, smiles_metal_bond_node_idx_groups = get_smiles_and_smiles_metal_connecting_index_group(xyz, metal_bond_node_idx_groups)

                # add ligand to dictionary
                ligands[subgraph.id] = {
                    
                    # topological info
                    'xyzs': {subgraph.id: xyz},
                    'smiles': smiles,
                    'stoichiometry': get_stoichiometry(xyz),
                    'smiles_metal_bond_node_idx_groups': smiles_metal_bond_node_idx_groups,
                    'metal_bond_node_idx_groups': {subgraph.id: metal_bond_node_idx_groups},
                    'parent_metal_occurrences': {graph.nodes[metal_node_idx].label: [subgraph.id]},
                    'spectrum_hash': spectrum_hash,

                    # electronic info
                    'charge_occurrences': {predicted_charge: [subgraph.id]},

                    # graph data
                    'networkx': subgraph.get_networkx_graph_object()
                }

    # cast to pandas DF
    ligands = pd.DataFrame(ligands).T
    # second pass over dataset to split based on charge
    # also gets rid of some unnecessary data
    new_ligands = {}
    for i, ligand in enumerate(ligands.to_dict(orient='records')):

        # determine the main charge for this ligand
        for j, item in enumerate(sorted(ligand['charge_occurrences'].items(), key=lambda x: len(x[1]), reverse=True)):

            # initialise main charge
            if j == 0:
                main_charge = item
                continue
            
            # check if another charge with same occurrence exists
            if item[1] == main_charge[1]:
                # replace if a lower charge for that occurrence is found
                if item[0] < main_charge[0]:
                    main_charge = item
                continue
            else:
                break

        # list to hold all charges to consider for this ligand
        charges_to_consider = []
        # loop through ligands to find alternative charge assignments
        last = None
        for j, item in enumerate(sorted(ligand['charge_occurrences'].items(), key=lambda x: x[1], reverse=True)):

            # skip if main charge is encountered
            if item[0] != main_charge[0]:

                # ckeck if alternative charge assignment occupies at least 25% of the occurrences
                if len(item[1]) / sum([len(ligand['charge_occurrences'][key]) for key in ligand['charge_occurrences'].keys()]) > 0.25:
                    
                    # check for same occupations
                    if last is not None and item[1] == last[1]:
                        # if a lower charge for that occupation is found, assign it
                        if item[0] < last[0]:
                            charges_to_consider[-1] = item[0]

                    # check that charge is 0 or negative
                    # check that difference to main charge is 2
                    if item[0] <= 0 and np.abs(main_charge[0] - item[0]):
                        charges_to_consider.append(item[0])
            
            # store last item
            last = item

        charges_to_consider.insert(0, main_charge[0])
        # iterate through all charges to consider and add individual ligands
        for j, charge_to_consider in enumerate(charges_to_consider):

            # add ligand to dictionary
            new_ligands['ligand' + str(i) + '-' + str(j)] = {
                
                # topological info
                'xyzs': {ligand_name: ligand['xyzs'][ligand_name] for ligand_name in ligand['charge_occurrences'][charge_to_consider]},
                'smiles': ligand['smiles'],
                'stoichiometry': ligand['stoichiometry'],
                'smiles_metal_bond_node_idx_groups': ligand['smiles_metal_bond_node_idx_groups'],
                'metal_bond_node_idx_groups': {ligand_name: ligand['metal_bond_node_idx_groups'][ligand_name] for ligand_name in ligand['charge_occurrences'][charge_to_consider]},
                # stores only the parent metal occurrences that occur with that particular charge
                'parent_metal_occurrences': {parent_metal: list(set(ligand['parent_metal_occurrences'][parent_metal]) & set(ligand['charge_occurrences'][charge_to_consider])) for parent_metal in ligand['parent_metal_occurrences'].keys() if len(list(set(ligand['parent_metal_occurrences'][parent_metal]) & set(ligand['charge_occurrences'][charge_to_consider]))) != 0},
                # electronic info
                'charge': charge_to_consider,
                'is_alternative_charge': 0 if not j else 1,
                'occurrence': len(set([_.split('-')[0] for _ in ligand['charge_occurrences'][charge_to_consider]]))
            }

    new_ligands = pd.DataFrame(new_ligands).T
    
    return new_ligands


# - - - entry point - - - #
if __name__ == "__main__":

    nbo_ligands = get_nbo_extracted_ligands('/home/hkneiding/Documents/UiO/Data/tmQMg/json/')

    dataset_name = 'ligands'
    dataset_directory = './'

    xyzs = []
    fingerprints = []
    misc = []

    # extract ligands
    for i, ligand in enumerate(nbo_ligands.to_dict(orient='records')):

        # make xyz
        for key in ligand['xyzs'].keys():
            xyz_lines = ligand['xyzs'][key].split('\n')
            xyz_lines[1] = xyz_lines[1].replace('id: ', '')
            xyzs.append('\n'.join(xyz_lines))

        # make fingerprint csv
        pivot_ligand_key = list(ligand['xyzs'].keys())[0]
        fingerprint_dict = get_fingerprint_dict(ligand['xyzs'][pivot_ligand_key], ligand['charge'], ligand['metal_bond_node_idx_groups'][pivot_ligand_key], is_alternative_charge=ligand['is_alternative_charge'])
        fingerprints.append({'name': nbo_ligands.index[i]} | fingerprint_dict)
        
        # make misc csv
        misc_keys = ['smiles', 'stoichiometry', 'occurrence', 'parent_metal_occurrences', 'smiles_metal_bond_node_idx_groups', 'metal_bond_node_idx_groups']
        misc_dict = {key: ligand[key] for key in misc_keys}
        misc.append({'name': nbo_ligands.index[i]} | misc_dict)

    # save raw dataframe
    nbo_ligands.to_csv(dataset_directory + '.' +  dataset_name + '_raw.csv', sep=';')

    # write xyz to file
    with open(dataset_directory + dataset_name + '_xyzs.xyz', 'w') as fh:
        fh.write('\n\n'.join(xyzs))

    # save fingerprints df to csv
    fingerprints_df = pd.DataFrame(fingerprints)
    fingerprints_df.set_index('name')
    fingerprints_df.to_csv(dataset_directory + dataset_name + '_fingerprints.csv', sep=';')

    # save misc info df to csv
    misc_df = pd.DataFrame(misc)
    misc_df.set_index('name')
    misc_df.to_csv(dataset_directory + dataset_name + '_misc_info.csv', sep=';')
