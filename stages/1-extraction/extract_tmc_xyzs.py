import os

from tqdm import tqdm

from HyDGL import GraphGenerator, GraphGeneratorSettings, QmData, FileHandler
from HyDGL.enums import NodeFeature, EdgeFeature, HydrogenMode, BondOrderType, EdgeType


def get_tmc_xyzs(data_dir: str):

    """Extracts the xyz structures of TMCs.

    Arguments:
        data_dir (str): The path to the directory containing the json formatted raw data of TMCs.

    Returns:
        list[str]: A list of TMC xyz strings.
    """

    files = [x for x in os.listdir(data_dir)]

    # set up graph generator
    ggs = GraphGeneratorSettings.baseline([])
    gg = GraphGenerator(ggs)

    xyzs = []
    for file in tqdm(files):

        # skip files that are not in JSON format
        if file.split('.')[-1] != 'json':
            continue

        # read qm data
        qm_data = QmData.from_dict(FileHandler.read_dict_from_json_file(data_dir + file))
        # generate graph
        graph = gg.generate_graph(qm_data)

        xyz = graph.get_xyz_data().strip()
        xyz_lines = xyz.split('\n')
        xyz_lines[1] = xyz_lines[1].split()[1]
        xyz = '\n'.join(xyz_lines)

        xyzs.append(xyz)

    return xyzs

# specify directory containing the extracted QM data in json format
tmc_xyzs = get_tmc_xyzs('./json/')
with open('tmc_xyzs.xyz', 'w') as fh:
    fh.write('\n\n'.join(tmc_xyzs))
