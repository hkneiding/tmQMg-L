import numpy as np


def euclidean_distance(v1: list[float], v2: list[float]):

    """Calculates the euclidean distance between two points.

    Arguments:
        v1 (list[float]): The first point.
        v2 (list[float]): The second point.

    Returns:
        float: The euclidean distance.
    """

    return np.linalg.norm(np.array(v1) - np.array(v2))


def get_neighbour_list(xyz: str, threshold=1.6):

    """Gets the neighbour list of a given xyz structure.

    Arguments:
        xyz (str): The xyz structure.

    Returns:
        list[list[int]]: The neighbour list.
    """

    xyz_split = xyz.split('\n')
    xyz_split = xyz_split[2:]
    xyz_coord = []
    for i in range(len(xyz_split)):

        line_split = xyz_split[i].split(' ')
        line_coord = [float(line_split[1]), float(line_split[2]), float(line_split[3])]
        xyz_coord.append(line_coord)

    assert len(xyz_coord) == len(xyz.split('\n')) - 2

    adj_list = []
    for i in range(0, len(xyz_split) - 1, 1):
        for j in range(i + 1, len(xyz_split), 1):
            if euclidean_distance(xyz_coord[i], xyz_coord[j]) < threshold:
                adj_list.append([i,j])

    return adj_list

def get_neighbours(idx, adj_list):

    """Gets the neighbours of a given point based on a neighbourlist.

    Arguments:
        idx (int): The index of the query point.
        adj_list (list[list[int]]): The neighbourlist.
    Returns:
        list[int]: The list of neighbours.
    """

    neighbours = []
    for bond in adj_list:
        if idx in bond:
           neighbours.append(bond[(bond.index(idx) + 1) % 2])

    return neighbours
