import copy
from math import isclose

import numpy as np


def flatten_2d_list(list: list):

    """Flattens a 2d list to 1d.

    Arguments:
        list (list[list[]]): The 2d list.

    Returns:
        list[]: The flattend list.
    """

    return [item for sublist in list for item in sublist]

def contains_duplicates(list: list):

    """Checks for duplicates within a given list.

    Arguments:
        list (list[]): The list.

    Returns:
        bool: Flag indicating whether there duplicates or not.
    """

    if len(list) == len(set(list)):
        return False
    return True

def merge_lists_with_common_elements(lists):

    """Merges elements of a list together if they contain common elements. Removes duplicates.

    Arguments:
        lists (list[list[]]): A 2d list.

    Returns:
        list[list[]]: The 2d list with sublist containing common elements merged.
    """

    flat_list = flatten_2d_list(lists)
    while contains_duplicates(flat_list):

        for element in flat_list:
            if flat_list.count(element) > 1:
                duplicate_element = element

        duplicate_list_idx = []
        for i, li in enumerate(lists):
            if duplicate_element in li:
                duplicate_list_idx.append(i)

        merged = lists[duplicate_list_idx[0]]
        for idx in duplicate_list_idx[1:]:
            merged = list(set(merged + lists[idx]))

        new_lists = [merged]
        for i, li in enumerate(lists):
            if i not in duplicate_list_idx:
                new_lists.append(li)
        
        lists = new_lists
        flat_list = flatten_2d_list(lists)

    return lists
