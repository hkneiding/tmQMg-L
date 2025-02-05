def get_stoichiometry(xyz: str):

    """Gets the stoichiometry string for a given xyz structure.

    Arguments:
        xyz (str): The xyz structure.

    Returns:
        str: The stoichiometry string.
    """

    lines = xyz.strip().split('\n')

    # count elements in xyz
    element_counts = {}
    for i in range(2, len(lines)):
        element = lines[i].split(' ')[0]
        if element in element_counts.keys():
            element_counts[element] += 1
        else:
            element_counts[element] = 1

    # setup chemical formula string
    stoichiometry = ''
    for key in sorted(element_counts.keys()):

        if key == 'H' or key == 'C':
            continue

        stoichiometry += key 
        if element_counts[key] > 1:
            stoichiometry += str(element_counts[key])

    for key in ['H', 'C']:
        if key in element_counts.keys():
            if element_counts[key] > 1:
                stoichiometry = str(element_counts[key]) + stoichiometry
            stoichiometry = key + stoichiometry

    return stoichiometry
