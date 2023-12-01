
# tmQMg-L

This repository contains the data files of the ligands dataset tmQMg-L containing 29k ligands extracted from the Cambridge Structural Database. The ligands come with their atomic positions, metal-coordinating atom indices and corresponding formal charges. Electronic, steric and cheminformatics descriptors have been calculated for each ligand and are included as well. Details on how the data was compiled can be found in the corresponding publication [Directional Multiobjective Optimization of Metal Complexes at the Billion-Scale with the tmQMg-L Dataset and PL-MOGA Algorithm](https://chemrxiv.org/engage/chemrxiv/article-details/651051d4ed7d0eccc32252ea).

![tmQMg-L_Figure](tmQMg-L.png)

## Data

###### [ligands_misc_info.csv](ligands_misc_info.csv)
- The main ligand file containing information about IDs, SMILES, stoichiometry, occurrence and metal-coordinating atom indices.
- The column `parent_metal_occurrences` contains for all ligands a serialized Python dictionary with the different metal elements as keys and lists of occurrence names as values. For example, an entry might read ```{'Ni': ['XXYYZZ-subgraph-0', 'XXYYZZ-subgraph-2'], 'Pt': ['ZZYYXX-subgraph-0']}```, which would mean that this particular ligand has in total three occurrences where two of them come from the TMC with CSD code XXYYZZ with Nickel as the metal center and one of them comes from the TMC with CSD code ZZYYXX with Platinum as the metal center.
- The column `metal_bond_node_idx_groups` contains for all ligands a serialized Python dictionary with the different occurrences as keys and lists of the corresponding metal coordinating atom indices as values. This is necessary because the ordering of atoms in the xyz files of occurrences of the same ligand differ. For example, an entry might read ```{'XXYYZZ-subgraph-0': [[0], [2,3]], 'XXYYZZ-subgraph-2': [[1],[3,5]], 'ZZYYXX-subgraph-0': [[5],[1,4]]}``` which would mean that this particular ligand has three occurrences, each with different metal coordinating atom indices. Note that the coordinating atom indices are stored as a list of lists to denote denticity and hapticity. Metal coordinating atom indices in one sublist are haptic meaning that they are continuous in the ligand whereas indices in different sublists are dentic with respect to each other.

###### [ligands_fingerprints.csv](ligands_fingerprints.csv)
- The ligand fingerprints for all ligands containing information such as the charge, number of atoms and their type.
- The column `charge` contains the NBO derived formal charge as described in the publication.
- The columns `n_metal_bound`, `n_dentic_bound`, and `n_haptic_bound` refer to the number of atoms bound to the metal center, the number of which are dentic, and the number of which are haptic, respectively.
- The column `n_atoms` contains the total number of atoms in the ligand. For each element, the number of occurrences in the ligand is listed in column `X` where X denotes the element symbol.
- The columns `dentic_X` and `haptic_X` where X denotes the element symbol contain the number of dentic/haptic bound atoms of that element.

###### [ligands_descriptors.csv](ligands_descriptors.csv)
- The calculated RDKit, steric and electronic descriptors for all ligands as described in the publication.
- Some properties were specifically calculated based on either the geometry of the most stable occurrence or for the gas phase optimized structure. The column names reflect this with the prefix `L*` to refer to the most stable occurrence and the prefix `L_free` to refer to the relaxed structure. Properties without a prefix were simply calculated based on the ligands SMILES string.

###### [stable.csv](stable.csv)
- List of all ligands and their most stable occurrence.

###### [xyz/](xyz/)
- Directory containing the geometries of all ligands ([xyz/ligands_xyzs.xyz](xyz/ligands_xyzs.xyz)), only the stable ligands ([xyz/ligands_stable_xyzs.xyz](xyz/ligands_stable_xyzs.xyz)) and the optimized stable ligands ([xyz/ligands_stable_xyzs_opt.xyz](xyz/ligands_stable_xyzs_opt.xyz)).
- With Python the xyzs can easily be loaded as a dictionary with the occurrence names as keys and the xyzs as values using the following code snippet:
```
xyzs = {}
with open('./xyz/ligands_xyzs.xyz)', 'r') as fh:
	for xyz in fh.read().split('\n\n'):
		xyzs[xyz.split('\n')[1]] = xyz
```
###### [descriptors/](descriptors/)
- Directory containing the RDKit, steric and electronic descriptors for all ligands in separate files, the scripts to create them, and a script to merge them into one.
---

[![CC BY NC 4.0][cc-by-nc-image]][cc-by-nc]

This work is licensed under a
[Creative Commons Attribution-NonCommercial 4.0 International License][cc-by-nc].

[cc-by-nc]: http://creativecommons.org/licenses/by-nc/4.0/
[cc-by-nc-image]: https://i.creativecommons.org/l/by-nc/4.0/88x31.png
