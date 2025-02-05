
# tmQMg-L

This repository contains the data files of the ligands dataset tmQMg-L containing 35k ligands extracted from the tmQMg dataset containing transition metal complexes reported in the Cambridge Structural Database. The ligands come with their atomic positions, metal-coordinating atom indices and corresponding formal charges. Electronic, steric and cheminformatics descriptors have been calculated for each ligand and are included as well. Details on how the data was compiled can be found in the corresponding publication [Directional multiobjective optimization of metal complexes at the billion-system scale](https://www.nature.com/articles/s43588-024-00616-5).

![tmQMg-L_Figure](tmQMg-L.png)

## Data

###### [ligands_misc_info.csv](ligands_misc_info.csv)
- The main ligand file containing information about IDs, SMILES, stoichiometry, occurrence and metal-coordinating atom indices.
- The column `parent_metal_occurrences` contains for all ligands a serialized Python dictionary with the different metal elements as keys and lists of occurrence names as values. For example, an entry might read ```{'Ni': ['XXYYZZ-subgraph-0', 'XXYYZZ-subgraph-2'], 'Pt': ['ZZYYXX-subgraph-0']}```, which would mean that this particular ligand has in total three occurrences where two of them come from the TMC with CSD code XXYYZZ with Nickel as the metal center and one of them comes from the TMC with CSD code ZZYYXX with Platinum as the metal center.
- The column `smiles_metal_bond_node_idx_groups` contains for all ligands the indices of metal coordinating atoms in the corresponding SMILES string. Note that they are stored as a list of lists to denote denticity and hapticity. Metal coordinating atom indices in one sublist are haptic meaning that they are continuous in the ligand whereas indices in different sublists are dentic with respect to each other.
- The column `metal_bond_node_idx_groups` contains for all ligands a serialized Python dictionary with the different occurrence names as keys and the corresponding metal coordinating atom indices as values. This distinction is necessary because the ordering of atoms in the xyz files of different occurrences of the same ligand differ. For example, an entry might read ```{'XXYYZZ-subgraph-0': [[0], [2,3]], 'XXYYZZ-subgraph-2': [[1],[3,5]], 'ZZYYXX-subgraph-0': [[5],[1,4]]}``` which would mean that this particular ligand has three occurrences, each with different metal coordinating atom indices due to different ordering of atoms in the xyzs. Note again that the coordinating atom indices are stored as a list of lists to denote denticity and hapticity. Metal coordinating atom indices in one sublist are haptic meaning that they are continuous in the ligand whereas indices in different sublists are dentic with respect to each other.

###### [ligands_fingerprints.csv](ligands_fingerprints.csv)
- The ligand fingerprints for all ligands containing information such as the charge, number of atoms and their type.
- The column `charge` contains the NBO derived formal charge as described in the publication.
- The columns `n_metal_bound`, `n_dentic_bound`, and `n_haptic_bound` refer to the number of atoms bound to the metal center, the number of which are dentic, and the number of which are haptic, respectively.
- The column `n_atoms` contains the total number of atoms in the ligand. For each element, the number of occurrences in the ligand is listed in column `X` where X denotes the element symbol.
- The columns `dentic_X` and `haptic_X` where X denotes the element symbol contain the number of dentic/haptic bound atoms of that element.
- The column `is_alternative_charge` is a boolean flag that denotes if a ligand is occuring with an alternative charge. For some of the ligands, the charge determination algorithm gave multiple different charges for different occurrences of the same ligand (same topology and connection atoms). Usually one charge was in the majority and all others were discarded as outliers/errors. However, in cases where another charge was present in a significant amount (>25%) it was also recorded with the flag `is_alternative_charge` set to `True`.

###### [ligands_descriptors.csv](ligands_descriptors.csv)
- The calculated RDKit, steric and electronic descriptors for all ligands as described in the publication.
- Some properties were specifically calculated based on either the geometry of the most stable occurrence or for the gas phase optimized structure. The column names reflect this with the prefix `L*` to refer to the most stable occurrence and the prefix `L_free` to refer to the relaxed structure. Properties without a prefix were simply calculated based on the ligands SMILES string.

###### [weisfeiler_lehman_graph_hashes.csv](weisfeiler_lehman_graph_hashes.csv)
- List of Weisfeiler-Lehman graph hashes for all ligands.
- The column `base_hash` refers to the hashes only considering connectivity.
- The columns `atom_attribution_hash`, `bond_attribution_hash`, and `atom_bond_attribution_hash` refer to the hashes also considering atom element, bond order, and both atom element as well as bond order, respectively.

###### [benchmarks/OctLig_tmQMg-L_charge_agreement.csv](benchmarks/OctLig_tmQMg-L_charge_agreement.csv)
- A list of ligands included in the OctLig (Kulik and co-workers, DOI: 10.1021/acs.jctc.2c00468) and tmQMg-L datasets and their determined charges. The columns `OctLig_name` and `tmQMg-L_name` provide the identifiers from both the OctLig and tmQMg-L datasets, respectively, and the columns `OctLig_charge` and `tmQMg-L_charge` provide their corresponding predicted charges. The column `charge_agreement` contains TRUE if the two charges agree and FALSE otherwise. The last column provides the SMILES string for each ligand. Graph matching was done by generating cutoff radius graphs for all ligands of both datasets and performing node attributed graph isomorphy tests to determine equivalent ligands in the two datasets. If a ligand is only found in the OctLig dataset but not in tmQMg-L, only the columns `OctLig_name` and `smiles` will contain entries.

###### [benchmarks/1.37M_space_ground_truth.csv](benchmarks/1.37M_space_ground_truth.csv)
- The ground truth values for Polarizability and HOMO-LUMO gap calculated with GFN2-xTB of a chemical space of square-planar palladium (II) complexes with 50 monodentate ligands (25 neutral, 25 mono-anionic). This space was used to benchmark the [PL-MOGA](https://github.com/hkneiding/PL-MOGA) approach presented in the associated paper.

###### [xyz/](xyz/)
- Directory containing the geometries of all ligands ([xyz/ligands_xyzs.xyz](xyz/ligands_xyzs.xyz)), only the stable ligands ([xyz/ligands_stable_xyzs.xyz](xyz/ligands_stable_xyzs.xyz)) and the optimized stable ligands ([xyz/ligands_stable_xyzs_opt.xyz](xyz/ligands_stable_xyzs_opt.xyz)).
- With Python the xyzs can easily be loaded as a dictionary with the occurrence names as keys and the xyzs as values using the following code snippet:
```
xyzs = {}
with open('./xyz/ligands_xyzs.xyz)', 'r') as fh:
	for xyz in fh.read().split('\n\n'):
		xyzs[xyz.split('\n')[1]] = xyz
```

###### [stages/](stages/)
- Directory containing subdirectories with the scripts for each processing step in the ligand extraction pipeline.

---

[![CC BY NC 4.0][cc-by-nc-image]][cc-by-nc]

This work is licensed under a
[Creative Commons Attribution-NonCommercial 4.0 International License][cc-by-nc].

[cc-by-nc]: http://creativecommons.org/licenses/by-nc/4.0/
[cc-by-nc-image]: https://i.creativecommons.org/l/by-nc/4.0/88x31.png
