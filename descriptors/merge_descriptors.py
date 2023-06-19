import pandas as pd
from functools import reduce


# This script merges the four separate descriptor files.
# Namely these are: RDKit, steric, and SP and OPT electronic
# descriptors.

# load individual descriptor files
rdkit_df = pd.read_csv('./rdkit_descriptors.csv', sep=';')
steric_df = pd.read_csv('./steric_descriptors.csv', sep=';')
sp_df = pd.read_csv('./electronic_descriptors_sp.csv', sep=';').drop('occurrence_name', axis=1).add_prefix('L*-')
opt_df = pd.read_csv('./electronic_descriptors_opt+sp.csv', sep=';').drop('occurrence_name', axis=1).add_prefix('L_free-')

# merge rdkit and steric descriptors
base_df = reduce(lambda  left, right:
                     pd.merge(
                         left,
                         right,
                         on=['name'],
                         how='inner'
                        ), [rdkit_df, steric_df])

# merge electronic descriptors
electronic_df = reduce(lambda  left, right:
                           pd.merge(
                               left,
                               right,
                               left_on='L*-name',
                               right_on='L_free-name',
                               how='inner'
                            ), [sp_df, opt_df])

full_df = reduce(lambda  left,right: pd.merge(left, right, left_on='name', right_on='L*-name', how='inner'), [base_df, electronic_df])

# remove superfluous columns
full_df = full_df.drop(['L*-name', 'L_free-name'], axis=1)

# reorder
col_names = full_df.columns.to_list()
del col_names[col_names.index('stable_occurrence_name')]
col_names.insert(1, 'stable_occurrence_name')
full_df = full_df[col_names]

full_df.to_csv('./descriptors.csv', sep=';', index=False)
