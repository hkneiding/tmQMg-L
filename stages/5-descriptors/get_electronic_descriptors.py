import pandas as pd


# load dfs
stable_ligands = pd.read_csv('../3-stable/stable.csv', sep=';')
sp_df = pd.read_csv('../2-singlepoints/sp_summary.csv', sep=';').add_prefix('L*-')
opt_sp_df = pd.read_csv('../4-optimizations/opt+sp_summary.csv', sep=';').add_prefix('L_free-')

# merge electronic descriptors
electronic_df = pd.merge(
                    sp_df,
                    opt_sp_df,
                    left_on='L*-occurrence',
                    right_on='L_free-occurrence',
                    how='inner'
)

# merge with stable ligands
electronic_df = pd.merge(
                    stable_ligands,
                    electronic_df,
                    left_on='stable_occurrence',
                    right_on='L*-occurrence',
                    how='inner'
).drop(['L*-occurrence', 'L_free-occurrence', 'stable_occurrence'], axis=1)

electronic_df.to_csv('electronic_descriptors.csv', index=False, sep=';')
