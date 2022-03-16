import pandas as pd
import padelpy
import wget

df3 = pd.read_csv('bioactivitydata.csv')
selection = ['canonical_smiles','molecule_chembl_id']
df3_selection = df3[selection]
wget()