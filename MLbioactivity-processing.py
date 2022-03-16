import pandas as pd
from chembl_webresource_client.new_client import new_client
target = new_client.target
target_query = target.search('CYP17')
targets = pd.DataFrame.from_dict(target_query)
# print(targets)
selected_target = targets.target_chembl_id[1]
# print(selected_target)

activity = new_client.activity
res = activity.filter(target_chembl_id=selected_target).filter(standard_type='IC50')
df = pd.DataFrame.from_dict(res)
df.standard_type.unique()
#df.to_csv('bioactivitydata.csv', index=False)
df2 = df.dropna(subset=['standard_value', 'canonical_smiles'])
df3 = df2.drop_duplicates(['canonical_smiles'])

bioactivity_class = []
for i in df3.standard_value:
    if float(i) >= 10000:
        bioactivity_class.append('inactive')
    elif float(i) <= 1000:
        bioactivity_class.append('active')
    else:
        bioactivity_class.append('intermediate')
print(bioactivity_class)
print(len(bioactivity_class))
selection = ['molecule_chembl_id', 'canonical_smiles', 'standard_value']
df4 = df3[selection]
df4 = df4.assign(bioactivityclass=bioactivity_class)
print(df4)
df4.to_csv('bioactivitydata_preprocessed.csv', index=False)







