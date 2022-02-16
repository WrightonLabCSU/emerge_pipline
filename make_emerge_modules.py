import pandas as pd
import re

rawpaths = pd.read_excel('./from_emerge_team/metabolic_reactions_updatedwithDRAMinfo.ods', sheet_name='metabolic_reactions_updatedwith', engine='odf')
ko_ids = set([j for i in rawpaths['definition'].values for j in re.findall('K\d\d\d\d\d', i)])
# rawdram = pd.read_excel('./metabolic_reactions_updatedwithDRAMinfo.ods', sheet_name='DRAM_distillate', engine='odf')
# dram_kos = set(rawdram['gene_id'][rawdram['gene_id'].str.startswith('K')])

kopaths = pd.DataFrame({j:i for _, i in rawpaths.iterrows() for j in re.findall('K\d\d\d\d\d', i['definition'])}).T
#Take the fist instance where the ko is in the path
kopaths = kopaths.groupby(by=kopaths.index).first()

emerge_distilate = pd.DataFrame({
    "gene_id":kopaths.index,
    "gene_description": ['EMERGE' for i in kopaths.index],
    "module": kopaths["pathway"].astype(str) + '_' + kopaths["reaction"].astype(str),
    "sheet": ['EMERGE' for i in kopaths.index],
    "header": ['Emerge' for i in kopaths.index],
    "subheader": kopaths["substrate_names"] + '->' + kopaths["product_names"],
    "potential_amg": [False for i in kopaths.index]
})

emerge_distilate.to_csv("./custom_input_modules/EMERGE_distillate_module.tsv", sep='\t', index=False)
emerge_paths = pd.DataFrame({
    "module_id": ['EMERGE' + str(i) for i in rawpaths.index],
    "module_name": rawpaths["substrate_names"] + '->' + rawpaths["product_names"] + ' -' + rawpaths['reaction'].astype(str),
    "complex": rawpaths["pathway"],
    "definition": rawpaths['definition']
})

len(emerge_paths['module_id'])
len(emerge_paths['module_id'].unique())
emerge_paths.to_csv("./custom_input_modules/EMERGE_pathways_module.tsv", sep='\t', index=False)
