import pandas as pd 
from datetime import datetime, date
import os

print(datetime.strptime(, '%y-%m-%d')).date()

xl = pd.ExcelFile('./emerge_distilled_output/metabolism_summary.xlsx')
emerge = xl.parse('EMERGE')
prod = pd.read_csv('./emerge_distilled_output/product.tsv', sep='\t')
path = pd.read_csv('./custom_input_modules/EMERGE_pathways_module.tsv', sep='\t')
cols = ['genome'] + [f"Complex {i}: {j}" for i, j in 
                     zip(path['complex'].values, path['module_name'].values)]

assert len(set(cols) - set([i for i in prod.columns])) == 0, "Your column names are wrong"

# for i in energy.columns[0:10]:
#     print(i, sum(energy[i].apply(lambda x: x in ect_mods)))

os.system(f"git mv"
          f"./emerge_filtered_output/EMERGE_*_distillate.tsv"
          f" ./emerge_filtered_output/EMERGE_{date.today().strftime('%d%m%y')}_distillate.tsv")
emerge.to_csv(f"./emerge_filtered_output/EMERGE_{date.today().strftime('%d%m%y')}_distillate.tsv", sep='\t', index=False)
os.system(f"git mv"
          f"./emerge_filtered_output/EMERGE_*_product.tsv"
          f" ./emerge_filtered_output/EMERGE_{date.today().strftime('%d%m%y')}_product.tsv")
prod[cols].to_csv(f"./emerge_filtered_output/EMERGE_{date.today().strftime('%d%m%y')}_product.tsv", sep='\t', index=False)
