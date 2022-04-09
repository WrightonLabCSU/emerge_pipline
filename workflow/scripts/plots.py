import altair as alt
import pandas as pd


product_refied = pd.read_csv('results/EMERGE_032822_product_refined.tsv', 
                             sep='\t', index_col=0).T
product_refiedold = pd.read_csv('results_old_kegg/EMERGE_012722_product_refined.tsv', 
                             sep='\t', index_col=0).T

product = pd.read_csv('results/EMERGE_032822_product.tsv', 
                             sep='\t', index_col=0).T
distillate= pd.read_csv('results/EMERGE_032822_distillate.tsv', 
                             sep='\t', index_col=0).T

anammox = {
    'nitrogen_key':  ['K20934'],
    'nitrogen_redox-8': ['K20932','K20933','K20934'],
    'nitrogen_redox-9': ['K20935']}
distillate[set([j for i in anammox.values() for j in i])].T.sum(axis=1)

annotations = pd.read_csv('results/all_annotations.tsv', 
                             sep='\t', index_col=0).T
annotationsold = pd.read_csv('results_old_kegg/all_annotations.tsv', 
                             sep='\t', index_col=0).T
annotations = annotations.T
annotationsold = annotationsold.T
kos = annotations[['ko_id']]
koso = annotationsold[['kegg_id']]
all_kos = pd.merge(kos, koso, how='left', right_index=True, left_index=True)

all_kos.head()
annotations.head()
len(all_kos) - 
(all_kos.head()['kegg_id'].equals(all_kos.head()['ko_id']))
sum(all_kos['kegg_id'].diff(all_kos['ko_id']))
all_kos['kegg_id'].difference(all_kos['ko_id'])
all_kos.apply(lambda x: x['kegg_id'] is None, axis=1)
annotations.head()
for i 
roduct_refied.shape
product_refiedold.shape
product_refied_sum = pd.DataFrame(product_refied.apply(sum, axis=1), columns=['sum']).reset_index()
product_refied_sum.rename(columns={"index":"reaction"}, inplace=True)
print(product_refied_sum)
source = product_refied_sum

bars = alt.Chart(source).mark_bar().encode(
    x='sum',
    y="reaction"
).properties(
    width=700,
    height=900
)


text = bars.mark_text(
    align='left',
    baseline='middle',
    dx=40  # Nudges text to right so it doesn't appear on top of the bar
).encode(
    text='sum'
)

bars.save('results/product_refined_prevalance_bars.html')

