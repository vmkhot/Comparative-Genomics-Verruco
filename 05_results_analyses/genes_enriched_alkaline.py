import os
import sys
import glob, csv
import pandas as pd
import numpy as np
from pathlib import Path

alkaline_node_list = [112, 79, 113, 80, 129, 93, 97, 136, 123, 133, "g0009", "g0011", "g0012", "g0013", "g0014", "g0015", "g0016", "g0017", "g0029", "g0030", "g0031", "g0032", "g0033", "g0034", "g0035", 'g0036']

df_full = pd.read_csv('Events_at_all_nodes_above_0.5.csv',header=0)
df_copies = df_full[['Node','Gene Family','Copies']]
df_copies_filt = df_copies.loc[df_copies['Copies'] >= 0.5]
df_copies_filt['alkaline_or_not'] = np.where(np.isin(df_copies_filt['Node'],alkaline_node_list), "alkaline", "not_alkaline")
    
# print(df_copies_filt)

df_sum = df_copies_filt.groupby(['Gene Family', 'alkaline_or_not']).sum().reset_index()
# print(df_sum)

df_sum_pivot = df_sum.pivot(index='Gene Family', columns='alkaline_or_not', values='Copies')

def add_alkaline_presence(row):
    if (row['alkaline'] > 0) & (row['not_alkaline'] > 0):
        result = "both"
    elif (np.isnan(row['alkaline'])) & (row['not_alkaline'] > 0):
        result = "non-alkaline-unique"
    elif (row['alkaline'] > 0) & (np.isnan(row['not_alkaline'])):
        result = "alkaline-unique"
    return result

df_sum_pivot['present_in_multiple'] =df_sum_pivot.apply(lambda row: add_alkaline_presence(row), axis=1)
# print(df_sum_pivot)

df_sum_pivot['alkaline_enriched'] = df_sum_pivot['alkaline'] / df_sum_pivot['not_alkaline']
# print(df_sum_pivot)

df_sum_pivot['alkaline_enriched'] = np.where(df_sum_pivot['present_in_multiple'] == "alkaline-unique", 100,df_sum_pivot['alkaline_enriched'])
df_sum_pivot['alkaline_enriched'] = np.where(df_sum_pivot['present_in_multiple'] == "non-alkaline-unique", 0,df_sum_pivot['alkaline_enriched'])
df_sum_pivot.reset_index(inplace=True)

df_sum_pivot[['cluster_id']] = df_sum_pivot['Gene Family'].str.split("_", expand=True)[1]
df_sum_pivot.set_index('cluster_id')
# print(df_sum_pivot)


df_annotations = pd.read_csv('homologues.tsv', sep='\t').set_index('cluster id')
# print(df_annotations)

df_merged = df_sum_pivot.join(df_annotations, how='left')
# print(df_merged)

df_merged = df_merged[['Gene Family','alkaline','not_alkaline','present_in_multiple','alkaline_enriched','cluster_id','type','annotation', 'length','selective pressure (low = purifying)','codon bias',"representation", 'count','% id']]
df_merged.rename({'Gene Family':'Gene_Family','alkaline':'alkaline_gene_count','not_alkaline':'not_alkaline_gene_count','selective pressure (low = purifying)':'selective_pressure','codon bias':'codon_bias','% id':'pid'},inplace=True)

print(df_merged)

# df_merged_sorted = df_merged.sort_values(['alkaline_enriched','alkaline'], ascending=[False,False], inplace=True)

# print(df_merged.head(n=100))