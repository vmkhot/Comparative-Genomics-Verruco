# %%
import os
import sys
import glob, csv
import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import fisher_exact


# %%
alkaline_node_list = [79,112,123,129,133,113,80,93,97,136,138,120, "g0009", "g0011", "g0012", "g0013", "g0014", "g0015", "g0016", "g0017", "g0029", "g0030", "g0031", "g0032", "g0033", "g0034", "g0035", 'g0036']

# %%
cutoff = 0.5        # adjust for stringency

df_full = pd.read_csv('Events_at_all_nodes.csv',header=0)
df_copies = df_full[['Node','Gene Family','Copies']]
df_copies_filt = df_copies.loc[df_copies['Copies'] >= cutoff]      
df_copies_filt['alkaline_or_not'] = np.where(np.isin(df_copies_filt['Node'],alkaline_node_list), "alkaline", "not_alkaline")

# %%
df_copies_filt

# %%
# Creating a gene matrix

gene_matrix = df_copies_filt.pivot(index='Node', columns='Gene Family', values='Copies').fillna(0)
print(gene_matrix)


# %%

# convert to binary data
binary_matrix = (gene_matrix > 0).astype(int)
# print(binary_matrix)

# %%
# Creating a labels series

labels = df_copies_filt[['Node','alkaline_or_not']].drop_duplicates(['Node','alkaline_or_not']).drop('Node', axis=1)
labels = labels.squeeze()
labels.index = binary_matrix.index
print(labels)

# %% [markdown]
# ## Fisher's Exact Test
# 
# For testing significance of gene enrichment/depletion
# - Only binary data - have to turn counts to binary

# %% [markdown]
# ### Run the Fisher's test
# 
# For each gene family in a loop and then convert to dataframe
# 
# #### Interpretation of the Odds Ratio (OR)
# **OR = 1:** This indicates no association between the gene family and the genome type (alkaline vs. non-alkaline). The odds of presence in alkaline genomes are the same as the odds of presence in non-alkaline genomes.
# 
# **OR > 1:** This indicates a positive association between the gene family and alkaline genomes. The gene family is more likely to be present in alkaline genomes compared to non-alkaline genomes.
# The further the OR is greater than 1, the stronger the association. For example, an OR of 2 means that alkaline genomes are twice as likely to have the gene family present compared to non-alkaline genomes.
# 
# **OR < 1:** This indicates a negative association between the gene family and alkaline genomes, meaning the gene family is more likely to be absent in alkaline genomes compared to non-alkaline genomes.
# The further the OR is less than 1, the stronger the depletion. For example, an OR of 0.5 means that alkaline genomes are half as likely to have the gene family present compared to non-alkaline genomes.
# 
# #### Interpretation of the p-value
# p-value < 0.05 : significant
# p-value > 0.05 : non-significant
# 
# A p-value > 0.05 means that the test did not find a statistically significant difference. It does not imply that the gene family is depleted or enriched.
# 
# To determine depletion, check the odds ratio (OR). If the OR is less than 1, it suggests depletion, but the p-value must be low (typically < 0.05) for it to be statistically significant.
# 

# %%
# Initialize results
results = []

# Loop through each gene family
for gene_family in binary_matrix.columns:
    # Create contingency table
    alkaline_present = binary_matrix.loc[labels == 'alkaline', gene_family].sum()
    non_alkaline_present = binary_matrix.loc[labels == 'not_alkaline', gene_family].sum()
    alkaline_absent = len(labels[labels == 'alkaline']) - alkaline_present
    non_alkaline_absent = len(labels[labels == 'not_alkaline']) - non_alkaline_present

    contingency_table = [
        [alkaline_present, non_alkaline_present],
        [alkaline_absent, non_alkaline_absent]
    ]
    # Perform Fisher's Exact Test
    odds_ratio, p_value = fisher_exact(contingency_table, alternative='two-sided')
     # Store results
    results.append({
        'GeneFamily': gene_family,
        'OddsRatio': odds_ratio,
        'PValue': p_value
    })

# Convert results into a DataFrame for easy inspection
results_df = pd.DataFrame(results)

# Display the results
print(results_df)


# %% [markdown]
# ### ADJUSTING THE PVALUE USING BENJAMINI-HOCHBERG PROCEDURE for FDR (LESS CONSERVATIVE THAN BONFERRONI)
# 

# %%
# Assume `results_df` is the DataFrame containing the original p-values
# that you obtained from Fisher's Exact Test (with 'PValue' column)

# Step 1: Sort p-values in ascending order and assign ranks
results_df['Rank'] = results_df['PValue'].rank(method='min')

# Step 2: Calculate the adjusted p-values using the BH procedure formula
m = len(results_df)  # Total number of tests (gene families)
results_df['BH_Adjusted_PValue'] = results_df['PValue'] * m / results_df['Rank']

# Step 3: Ensure adjusted p-values do not exceed 1
results_df['BH_Adjusted_PValue'] = np.minimum(results_df['BH_Adjusted_PValue'], 1)

# Step 4: Determine which tests are significant after BH correction
results_df['Significant'] = results_df['BH_Adjusted_PValue'] < 0.05

# Display the results
print(results_df)

results_df_filtered = results_df.loc[results_df['PValue'] < 0.05]
print(results_df_filtered)

# %% [markdown]
# Dunno, I am not convinced by this as it is classifying some genes which I know are unique in alkaline as non-significant

# %% [markdown]
# ### Plotting the significant enriched/depleted genes in alkaline genomes

# %%
import matplotlib.pyplot as plt

# Assume 'results_df' contains the adjusted p-values and log-fold change (or odds ratio)
# If you don't have log-fold change, you can calculate the odds ratio and take the log2

# Example of adding log-fold change (using Odds Ratio or other measure of effect size)
# results_df['log2_fold_change'] = np.log2(results_df['OddsRatio'])

# Volcano plot
plt.figure(figsize=(8, 6))
plt.scatter(
    results_df['OddsRatio'],  # x-axis: log-fold change (enrichment or depletion)
    -np.log10(results_df['BH_Adjusted_PValue']),  # y-axis: -log10(p-value)
    c=results_df['BH_Adjusted_PValue'] < 0.05,  # Color points by significance (adjusted p-value < 0.05)
    cmap='coolwarm',  # Color map
    edgecolor='k',
    alpha=0.7
)

# Add significance threshold line
plt.axhline(y=-np.log10(0.05), color='red', linestyle='--', label='Significance Threshold')

# Labels and title
plt.title('Volcano Plot: Gene Family Enrichment/Depletion')
plt.xlabel('Odds Ratio (Enrichment/Depletion)')
plt.ylabel('-Log10 Adjusted P-value')
plt.legend()

# Show plot
plt.show()


# %% [markdown]
# ## Poisson regression
# Takes into account count data - non binary
# 

# %%
# For one gene family only

import statsmodels.api as sm

# Example: Poisson regression for count data
X = gene_matrix['pmsf_0_MSA_trimmed.ufboot.ale.uml_rec']  # Count data (number of occurrences of a gene family)
y = (labels == 'alkaline').astype(int)  # Binary outcome: Alkaline vs. Non-Alkaline


# Add constant (intercept)
X = sm.add_constant(X)


# Poisson regression model
poisson_model = sm.GLM(y, X, family=sm.families.Poisson())
poisson_results = poisson_model.fit()

print(poisson_results.summary())


# %%
import statsmodels.stats.multitest as mt

# Assuming `gene_matrix` is your dataframe of gene family counts
# and `labels` contains the corresponding 'Alkaline' or 'Non-Alkaline' labels

# Initialize a list to store p-values for each gene family
results = []

# Loop through each gene family
for gene_family in gene_matrix.columns:
    # Create the dependent variable (alkaline or not) 
    # assuming 'alkaline_or_not' is a binary column in `labels` (0 = Non-Alkaline, 1 = Alkaline)
    X = gene_matrix[gene_family]
    X = sm.add_constant(X)  # Add intercept to the model
    
    y = (labels == 'alkaline').astype(int)  # Convert labels to binary: 1 for Alkaline, 0 for Non-Alkaline
    
    # Fit Poisson regression model
    poisson_model = sm.GLM(y, X, family=sm.families.Poisson(), link=sm.families.links.log()).fit()
    
    # Get the p-value for the gene family
    results.append({
        'GeneFamily': gene_family,
        'OddsRatio': np.exp(poisson_model.params[gene_family]),
        'PValue': poisson_model.pvalues[gene_family]
    })


# %%

# Display the results
results_df_pois = pd.DataFrame(results)

# Step 1: Sort p-values in ascending order and assign ranks
results_df_pois['Rank'] = results_df_pois['PValue'].rank(method='min')

# Step 2: Calculate the adjusted p-values using the BH procedure formula
m = len(results_df_pois)  # Total number of tests (gene families)
results_df_pois['BH_Adjusted_PValue'] = results_df_pois['PValue'] * m / results_df_pois['Rank']

# Step 3: Ensure adjusted p-values do not exceed 1
results_df_pois['BH_Adjusted_PValue'] = np.minimum(results_df_pois['BH_Adjusted_PValue'], 1)

# Step 4: Determine which tests are significant after BH correction
results_df_pois['Significant'] = results_df_pois['BH_Adjusted_PValue'] < 0.05

# Display the results
print(results_df_pois)

# filtered cos it gave me some really weird odds ratios
results_df_pois_filt = results_df_pois.loc[(results_df_pois['OddsRatio'] < 50) & (results_df_pois['OddsRatio'] > 0.01)]#& ] (results_df['BH_Adjusted_PValue'] < 0.05)
print(results_df_filtered) 

# Optional: Identify significant gene families (adjusted p-value < 0.05)
# significant_gene_families = results[results['PValue'] < 0.05]
# print("Significant Gene Families:")
# print(significant_gene_families)


# %%
# Volcano plot
plt.figure(figsize=(8, 6))
plt.scatter(
    results_df_pois_filt['OddsRatio'],  # x-axis: log-fold change (enrichment or depletion)
    -np.log10(results_df_pois_filt['BH_Adjusted_PValue']),  # y-axis: -log10(p-value)
    c=results_df_pois_filt['BH_Adjusted_PValue'] < 0.05,  # Color points by significance (adjusted p-value < 0.05)
    cmap='coolwarm',  # Color map
    edgecolor='k',
    alpha=0.7
)

# Add significance threshold line
plt.axhline(y=-np.log10(0.05), color='red', linestyle='--', label='Significance Threshold')

# Labels and title
plt.title('Volcano Plot: Gene Family Enrichment/Depletion by Poisson Regression')
plt.xlabel('Odds Ratio (Enrichment/Depletion)')
plt.ylabel('-Log10 BH Adjusted P-value')
plt.legend()

# Show plot
plt.show()

# %%
# from Fisher's test

results_df_filt = results_df.loc[results_df['BH_Adjusted_PValue'] < 0.05]
results_df_filt

# from Poisson regression
results_df_pois_filt_2 = results_df_pois_filt.loc[results_df_pois_filt['BH_Adjusted_PValue'] < 0.05]
results_df_pois_filt_2

# %%
# merged the two tests

df_test_merged = pd.merge(results_df,results_df_pois_filt, how='outer', on='GeneFamily', suffixes=['_Fisher', '_Poisson'])
print(df_test_merged)
df_test_merged = df_test_merged.loc[(df_test_merged['BH_Adjusted_PValue_Fisher'] < 0.05) | (df_test_merged['BH_Adjusted_PValue_Poisson'] < 0.05)]
df_test_merged
# df_test_merged.drop(['Rank','BH_Adjusted_PValue','Significant'], inplace=True, axis=1)
df_test_merged
temp = df_test_merged['GeneFamily'].str.split("_", expand=True)[1]

# df_test_merged['cluster_id'] = temp['1']
df_test_merged = df_test_merged.join(temp)
df_test_merged.rename(columns={1 : 'cluster_id'},inplace=True)
df_test_merged.set_index('cluster_id', inplace=True)
df_test_merged.index = df_test_merged.index.astype(int)

df_test_merged



# %%
df_annotations = pd.read_csv('homologues.tsv', sep='\t')
df_annotations.rename(columns={'cluster id': 'cluster_id'},inplace=True)
df_annotations.set_index('cluster_id', inplace=True)
df_annotations

# %%
df_copies_filt_grouped = df_copies_filt.groupby('Gene Family').agg({'Node': lambda x: x.tolist()})
df_copies_filt_grouped

# %%
df_final = df_test_merged.join(df_annotations, how='left')
df_final

df_final.reset_index(inplace=True)
# print(df_final)
df_final.to_csv('genes_enriched_alkaline_statistics_incl_genomes.csv', index=False)

df_final = df_final[['GeneFamily','OddsRatio_Fisher','PValue_Fisher','BH_Adjusted_PValue_Fisher','OddsRatio_Poisson','PValue_Poisson','BH_Adjusted_PValue_Poisson','cluster_id','type','annotation', 'length','selective pressure (low = purifying)','codon bias',"representation", 'count','% id']]

df_final.rename({'selective pressure (low = purifying)':'selective_pressure','codon bias':'codon_bias','% id':'pid'},inplace=True)
df_final = df_final.merge(df_copies_filt_grouped, left_on='GeneFamily', right_on='Gene Family', how='left')
print(df_final)
df_final.to_csv('genes_enriched_alkaline_statistics_short.csv', index=False)