import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import pickle as pk

def dotplot_gseapy(df, top_n=10, size_col='Overlap', pval_col='FDR'):
    """
    Generate a dot plot similar to gseapy.dotplot() without using gseapy.

    Parameters:
        df (pd.DataFrame): DataFrame containing enrichment results with at least:
            - 'Term': Enriched pathway names
            - 'Adjusted P-value' (or similar): Adjusted p-values for enrichment
            - 'GeneRatio': Ratio of enriched genes (or an alternative like 'Overlap')
        top_n (int): Number of top enriched terms to display (sorted by p-value)
        size_col (str): Column representing the gene ratio (or number of genes)
        pval_col (str): Column representing the adjusted p-value (or FDR)
    """
    # Sort and take the top N terms
    df = df.sort_values(by=pval_col).head(top_n)

    # Convert p-values to -log10 scale for better visualization
    df['log_pval'] = -np.log10(df[pval_col])

    # Ensure GeneRatio is numeric (or convert if needed)
    if isinstance(df[size_col].iloc[0], str) and '/' in df[size_col].iloc[0]:
        df[size_col] = df[size_col].apply(lambda x: eval(x))  # Convert '5/100' to float

    # Create the dot plot
    plt.figure(figsize=(8, 6))
    scatter = sns.scatterplot(
        data=df,
        x='log_pval',
        y='Term',
        size=size_col,
        hue='Combined Score',
        palette='coolwarm',
        sizes=(20, 200),
        edgecolor='black'
    )

    # Labels and formatting
    plt.xlabel('-log10 Adjusted P-value')
    plt.ylabel('Enriched Terms')
    plt.title('Enrichment Analysis Dot Plot')
    plt.legend(title='', loc='best')
    plt.grid(True, linestyle='--', alpha=0.5)

    # Show the plot
    plt.show()

#Example
with open('output/G_elma_icd_male_genetic_fdr_dict.pickle', 'rb') as handle:
    G_elma_icd_male_genetic_fdr_dict = pk.load(handle)

dotplot_gseapy(G_elma_icd_male_genetic_enrichment_dict[('K35','K37')])
