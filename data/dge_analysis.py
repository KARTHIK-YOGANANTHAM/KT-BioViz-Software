import pandas as pd
import numpy as np
import seaborn as sns
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import streamlit as st

# Function to perform DGE analysis
def perform_dge(expression_data, control_cols, disease_cols):
    """ perform Differential Gene Expression Analysis """
    results = []

    for gene in expression_data.index:
        control_values = expression_data.loc[gene, control_cols]
        disease_values = expression_data.loc[gene, disease_cols]

        # Unpacking three values correctly
        t_stat, p_value, _ = sm.stats.ttest_ind(control_values, disease_values, alternative="two-sided")

        # compute log fold change
        mean_control = np.mean(control_values)
        mean_disease = np.mean(disease_values)

        if mean_control == 0:
            log_fc = np.log2(mean_disease / mean_control)
        else:
            log_fc = np.log2(mean_disease / mean_control)
        
        results.append([gene, log_fc, p_value])

    # convert results to DataFrame
    df_results = pd.DataFrame(results, columns=["Gene", "LogFC", "P_value"])

   # Adjust P-values using Benjamini-Hochberg FDR correction
    df_results["Adj_P_value"] = multipletests(df_results["P_value"], method="fdr_bh")[1]
    
    return df_results


# Function to generate a Volcano plot
def plot_volcano(dge_results):
    """ Generates a Volcano Plot from DGE results """
    plt.figure(figsize=(8, 6))
    
    # Define thresholds for significance
    threshold = 0.05  # Adjusted P-value threshold
    logfc_cutoff = 1  # Log2 fold change threshold

    # Define colors
    dge_results["color"] = "grey"
    dge_results.loc[(dge_results["P_value"] < threshold) & (dge_results["LogFC"] > logfc_cutoff), "color"] = "red"
    dge_results.loc[(dge_results["P_value"] < threshold) & (dge_results["LogFC"] < -logfc_cutoff), "color"] = "blue"

    # Create scatter plot
    sns.scatterplot(x="LogFC", y=-np.log10(dge_results["P_value"]), hue=dge_results["color"], 
                    palette={"red": "red", "blue": "blue", "grey": "grey"}, edgecolor="black", alpha=0.7, data=dge_results)

    plt.axhline(-np.log10(threshold), linestyle="--", color="black", linewidth=1)
    plt.axvline(logfc_cutoff, linestyle="--", color="black", linewidth=1)
    plt.axvline(-logfc_cutoff, linestyle="--", color="black", linewidth=1)
    
    plt.xlabel("Log2 Fold Change")
    plt.ylabel("-Log10 Adjusted P-value")
    plt.title("Volcano Plot")
    plt.legend(["Significance Threshold", "Upregulated", "Downregulated"], loc="upper right")
    st.pyplot(plt)


# Function to generate a Heatmap
# Function to generate a Heatmap
def plot_heatmap(expression_data, dge_results):
    """ Generates a Heatmap for significantly different genes """
    
    # Select significant genes (P < 0.05)
    significant_genes = dge_results[dge_results["P_value"] < 0.05]["Gene"].tolist()
    
    if len(significant_genes) == 0:
        st.warning("No significant genes found for heatmap. Try adjusting the threshold.")
        return

    # Extract expression values for significant genes
    subset_data = expression_data.loc[significant_genes]
    subset_data = subset_data.apply(lambda x: (x-x.mean())/x.std(), axis=1)

    # Plot heatmap
    plt.figure(figsize=(12, 10))
    ax=sns.heatmap(subset_data, cmap="RdBu_r", xticklabels=True, yticklabels=True, linewidths=0.5, cbar_kws={"shrink":0.6})
    plt.xticks(rotation = 45, fontsize= 10, ha="right")
    plt.yticks(fontsize=8)
    plt.title("Heatmap of Significant Genes")
    st.pyplot(plt)


# Function to perform PCA Analysis
def perform_pca(expression_data):
    """ Performs PCA on the expression data """
    pca = PCA(n_components=2)
    transformed_data = pca.fit_transform(expression_data.T)  # Transpose to make samples rows

    pca_df = pd.DataFrame(transformed_data, columns=["PC1", "PC2"])
    pca_df["Sample"] = expression_data.columns

    # Plot PCA
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x="PC1", y="PC2", data=pca_df, hue="Sample", palette="Set1", s=100)
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.2f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.2f}%)")
    plt.title("PCA Plot")
    st.pyplot(plt)
