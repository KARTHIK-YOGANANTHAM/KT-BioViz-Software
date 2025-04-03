import streamlit as st
import pandas as pd
from dge_analysis import perform_dge, plot_volcano, plot_heatmap, perform_pca

# --- APP SETTINGS ---
st.set_page_config(page_title="KT-BioViz🧬", page_icon="🧬", layout="wide")

# --- HEADER ---
st.image("logo.png", width=250)  # Add your logo file
st.title("KT-BioViz🧬: Advanced Differential Gene Expression Analysis Tool🔎")
st.markdown("""
### Welcome to Key Taxanomy-Biological Vizualization Tool (KT-BioViz)🧬
KT-BioViz is an interactive tool for **Differential Gene Expression (DGE) Analysis**, allowing you to:
- Upload expression data🔎🔽
- Identify significantly expressed genes🧑‍💻
- Generate **Volcano Plots🌋, Heatmaps🔥, and PCA analysis📈**
- Export results for downstream analysis📒

👉 **Follow the steps below to perform your analysis!**
""")

# --- SIDEBAR ---
st.sidebar.image("logo.png", width=100)  # Add your sidebar icon
st.sidebar.title("Navigation")
page = st.sidebar.radio("Select Analysis Step", ["Upload Data", "Run DGE Analysis", "Visualizations", "Documentation"])

# --- UPLOAD DATA ---
if page == "Upload Data":
    st.header("Upload Your Expression Data👨‍🔬")
    st.subheader("Upload a CSV file with Gene Expression Data⏬")

    uploaded_file = st.file_uploader("Upload CSV file", type=["csv"])
    if uploaded_file is not None:
        df = pd.read_csv(uploaded_file, index_col=0)  # Gene names as index
        st.write(df.head())
        st.session_state["expression_data"] = df
        st.success("Data uploaded successfully!")

# --- RUN DGE ANALYSIS ---
elif page == "Run DGE Analysis":
    st.header("Run Differential Gene Expression Analysis🕵️‍♂️")
    st.write("The DGE analysis in KT-BioViz is based on t-test and log2 Fold change between expression values")
    st.write("- **Method** : Uses t-test + log2 fold change")
    st.write("- **Similar to** : Basic limma analysis (without Bayesian shrinkage)")
    st.write("- **Limitations** : Doesn't model RNA-seq count distribution (assumes normality, which RNA-seq data lacks)")
    
    if "expression_data" in st.session_state:
        df = st.session_state["expression_data"]
        columns = list(df.columns)
        control_cols = st.multiselect("Select Control Samples", columns)
        disease_cols = st.multiselect("Select Disease Samples", columns)

        if st.button("Run DGE Analysis"):
            if not control_cols or not disease_cols:
                st.warning("Please select both control and disease samples!")
            else:
                dge_results = perform_dge(df, control_cols, disease_cols)
                st.session_state["dge_results"] = dge_results
                st.write("### DGE Results", dge_results.head())
                dge_results.to_csv("dge_results.csv", index=False)
                st.success("DGE Analysis complete! Proceed to visualizations.")
    else:
        st.warning("Please upload data first in the 'Upload Data' section.")

# --- VISUALIZATIONS ---
elif page == "Visualizations":
    st.header("Generate Visualizations")
    st.write("------------------------------------------------------------------------------------------------------------------------------------------")
    
    if "dge_results" in st.session_state:
        st.subheader("Create Volcano Plots🌋")
        st.write("- A scatter plot used to visualize differentially expressed genes based on log2 fold change (x-axis) and p-value (-log10 scale, y-axis).")
        st.write("- Highlights significantly upregulated (red) and downregulated (blue) genes, aiding in biomarker discovery.")
        st.write("- Helps quickly identify genes with both statistical significance and biological relevance. ")
        st.write("- Common threshold: **|log2FC| > 1** and **p-value < 0.05** (or adjusted p-value like FDR).")

        if st.button("Generate Volcano Plot"):
            plot_volcano(st.session_state["dge_results"])
        
        st.write("------------------------------------------------------------------------------------------------------------------------------------------")
        st.subheader("Create Heatmap🔥")
        st.markdown(""" 
- A color-coded matrix that visualizes gene expression patterns across multiple samples.

- Rows represent genes, columns represent samples, and colors indicate expression levels (red = high, blue = low).

- Often includes hierarchical clustering to group similar gene expression profiles.

- Useful for detecting expression trends, sample similarities, and potential biomarkers.""")
        st.write("SELECT THRESHOLD FOR HEATMAP⏬")
        p_value_threshold = st.slider("Select P-value Threshold for Heatmap", 0.01, 0.1, 0.05, 0.01)
        if st.button("Generate Heatmap"):
            filtered_dge_results = st.session_state["dge_results"][st.session_state["dge_results"]["P_value"] < p_value_threshold]
            if filtered_dge_results.empty:
                st.warning("No significant genes found for the selected threshold. Try adjusting the slider.")
            else:
                plot_heatmap(st.session_state["expression_data"], filtered_dge_results)
        
        st.write("------------------------------------------------------------------------------------------------------------------------------------------")
        st.subheader("Perform PCA Analysis📈")
        st.markdown(""" 
- A dimensionality reduction technique that identifies major variations in gene expression across samples.

- Converts high-dimensional gene expression data into principal components (PC1, PC2, etc.) to reveal patterns.

- Clusters similar samples together, distinguishing experimental conditions (e.g., disease vs. control).

- Helps detect outliers and batch effects, improving downstream analysis.""" )
        
        if st.button("Perform PCA Analysis"):
            perform_pca(st.session_state["expression_data"])
    else:
        st.warning("Run DGE analysis first!")

# --- DOCUMENTATION ---
elif page == "Documentation":
    st.header("Documentation & How-to Guide")
    st.markdown("""
   # KT-BioViz Documentation

## Table of Contents
- [Introduction](#introduction)
- [Options](#options)
- [Applications](#applications)
- [About the Developer](#about-the-developer)
- [Contact Us](#contact-us)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [License](#license)

---

## Introduction
KT-BioViz is an interactive gene expression visualization tool designed to facilitate biological data analysis. It enables researchers and bioinformaticians to explore and interpret gene expression data efficiently with user-friendly visualization features.

## Options
KT-BioViz provides the following options for seamless data exploration:

- **Data Upload**: Upload gene expression datasets in CSV, TSV, or other supported formats.
- **Filtering & Normalization**: Apply data filtering techniques and normalize expression values.
- **Visualization Tools**:
  - Volcano plots
  - Heatmaps
  - PCA (Principal Component Analysis)
- **Statistical Analysis**: Perform differential expression analysis, clustering.
- **Export & Sharing**: Export plots and results in multiple formats (PNG, PDF, CSV).
- **Customization**: Select control and disease columns as you preferred.

## Applications
KT-BioViz is useful in various bioinformatics and biotechnology applications, including:

- **Gene Expression Analysis**: Visualizing expression patterns across different conditions.
- **Disease Biomarker Identification**: Identifying differentially expressed genes linked to diseases.
- **Drug Discovery & Target Validation**: Exploring gene expression responses to treatments.
- **Pathway & Functional Analysis**: Understanding gene functions and regulatory pathways.
- **Comparative Genomics**: Analyzing expression differences between species or conditions.

## About the Developer
This software is developed by Karthik Yoganantham, an MSc Biotechnology student specializing in bioinformatics, data analytics, and biological data visualization. With expertise in R and Python programming, their focus is on developing computational tools for biological research and advancing the field of bioinformatics.

## Contact Us
For inquiries, support, or collaboration, reach out to us:

- **Email**: karthikyoganantham@gmail.com
- **GitHub**: [github.com/KARTHIK-YOGNANTHAM](https://github.com/KARTHIK-YOGANANTHAM)
- **LinkedIn**: [linkedin.com/in/karthik-yoganantham](www.linkedin.com/in/karthik-yoganantham)

## System Requirements
To run KT-BioViz efficiently, ensure your system meets the following requirements:

- **Operating System**: Windows, macOS, or Linux
- **Processor**: Intel Core i3 (or higher)
- **RAM**: Minimum 8GB (Recommended 16GB for large datasets)
- **Storage**: At least 5GB of free space
- **Software Dependencies**:
  - Python (version X.X or later)
  - R (if applicable)
  - Required libraries: pandas, matplotlib, seaborn, plotly, etc.


## License
KT-BioViz is licensed under the [MIT License](LICENSE). You are free to use, modify, and distribute this software following the license terms.

    """)
