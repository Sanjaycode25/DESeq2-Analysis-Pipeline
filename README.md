
# DESeq2-Based Transcriptomic Analysis (POLK_KO vs Control)

**📌 Project Overview**

This project demonstrates a differential gene expression analysis workflow using DESeq2 in R. It includes raw RNA-seq counts, statistical analysis, and visualization of results through MA plots, PCA plots, Volcano plots, and Heatmaps.
The dataset used in this analysis comes from GSE145313 and compares Control vs POLK_KO samples.


## 📂 Repository Structure

DESeq2-Analysis-Pipeline/
│
├── Dataset/       # Contains the raw counts dataset (TSV format)
├── Script/        # Contains the R script for DESeq2 analysis
├── Results/       # Contains CSV files of DEGs and annotated results
└── Plots/         # Contains PCA, Volcano, Heatmap, and MA plots

## 🛠️ Tools & Libraries Used
**R Packages:**

DESeq2 → Differential expression analysis

ggplot2 → Visualization

pheatmap → Heatmap generation

EnhancedVolcano → Volcano plots

biomaRt & org.Hs.eg.db → Gene annotation

Dataset: GSE145313 raw counts data (GRCh38.p13)
## ✅ Steps in the Workflow
**1.Load Raw Counts Data**

--Import the count matrix from the dataset.

**2.Prepare Metadata**

--Define sample conditions (Control vs POLK_KO).

**3.Create DESeq2 Dataset**

--Build DESeqDataSetFromMatrix using count data and metadata.

**4.Filter Low Counts**

--Remove genes with very low expression.

**5.Run DESeq Analysis**

--Perform normalization and differential expression analysis.

**6.Extract Results**

--Get log2 fold changes, p-values, adjusted p-values.

**7.Generate Plots**

--MA plot, PCA plot, Volcano plot, Heatmap of top 30 genes.

**8.Annotate Genes**

--Add gene symbols and descriptions using Ensembl/OrgDb.

**9.Save Outputs**

--Export results as CSV and plots as PDF.

#**📊 Visualizations Included**

MA Plot: Shows log2 fold changes vs mean expression.

PCA Plot: Displays sample clustering based on variance.

Volcano Plot: Highlights significant up/down-regulated genes.

Heatmap: Top 30 DEGs with hierarchical clustering.
## 🚀 How to Run the Analysis

**1.Clone the repository:**

git clone https://github.com/Sanjaycode25/DESeq2-Analysis-Pipeline.git


**2.Open the R script:**

"Script/deseq2_analysis.R"

**3.Install required R packages:**

"install.packages(c("pheatmap"))
BiocManager::install(c("DESeq2", "EnhancedVolcano", "biomaRt", "org.Hs.eg.db"))"


**4.Update file paths in the script if needed.**

**5.Run the script in R or RStudio.**
## 📁 Outputs Generated

**Results:**

DESeq2_results_raw.csv → Raw DESeq2 output.

DESeq2_results_annotated.csv → Annotated results with gene symbols.

DEG_POLK_KO_vs_Control.csv → Significant DEGs.

**Plots:**

PCA plot, MA plot, Volcano plot, Heatmap.
## 📜 License

This project is licensed under the MIT License – free to use and modify.
## 🤝 Contributions

Feel free to fork the repo, open issues, or submit pull requests to improve the analysis.