# Load required libraries
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)

# Step 1: Load count data
counts <- read.delim("C:/Users/skutr/Desktop/DESeq2/Dataset/GSE145313_raw_counts_GRCh38.p13_NCBI.tsv",
                     header = TRUE, row.names = 1, check.names = FALSE)

# Step 2: Create metadata
sample_info <- data.frame(
  row.names = colnames(counts),
  condition = c("Control","Control","Control","POLK_KO","POLK_KO","POLK_KO")
)

# Step 3: Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_info,
                              design = ~ condition)

# Step 4: Filter low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Step 5: Run DESeq
dds <- DESeq(dds)

# Step 6: Get results
res <- results(dds, contrast = c("condition", "POLK_KO", "Control"))

# Save raw DESeq2 results
write.csv(as.data.frame(res), "DESeq2_results_raw.csv", row.names = TRUE)

# Step 7: QC Plots
plotMA(res, ylim = c(-5,5))
plotPCA(vst(dds), intgroup = "condition")

# Step 8: Significant genes
resSig <- subset(res, padj < 0.05)
head(resSig[order(resSig$log2FoldChange, decreasing = TRUE),])

# Save significant DEGs
write.csv(as.data.frame(resSig), "DEG_POLK_KO_vs_Control.csv", row.names = TRUE)

# Step 9: Annotate results
resSig$symbol <- mapIds(org.Hs.eg.db,
                        keys = rownames(resSig),
                        column = "SYMBOL",
                        keytype = "ENTREZID",
                        multiVals = "first")

resSig$genename <- mapIds(org.Hs.eg.db,
                          keys = rownames(resSig),
                          column = "GENENAME",
                          keytype = "ENTREZID",
                          multiVals = "first")

resSig$ensembl <- mapIds(org.Hs.eg.db,
                         keys = rownames(resSig),
                         column = "ENSEMBL",
                         keytype = "ENTREZID",
                         multiVals = "first")

# Save annotated results
write.csv(as.data.frame(resSig), "DESeq2_results_annotated.csv", row.names = TRUE)

# Step 10: Convert to dataframe for custom plots
res_df <- as.data.frame(res)

# Volcano Plot using EnhancedVolcano
EnhancedVolcano(res_df,
                lab = rownames(res_df),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'POLK_KO vs Control',
                pCutoff = 0.05,
                FCcutoff = 2,
                pointSize = 3.0,
                labSize = 3.0)

# Step 11: Heatmap of top 30 DEGs
top_genes <- rownames(res[order(res$padj), ])[1:30]
vsd <- vst(dds, blind = FALSE)
pheatmap(assay(vsd)[top_genes, ],
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         annotation_col = sample_info["condition", drop = FALSE])

