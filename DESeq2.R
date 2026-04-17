# RNASeq project

#1 COLDATA
setwd("~/Desktop/DESeq2")
#a) finding all sample names
names <- list.files(path = "./quant")
#b) creating paths to the cont files "quant.sf"
files <- file.path("./quant", names, "quant.sf")
#c) constructing coldata
coldata <- data.frame(names,files)
#d) Adding Condition to the coldata
coldata$condition <- factor(ifelse(grepl("N$", coldata$names), "Normal", "Tumor"))
#e) Adding the patientID
coldata$patientID <- substr(coldata$names, 1, nchar(coldata$names)-1)
coldata$patientID <- factor(coldata$patientID)

#2 IMPORT QUANT FILES INTO R
#a) load library
library("tximeta")
#b)  Read the data into a SummarizedExperiment object
se <- tximeta::tximeta(coldata = coldata, type = "salmon")
#c) Summarize transcript-level quantifications to the gene level.
gse <- tximeta::summarizeToGene(se)

#3 CONSTRUCT DESEQ2 DATASET
# a) Set the control group
gse$condition <- stats::relevel(factor(gse$condition), ref = "Normal")
# b) Create the DESeq2 object with the design formula
library("DESeq2")
dds <- DESeq2::DESeqDataSet(se = gse, design = ~ patientID + condition)

#4 DATA FILTERING
# a) Filter for at least 10 counts in 
keep <- rowSums(counts(dds) >= 10) >= 4
# b) Apply the filter
dds <- dds[keep, ]

#5 EXploratory Data Analysis
# a) Transform counts to stabilize variance for Visualization
vsd <- DESeq2::vst(object = dds, blind = FALSE )
# b) Generate PCA plot
DESeq2::plotPCA(object = vsd, intgroup = "condition")

#6 Differential Expression Testing
# a) Run the DESeq pipeline
dds2 <- DESeq2::DESeq(object = dds)
# b) Extract the results table
res <- DESeq2::results(object = dds2)
# c) Summary of the analysis
summary(res)

# 7) Annotation
# a) Load library for annotation
library(AnnotationDbi)
library(org.Hs.eg.db) 
# b) Annotate with symbols
res$symbol <- AnnotationDbi::mapIds(x = org.Hs.eg.db, keys = row.names(res), column = "SYMBOL", keytype ="ENSEMBL" )

#8 QC plots 
# a) Dispersion plot
DESeq2::plotDispEsts(object = dds2, legend = TRUE, legendpos = "bottomright")
# b) Sparsity plot
DESeq2::plotSparsity(x = dds2, normalized = TRUE)
# c) Histogram of pvalues
hist(x= res$pvalue)

#9 LFC Shrinkage and Visualization
# a) Shrink LFC estimates using the apeglm method
resultsNames(dds2) # see the different coefficients
resape <- DESeq2::lfcShrink(dds = dds2, coef = "condition_Tumor_vs_Normal",res = res, type = "apeglm")
resape$symbol <- res$symbol # Adding symbols from the res object
resape <- resape[!is.na(resape$symbol), ] #Remove genes without names
resape <- resape[!duplicated(resape$symbol), ] #Remove duplicated reads
# b) Subsetting for the significant genes
degs <- subset(x= resape, abs(log2FoldChange)>1 & padj< 0.05)
# c) Generate MA plot
DESeq2::plotMA(object = resape, alpha = 0.05, ylim = c(-10,10), main = "MA plot", colNonSig = "blue", colSig = "red", colLine = "black")
# d) Generate volcano plot
library(EnhancedVolcano) #load library
top30_genes <- head(resape[order(resape$padj), "symbol"], 30) # Extracting top 30 genes
EnhancedVolcano::EnhancedVolcano(toptable = resape, x = "log2FoldChange", y = "padj", title = "Volcano Plot", pCutoff = 0.05, FCcutoff = 1, lab = resape$symbol,selectLab = top30_genes,ylim = c(0, 120), drawConnectors = TRUE)

#10 Variance stabilization, pheatmap and sample to sample distance matrix
vsd2 <- DESeq2::vst(object = dds2, blind = FALSE)
vsd2_mat <- assay(vsd2)
# a) Create a heatmap of the top 30 genes
library(pheatmap)
# b) Get the IDs for the top 30 genes 
top30_ids <- rownames(resape[order(resape$padj), ][1:30, ])
# c) Subset the matrix using those IDs
vsd2_mat_top30 <- vsd2_mat[top30_ids, ]
# d) REplace the IDs with symbols
rownames(vsd2_mat_top30) <- resape[top30_ids, "symbol"]
# e) Heatmap
pheatmap::pheatmap(mat = vsd2_mat_top30, legend = TRUE,)
# f) Sample to sample distance heatmap
transposed_vsd2 <- t(vsd2_mat)
sample_dists <- dist(x= transposed_vsd2, method = "euclidean")# calclate distance from transposed matrix
sample_dists <- as.matrix(sample_dists) # transponse and change  to matrix
pheatmap::pheatmap(mat = sample_dists)

#11 Exporting Results
# a) Order results by adjusted p-value: 
degs<- degs[order(degs$padj), ]
# b) Write to a text/CSV file: 
write.table(as.data.frame(degs), file = "degs.txt", sep = "\t")

#12 Fisher's Exact Test (Over Representation Analysis; ORA)
# a) Mapping significant DEGs to Entrez IDs
degs_mapped <- clusterProfiler::bitr(geneID = degs$symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
universe <- clusterProfiler::bitr(geneID = resape$symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
# b) Kegg ORA
enrich_go <- clusterProfiler::enrichGO(gene = degs_mapped$ENTREZID, OrgDb = org.Hs.eg.db, universe = universe, readable = TRUE)
enrich_kegg <- clusterProfiler::enrichKEGG(gene = degs_mapped$ENTREZID, universe = universe)
# c) Visualization
library(enrichplot)
clusterProfiler::dotplot(enrich_go)
clusterProfiler::dotplot(enrich_kegg)

#13  Gene Set Enrichment Analysis (GSEA)
# a) Getting the dataframe of the results
res_df <- as.data.frame(x=res)
# b) Cleaning the dataframe
library(dplyr)
res_ranked <- res_df %>% 
  dplyr::filter(!is.na(symbol)) %>% 
  dplyr::filter(!is.na(stat)) %>% 
  dplyr::group_by(symbol) %>%                              # Then group...
  dplyr::slice_max(order_by = abs(stat), 
                   n = 1, 
                   with_ties = FALSE) %>%                  # Then pick the best...
  dplyr::ungroup() %>%                                     # Then ungroup 
  dplyr::arrange(desc(stat))                               # Then sort.
# c) # 1. Sort the dataframe by the Wald statistic (highest to lowest)
gene_list <- res_ranked$stat
# d) Name the vector with the symbols
names(gene_list) <- res_ranked$symbol
# e) Sort in decreasing order (Required for GSEA)
gene_list <- sort(gene_list, decreasing = TRUE)
# f) GSEA GO for Biological Process (BP)
gse_go <- clusterProfiler::gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db, keyType = "SYMBOL")
enrichplot::dotplot(gse_go)
enrichplot::gseaplot(x = gse_go, geneSetID = 1)
# g) GSEA KEGG 
id_map <- clusterProfiler::bitr(geneID= res_ranked$symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
kegg_df <- res_ranked %>% dplyr::inner_join(id_map, by = c("symbol" = "SYMBOL"))
kegg_list <- kegg_df$stat
names(kegg_list) <- kegg_df$ENTREZID
kegg_list <- sort(x = kegg_list, decreasing = TRUE)
gse_kegg <- clusterProfiler::gseKEGG(geneList = kegg_list)
enrichplot::dotplot(gse_kegg)
enrichplot::gseaplot(x = gse_kegg, geneSetID = 1)
