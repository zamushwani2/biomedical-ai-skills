# Cancer Multi-Omics Analysis

Integrated analysis of expression, mutation, copy number, and methylation data from TCGA and GEO for solid tumor characterization.

## When to Use This Skill

Activate when the user requests:
- TCGA data download and processing for any cancer type
- Differential expression analysis between tumor conditions
- Pathway and gene set enrichment analysis
- Mutation landscape analysis (oncoplot, signatures, co-occurrence)
- Copy number variation analysis
- DNA methylation analysis (450K/EPIC arrays)
- Integration of two or more omics layers

## Inputs

| Data Type | Format | Source |
|-----------|--------|--------|
| Expression | Raw counts (STAR - Counts) | TCGA via TCGAbiolinks (v2.38+), GEO via GEOquery |
| Mutations | MAF | TCGA GDC, cBioPortal |
| Copy number | Segment files, GISTIC2.0 output | TCGA GDC, GDAC Firehose |
| Methylation | IDAT files or beta-value matrices | TCGA GDC, GEO |
| Clinical | Tabular | TCGA GDC, cBioPortal |

---

## Expression Analysis

### Data Retrieval

```r
library(TCGAbiolinks)  # v2.38.0+, Bioconductor 3.22

query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query, directory = "GDCdata")
se <- GDCprepare(query, directory = "GDCdata")
# Returns SummarizedExperiment; raw counts in assay(se, "unstranded")
# TPM in assay(se, "tpm_unstrand") — use for visualization only, never for DE

# Extract clinical data
clinical <- as.data.frame(colData(se))
```

### GEO Data Retrieval (Alternative)

```r
library(GEOquery)
gse <- getGEO("GSE72094", GSEMatrix = TRUE)[[1]]
# For count data from GEO Supplementary files:
# Download counts matrix manually, read with read.csv/read.delim
```

### Normalization Decision Tree

```
Input type?
  Raw counts (HTSeq, STAR)
    -> DE analysis: feed directly to DESeq2 (handles normalization internally)
    -> Visualization (PCA, heatmap): apply vst() or rlog() from DESeq2
    -> Cross-sample comparison: vst() preferred for n > 30 samples (faster than rlog)
  TPM/FPKM (already normalized)
    -> DE analysis: STOP. Go back, get raw counts. TPM/FPKM invalid for DESeq2/edgeR.
    -> Correlation/visualization: log2(TPM + 1), acceptable
    -> Gene set scoring (ssGSEA, GSVA): TPM acceptable as input
  RSEM expected counts (non-integer)
    -> round() before DESeq2: DESeqDataSetFromMatrix(countData = round(counts), ...)
```

### Differential Expression (DESeq2 v1.50+)

```r
library(DESeq2)  # v1.50.2+, Bioconductor 3.22
set.seed(42)

counts_mat <- assay(se, "unstranded")
col_data <- colData(se)

# Ensure consistent sample ordering
stopifnot(identical(colnames(counts_mat), rownames(col_data)))

dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData = col_data,
  design = ~ sample_type  # "Primary Tumor" vs "Solid Tissue Normal"
)

# Pre-filter: keep genes with >= 10 counts in at least 5 samples
# More stringent than rowSums >= 10; reduces false discoveries
keep <- rowSums(counts(dds) >= 10) >= 5
dds <- dds[keep, ]

# Set reference level explicitly
dds$sample_type <- relevel(factor(dds$sample_type), ref = "Solid Tissue Normal")

dds <- DESeq(dds)

# Check available contrasts before extracting results
resultsNames(dds)

res <- results(dds,
  name = "sample_type_Primary.Tumor_vs_Solid.Tissue.Normal",
  alpha = 0.05
)

# LFC shrinkage for ranking and visualization (not for significance testing)
# apeglm is recommended: Cauchy prior, uses likelihood directly
res_shrunk <- lfcShrink(dds,
  coef = "sample_type_Primary.Tumor_vs_Solid.Tissue.Normal",
  type = "apeglm"
)

# Export significant genes
sig <- subset(as.data.frame(res_shrunk), padj < 0.05 & abs(log2FoldChange) > 1)
sig <- sig[order(sig$padj), ]
write.csv(sig, "results/significant_degs.csv")
message(sprintf("DEGs (padj<0.05, |LFC|>1): %d up, %d down",
  sum(sig$log2FoldChange > 0), sum(sig$log2FoldChange < 0)))
```

### Multi-Factor Design (Batch + Condition)

```r
# When batch information is available, include in design
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat, colData = col_data,
  design = ~ batch + sample_type  # batch first, condition last
)
# DESeq2 accounts for batch in the model; no separate correction needed for DE
```

### Batch Correction for Visualization

```
Batch correction decision:
  For DE testing -> Include batch in DESeq2 design formula. Do NOT pre-correct.
  For PCA/heatmap (log-scale data) -> limma::removeBatchEffect on vst/rlog values
  For downstream tools requiring counts -> ComBat-seq (sva package) on raw counts
  NEVER: ComBat (original) on raw counts. NEVER: removeBatchEffect before DE.
```

```r
# Visualization with batch removal
vsd <- vst(dds, blind = FALSE)
library(limma)
assay(vsd) <- removeBatchEffect(assay(vsd), batch = vsd$batch)
# Use this corrected vsd for PCA and heatmaps only

# ComBat-seq for count-level correction (when other tools need corrected counts)
library(sva)
adjusted_counts <- ComBat_seq(
  counts = as.matrix(counts_mat),
  batch = col_data$batch,
  group = col_data$sample_type  # preserve biological signal
)
```

### Gene ID Conversion

```r
# TCGA uses Ensembl IDs (ENSG...). Convert to symbols for downstream analysis.
library(org.Hs.eg.db)

# Method 1: AnnotationDbi (fast, local)
gene_symbols <- mapIds(org.Hs.eg.db,
  keys = sub("\\..*", "", rownames(res)),  # strip version suffix
  keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")

# Remove unmapped and duplicates
res$symbol <- gene_symbols[sub("\\..*", "", rownames(res))]
res <- res[!is.na(res$symbol) & !duplicated(res$symbol), ]
```

### Pathway and Enrichment Analysis

#### Over-Representation Analysis (ORA)

```r
library(clusterProfiler)
library(org.Hs.eg.db)

# ORA on significant DEGs — use when you have a clear gene list
sig_genes <- res$symbol[res$padj < 0.05 & abs(res$log2FoldChange) > 1]

ego <- enrichGO(
  gene = sig_genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
  ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2
)

# KEGG enrichment (requires Entrez IDs)
entrez_ids <- mapIds(org.Hs.eg.db, keys = sig_genes,
  keytype = "SYMBOL", column = "ENTREZID")
ekegg <- enrichKEGG(gene = entrez_ids[!is.na(entrez_ids)],
  organism = "hsa", pvalueCutoff = 0.05)
```

#### Gene Set Enrichment Analysis (GSEA)

```r
# GSEA on the FULL ranked gene list — never pre-filter
# Rank by: sign(log2FC) * -log10(pvalue) or stat column from DESeq2
gene_list <- res_shrunk$log2FoldChange
names(gene_list) <- res_shrunk$symbol
gene_list <- sort(gene_list[!is.na(names(gene_list))], decreasing = TRUE)

# Remove duplicates (keep highest absolute value)
gene_list <- gene_list[!duplicated(names(gene_list))]

gsea_go <- gseGO(
  geneList = gene_list, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
  ont = "BP", minGSSize = 15, maxGSSize = 500,
  pvalueCutoff = 0.05, eps = 0, seed = TRUE
)
```

#### MSigDB Gene Sets (Hallmark, C2, C5)

```r
library(msigdbr)
# Hallmark gene sets (most interpretable for cancer biology)
h_sets <- msigdbr(species = "Homo sapiens", category = "H")
h_list <- split(h_sets$gene_symbol, h_sets$gs_name)

gsea_h <- GSEA(
  geneList = gene_list,
  TERM2GENE = h_sets[, c("gs_name", "gene_symbol")],
  minGSSize = 15, maxGSSize = 500,
  pvalueCutoff = 0.05, eps = 0, seed = TRUE
)
```

#### Sample-Level Pathway Scoring (GSVA/ssGSEA)

```r
library(GSVA)
# GSVA transforms expression matrix to pathway activity scores per sample
# Input: log2-transformed expression (TPM or vst counts), not raw counts
expr_mat <- assay(vst(dds, blind = FALSE))

gsva_params <- gsvaParam(expr_mat, h_list, kcdf = "Gaussian")
gsva_scores <- gsva(gsva_params)
# Output: pathway x sample matrix of enrichment scores
# Use for: clustering, survival association, heatmaps
```

### Visualization

```r
library(ggplot2); library(ggrepel); library(ComplexHeatmap)

# --- PCA ---
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "sample_type", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = sample_type)) +
  geom_point(size = 2) +
  labs(x = paste0("PC1: ", percentVar[1], "%"),
       y = paste0("PC2: ", percentVar[2], "%")) +
  theme_classic(base_size = 12)

# --- Volcano plot ---
df <- as.data.frame(res_shrunk)
df$sig <- ifelse(df$padj < 0.05 & abs(df$log2FoldChange) > 1,
  ifelse(df$log2FoldChange > 0, "Up", "Down"), "NS")
top10 <- head(df[order(df$padj), ], 10)
ggplot(df, aes(log2FoldChange, -log10(padj), color = sig)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c(Up = "#D55E00", Down = "#0072B2", NS = "grey70")) +
  geom_text_repel(data = top10, aes(label = rownames(top10)),
    size = 3, color = "black", max.overlaps = 15) +
  theme_classic(base_size = 12)
```

### Output Specification

| Output | Format | Description |
|--------|--------|-------------|
| `significant_degs.csv` | CSV | DEGs with log2FC, padj, gene symbol |
| `dds.rds` | RDS | DESeqDataSet object for downstream use |
| `vsd.rds` | RDS | VST-transformed data for visualization |
| `gsva_scores.rds` | RDS | Pathway x sample enrichment matrix |
| `pca_plot.pdf` | PDF | PCA colored by condition |
| `volcano_plot.pdf` | PDF | Volcano with top genes labeled |

### Validation Checks

```
After running expression analysis, verify:
  DEG count: Expect 1000-5000 DEGs for tumor vs normal (varies by cancer type)
    TCGA-LUAD: ~3000-4000 DEGs at padj<0.05, |LFC|>1
    < 100 DEGs: check design formula, sample labels, reference level
    > 10000 DEGs: check for confounders, tumor purity bias
  PCA: tumor and normal should separate on PC1 or PC2
  Known markers: TP53, EGFR, KRAS mutations reflected in expression?
  Direction: oncogenes should trend up in tumor, TSGs down
  Volcano symmetry: roughly symmetric unless strong biological bias expected
```

## Mutation Analysis

### Load and Visualize

```r
library(maftools)  # v2.18+
maf <- read.maf(maf = "TCGA.LUAD.maf.gz", clinicalData = clinical_df)

# Summary statistics
plotmafSummary(maf, rmOutlier = TRUE, addStat = "median")

# Oncoplot: top 20 mutated genes with clinical annotations
oncoplot(maf, top = 20, clinicalFeatures = c("ajcc_pathologic_stage"))

# Lollipop plot for specific gene
lollipopPlot(maf, gene = "TP53", AACol = "HGVSp_Short")

# Rainfall plot (inter-mutation distance)
rainfallPlot(maf, detectChangePoints = TRUE, pointSize = 0.4)
```

### Mutational Signatures

```r
# Trinucleotide context extraction
# Requires >= 30 samples for reliable signature extraction
tnm <- trinucleotideMatrix(maf, ref_genome = "BSgenome.Hsapiens.UCSC.hg38")

# Extract de novo signatures via NMF
# nTry: test range 2-8, select by cophenetic correlation
sig <- extractSignatures(tnm, nTry = 6, plotBestFitRes = FALSE)

# Compare to COSMIC SBS v3.4 signatures
cosim <- compareSignatures(sig, sig_db = "SBS")
# cosim$cosine_similarities: matrix of cosine similarity scores
# Assignment: map each sample to dominant COSMIC signature
```

### Co-occurrence and Mutual Exclusivity

```r
# Somatic interactions: pairwise Fisher's test
somaticInteractions(maf, top = 25, pvalue = c(0.05, 0.01))
```

### Driver Gene Detection

```r
# dN/dS based driver detection
oncodrive_res <- maftools::oncodrive(maf, minMut = 5, pvalMethod = "zscore")
plotOncodrive(oncodrive_res, fdrCutOff = 0.1, useFraction = TRUE)
```

## Copy Number Analysis

```r
# Read GISTIC2.0 output files
gistic <- readGistic(
  gisticAllLesionsFile = "all_lesions.conf_99.txt",
  gisticAmpGenesFile = "amp_genes.conf_99.txt",
  gisticDelGenesFile = "del_genes.conf_99.txt",
  gisticScoresFile = "scores.gistic"
)

gisticChromPlot(gistic)        # genome-wide significance
gisticBubblePlot(gistic)       # bubble plot of recurrent CNVs
gisticOncoPlot(gistic, top = 10)  # oncoplot with focal events

# Integration: add CNV to mutation oncoplot
oncoplot(maf, top = 20, additionalFeature = list(cnv_data))
```

## Methylation Analysis

### Processing (ChAMP Pipeline)

```r
library(ChAMP)

# Load IDAT files (450K or EPIC)
myLoad <- champ.load(directory = "/path/to/idats/", arraytype = "EPIC")

# Quality control
champ.QC(beta = myLoad$beta, pheno = myLoad$pd$Sample_Group)

# Normalization: BMIQ recommended for mixed-type arrays
myNorm <- champ.norm(beta = myLoad$beta, arraytype = "EPIC", method = "BMIQ")

# Singular value decomposition to detect batch effects
champ.SVD(beta = myNorm, pd = myLoad$pd)

# Differentially methylated positions
myDMP <- champ.DMP(beta = myNorm, pheno = myLoad$pd$Sample_Group)

# Differentially methylated regions (more robust than individual probes)
myDMR <- champ.DMR(beta = myNorm, pheno = myLoad$pd$Sample_Group,
  method = "Bumphunter")
```

### Integration: Methylation-Expression Correlation

```r
# Match samples across platforms using TCGA barcode (first 12 characters)
common_patients <- intersect(
  substr(colnames(meth_beta), 1, 12),
  substr(colnames(expr_matrix), 1, 12)
)

# For each gene: correlate promoter methylation (TSS +/- 1500bp) with expression
# Use Spearman correlation (non-parametric, robust to outliers)
# Expected: negative correlation for genes silenced by promoter methylation
# Apply BH correction for multiple testing across all genes
```

## Common Pitfalls

1. **DESeq2 input**: Using FPKM/TPM as input. DESeq2 requires raw integer counts.
2. **Genome build mismatch**: TCGA legacy archive uses hg19, GDC harmonized uses hg38. Never mix builds. Check with `GenomeInfoDb::genome(se)`.
3. **Signature sample size**: Running mutational signature extraction on < 30 samples yields unreliable results.
4. **Methylation cutoffs**: Using fixed delta-beta cutoffs without statistical testing. Use DMP/DMR methods.
5. **Tumor purity**: TCGA samples range 20-95% purity. Use ESTIMATE or ABSOLUTE purity scores to filter (> 60%) or include as covariate.
6. **Pre-filtering GSEA input**: GSEA requires the full ranked gene list. Never pre-filter to significant genes only (that is ORA, not GSEA).
7. **Gene ID mismatch**: TCGA uses Ensembl IDs with version suffixes (ENSG00000141510.18). Strip versions before mapping: `sub("\\..*", "", ensembl_ids)`.
8. **LFC shrinkage for hypothesis testing**: apeglm shrinkage is for ranking and visualization. Use unshrunken results for strict significance calls unless using `lfcThreshold > 0`.

## Related Skills

- `immune-deconvolution`: Estimate immune cell composition from expression data
- `survival-analysis`: Build prognostic models from multi-omics features
- `variant-annotation`: Annotate and classify somatic mutations
- `data-visualization`: Publication-quality volcano plots, heatmaps, KM curves

## Public Datasets for Testing

| Dataset | Samples | Use Case |
|---------|---------|----------|
| TCGA-LUAD | 585 | Lung adenocarcinoma, well-characterized drivers (EGFR, KRAS, ALK) |
| TCGA-BRCA | 1098 | Breast cancer, strong molecular subtypes (Luminal A/B, HER2, Basal) |
| TCGA-GBM | 617 | Glioblastoma, classic for survival and methylation subtyping |
| TCGA-COAD | 521 | Colorectal, MSI-H vs MSS comparison |
