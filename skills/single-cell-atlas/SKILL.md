# Single-Cell Atlas Construction

QC, preprocessing, and normalization of single-cell RNA-seq data using Seurat v5 (R) and scanpy (Python). Covers MAD-based quality filtering, doublet detection, normalization, and feature selection.

## When to Use This Skill

Activate when the user requests:
- Single-cell RNA-seq quality control and preprocessing
- 10x Chromium data processing (Cell Ranger output)
- Doublet detection and removal
- Normalization of UMI count matrices
- Feature selection for dimensionality reduction
- Processing scRNA-seq data before integration, clustering, or annotation

## Inputs

| Data Type | Format | Source |
|-----------|--------|--------|
| Count matrix | 10x HDF5 (`.h5`), MEX (matrix.mtx + barcodes + features) | Cell Ranger, STARsolo |
| Count matrix | AnnData (`.h5ad`) for scanpy, Seurat object (`.rds`) for R | Processed datasets |
| Reference | PBMC 3k, Tabula Muris, HCA datasets | 10x Genomics, HCA |

---

## Loading Data

### Seurat v5 (R)

```r
library(Seurat)  # v5.4+

# From Cell Ranger output (filtered_feature_bc_matrix/)
obj <- Read10X("path/to/filtered_feature_bc_matrix/") |>
  CreateSeuratObject(project = "sample1", min.cells = 3, min.features = 200)

# From HDF5
obj <- Read10X_h5("path/to/filtered_feature_bc_matrix.h5") |>
  CreateSeuratObject(project = "sample1", min.cells = 3, min.features = 200)

# Seurat v5 uses Assay5 (layers: counts, data, scale.data)
# Access counts: obj[["RNA"]]$counts or LayerData(obj, layer = "counts")
```

### scanpy (Python)

```python
import scanpy as sc  # v1.12+

# From Cell Ranger output
adata = sc.read_10x_mtx("path/to/filtered_feature_bc_matrix/", var_names="gene_symbols")

# From HDF5
adata = sc.read_10x_h5("path/to/filtered_feature_bc_matrix.h5")

# Store raw counts for later
adata.raw = adata.copy()
```

---

## Quality Control

### Computing QC Metrics

#### Seurat

```r
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
# Mouse: pattern = "^mt-"
# Ribosomal: pattern = "^RP[SL]" (optional, rarely used for filtering)

# QC stats are in obj@meta.data:
#   nCount_RNA = total UMI per cell
#   nFeature_RNA = genes detected per cell
#   percent.mt = mitochondrial fraction
```

#### scanpy

```python
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=True, inplace=True)

# Adds to adata.obs:
#   total_counts, n_genes_by_counts, pct_counts_mt
#   log1p_total_counts, log1p_n_genes_by_counts
```

### MAD-Based Filtering

```
Why MAD over fixed thresholds:
  Fixed thresholds (>200 genes, <5% MT) are arbitrary and dataset-specific.
  A threshold that works for PBMCs fails for tumor tissues (higher MT baseline)
  or deep-sequenced datasets (more genes detected per cell).

  MAD (median absolute deviation) adapts to the actual data distribution.
  Flag cells that deviate > N MADs from the median in the problematic direction.
  Default: 3 MADs (~99% of values under normality).
  Permissive: 5 MADs (for heterogeneous tissues).
```

#### Seurat (manual MAD implementation)

```r
# Compute MAD thresholds per metric
mad_outlier <- function(x, nmads = 3, direction = "both") {
  med <- median(x, na.rm = TRUE)
  mad_val <- mad(x, na.rm = TRUE)
  if (direction == "lower") return(x < med - nmads * mad_val)
  if (direction == "higher") return(x > med + nmads * mad_val)
  (x < med - nmads * mad_val) | (x > med + nmads * mad_val)
}

obj$outlier_mt <- mad_outlier(obj$percent.mt, nmads = 3, direction = "higher")
obj$outlier_ncount <- mad_outlier(log1p(obj$nCount_RNA), nmads = 3, direction = "both")
obj$outlier_nfeat <- mad_outlier(log1p(obj$nFeature_RNA), nmads = 3, direction = "lower")

obj$discard <- obj$outlier_mt | obj$outlier_ncount | obj$outlier_nfeat

cat(sprintf("Cells flagged: %d / %d (%.1f%%)\n",
  sum(obj$discard), ncol(obj), 100 * mean(obj$discard)))

obj <- subset(obj, subset = discard == FALSE)
```

#### scanpy (MAD filtering)

```python
import numpy as np

def mad_outlier(adata, metric, nmads=3, direction="higher"):
    """Flag outliers beyond N MADs from median."""
    vals = adata.obs[metric].values
    med = np.median(vals)
    mad = np.median(np.abs(vals - med)) * 1.4826  # scale to match SD
    if direction == "higher":
        return vals > med + nmads * mad
    elif direction == "lower":
        return vals < med - nmads * mad
    else:
        return (vals < med - nmads * mad) | (vals > med + nmads * mad)

adata.obs["outlier_mt"] = mad_outlier(adata, "pct_counts_mt", nmads=3, direction="higher")
adata.obs["outlier_counts"] = mad_outlier(adata, "log1p_total_counts", nmads=3, direction="both")
adata.obs["outlier_genes"] = mad_outlier(adata, "log1p_n_genes_by_counts", nmads=3, direction="lower")

adata.obs["discard"] = adata.obs["outlier_mt"] | adata.obs["outlier_counts"] | adata.obs["outlier_genes"]

print(f"Cells flagged: {adata.obs['discard'].sum()} / {adata.n_obs} "
      f"({100 * adata.obs['discard'].mean():.1f}%)")

adata = adata[~adata.obs["discard"]].copy()
```

#### Bioconductor (scater — dedicated function)

```r
library(scater)  # Bioconductor
sce <- as.SingleCellExperiment(obj)
sce <- addPerCellQCMetrics(sce, subsets = list(mt = grep("^MT-", rownames(sce))))

# isOutlier handles log-transformation and batch-aware thresholds
discard <- isOutlier(sce$subsets_mt_percent, nmads = 3, type = "higher") |
  isOutlier(sce$detected, nmads = 3, type = "lower", log = TRUE) |
  isOutlier(sce$sum, nmads = 3, type = "both", log = TRUE)

sce <- sce[, !discard]
```

### Gene Filtering

```r
# Remove genes detected in very few cells
# min.cells = 3 in CreateSeuratObject handles this for Seurat

# For scanpy:
# sc.pp.filter_genes(adata, min_cells=3)

# Also consider removing:
#   Mitochondrial genes (after using them for QC)
#   Ribosomal genes (optional, sometimes informative)
#   Sex-linked genes (if confounding — XIST, DDX3Y, etc.)
```

---

## Doublet Detection

### scDblFinder (R — top performer in benchmarks)

```r
library(scDblFinder)  # Bioconductor

sce <- as.SingleCellExperiment(obj)
sce <- scDblFinder(sce)

# Results in colData(sce)$scDblFinder.class: "singlet" or "doublet"
# Also: scDblFinder.score (continuous), scDblFinder.cxds_score
table(sce$scDblFinder.class)

# Expected doublet rate for 10x Chromium: ~0.8% per 1000 cells captured
# 5000 cells -> ~4% doublets, 10000 cells -> ~8% doublets

# Transfer back to Seurat
obj$doublet_class <- sce$scDblFinder.class
obj$doublet_score <- sce$scDblFinder.score
obj <- subset(obj, subset = doublet_class == "singlet")
```

### Scrublet (Python — integrated in scanpy)

```python
# Scrublet is now in the core scanpy API (not external)
sc.pp.scrublet(adata, expected_doublet_rate=0.06)

# Adds to adata.obs:
#   predicted_doublet (bool), doublet_score (float)
print(f"Detected doublets: {adata.obs['predicted_doublet'].sum()}")

adata = adata[~adata.obs["predicted_doublet"]].copy()
```

```
Doublet detection order:
  Run BEFORE normalization and feature selection.
  Doublets with mixed transcriptomes distort HVG selection
  and normalization size factors.

  For multi-sample experiments:
    Run doublet detection PER SAMPLE, not on merged data.
    Cross-sample doublets don't exist in droplet-based protocols.
```

---

## Normalization

### Decision Tree

```
Which normalization method?

  Seurat v5 workflow, plan to integrate multiple samples?
    -> SCTransform(obj) — v2 is default in Seurat v5
    -> Replaces NormalizeData + FindVariableFeatures + ScaleData in one call
    -> Returns 3000 variable features by default
    -> Requires glmGamPoi package for speed

  Seurat v5, single sample or exploratory analysis?
    -> NormalizeData(obj) + FindVariableFeatures(obj) + ScaleData(obj)
    -> Log-normalization: counts / total * 10000 -> log1p
    -> Fast, simple, adequate for most single-sample analyses

  scanpy workflow?
    -> sc.pp.normalize_total(adata, target_sum=1e4) + sc.pp.log1p(adata)
    -> Standard log-normalization, analogous to Seurat's LogNormalize

  Need principled size factor estimation (composition bias)?
    -> scran::computeSumFactors() on a SingleCellExperiment
    -> Requires pre-clustering with scran::quickCluster()
    -> Slowest method but handles composition bias best

  Dealing with very large datasets (>500k cells)?
    -> SCTransform with BPCells on-disk storage (Seurat v5)
    -> Or scanpy with backed mode: sc.read_h5ad("data.h5ad", backed="r")
```

### SCTransform v2 (Seurat)

```r
# SCTransform v2 is the default in Seurat v5
# Uses glmGamPoi for fast regularized negative binomial regression
# install.packages("glmGamPoi")  # or BiocManager::install("glmGamPoi")

obj <- SCTransform(obj, verbose = FALSE)
# Replaces NormalizeData + FindVariableFeatures + ScaleData
# Sets default assay to "SCT"
# Variable features stored in VariableFeatures(obj)

# For multiple samples, run SCTransform on each before integration
# Then use IntegrateLayers() or SelectIntegrationFeatures()
```

### Log-Normalization (Seurat)

```r
obj <- NormalizeData(obj)  # LogNormalize, scale.factor = 10000
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj)      # centers and scales, optionally regress out variables
# ScaleData(obj, vars.to.regress = "percent.mt")  # regress out MT% if needed
```

### scanpy Normalization

```python
# Store raw counts first
adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
# adata.X now contains log-normalized data
# adata.layers["counts"] retains raw counts for DE analysis later
```

### scran Deconvolution (Bioconductor)

```r
library(scran)
sce <- as.SingleCellExperiment(obj)

# Pre-cluster to avoid violating the non-DE majority assumption
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters = clusters)
# Size factors in sizeFactors(sce), centered around 1

sce <- logNormCounts(sce)
# logcounts(sce) now contains normalized, log-transformed values
```

---

## Feature Selection

### Seurat (VST)

```r
# After NormalizeData:
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
# After SCTransform: already done (default 3000 features)

# Inspect top variable genes
head(VariableFeatures(obj), 20)

# Visualize variance-mean relationship
VariableFeaturePlot(obj)
```

### scanpy

```python
# Flavor depends on input: seurat_v3 needs raw counts, seurat needs log-normalized
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000)
# Adds adata.var["highly_variable"], adata.var["means"], adata.var["variances"]

# Subset to HVGs for PCA (keep all genes in adata.raw)
adata = adata[:, adata.var["highly_variable"]].copy()
```

### Deviance-Based (scry — more principled alternative)

```r
library(scry)  # Bioconductor
sce <- as.SingleCellExperiment(obj)

# Compute binomial deviance per gene on raw counts
sce <- devianceFeatureSelection(sce, assay = "counts", sorted = TRUE)
# Top genes by deviance are most informative
# Select top 2000-4000 genes
hvgs <- head(rownames(sce), 2000)

# Pairs with GLM-PCA instead of standard PCA
# glmpca <- GLMPCA(counts(sce[hvgs, ]), L = 30, fam = "poi")
```

```
Feature selection guidelines:
  2000-3000 HVGs works for most datasets.
  More HVGs (4000-5000) can help for heterogeneous tissues or rare cell types.
  No single method consistently outperforms others (2025 Nature Methods benchmark).
  Deviance-based selection (scry) is more principled than variance-based
  but less widely adopted.
```

---

## Output Specification

| Output | Format | Description |
|--------|--------|-------------|
| `qc_metrics.csv` | CSV | Per-cell QC: total counts, genes detected, MT%, outlier flags |
| `doublet_scores.csv` | CSV | Doublet probability per cell, singlet/doublet classification |
| `filtered_object.rds` | RDS | Seurat object after QC, doublet removal, normalization |
| `filtered_adata.h5ad` | H5AD | AnnData after QC, doublet removal, normalization |
| `hvg_list.csv` | CSV | Selected highly variable genes |
| `qc_violin.pdf` | PDF | Violin plots of nCount, nFeature, percent.mt |

## Validation Checks

```
After preprocessing PBMC 3k (10x reference dataset), verify:

Cell counts:
  Starting cells: ~2,700
  After QC (MAD 3): ~2,500-2,650 (lose ~2-7%)
  After doublet removal: ~2,400-2,600 (lose ~2-5%)

QC distributions:
  Median genes/cell: 800-1200
  Median UMI/cell: 2000-4000
  Median MT%: 2-5%
  Cells with > 10% MT: < 5% of total (PBMCs have low MT content)

Normalization:
  After log-normalization, expression values should be 0-8 range (log scale)
  SCTransform residuals centered near 0 with SD ~1

Feature selection:
  2000 HVGs should include known markers: IL7R, CD14, MS4A1, CD8A, GNLY, NKG7
  Housekeeping genes (ACTB, GAPDH) should NOT be in top HVGs

Doublet detection:
  Expected doublets for 3k cells: 20-60 (~1-2%)
  If > 10% flagged: check threshold or consider the data has quality issues
```

## Common Pitfalls

### QC
1. **Fixed MT% threshold on all datasets**: A 5% cutoff works for PBMCs but removes real cardiomyocytes (30-40% MT) or hepatocytes (15-20% MT). Use MAD-based filtering or at least check the distribution before setting a threshold.
2. **Not filtering genes**: Removing cells but keeping genes detected in 0-2 cells. These contribute noise to normalization and feature selection. Filter with `min.cells = 3` or `sc.pp.filter_genes(adata, min_cells=3)`.
3. **QC on merged data without batch awareness**: MAD thresholds should ideally be computed per sample/batch. A sample with lower quality shifts the median and hides outliers in better samples. Use `scater::isOutlier(..., batch = sample_id)` or compute MAD per sample.

### Doublets
4. **Doublet detection on merged multi-sample data**: Droplet-based doublets only form within a sample. Running scDblFinder or Scrublet on merged data creates false doublets between cell types that were never in the same droplet.
5. **Trusting doublet scores blindly for rare populations**: Rare cell types with intermediate profiles can be flagged as doublets. Cross-check flagged cells against known markers before removal.

### Normalization
6. **Using TPM/FPKM normalization for scRNA-seq**: Single-cell UMI counts are not length-biased. Standard bulk normalization (TPM, FPKM) is wrong for UMI data. Use library-size normalization or SCTransform.
7. **Running SCTransform without glmGamPoi**: SCTransform v2 requires glmGamPoi for the fast GLM fitting. Without it, the function falls back to a much slower implementation. Install it from Bioconductor.
8. **ScaleData on all genes**: `ScaleData()` only needs to run on variable features for PCA. Scaling all genes wastes memory and time. `ScaleData(obj, features = VariableFeatures(obj))`.

### Feature Selection
9. **Running HVG selection on normalized data with seurat_v3 flavor**: The `seurat_v3` flavor in scanpy expects raw counts, not log-normalized data. Using the wrong flavor silently produces bad HVG rankings.

## Related Skills

- `cancer-multiomics`: Bulk RNA-seq analysis, TCGA data retrieval
- `immune-deconvolution`: Bulk-level immune estimation (compare with single-cell ground truth)

## Public Datasets for Testing

| Dataset | Cells | Use Case |
|---------|-------|----------|
| PBMC 3k (10x) | 2,700 | Standard tutorial dataset, 9 cell types, fast processing |
| PBMC 10k (10x) | 10,000 | Larger PBMC, tests scalability |
| Tabula Muris | 100k+ | Multi-organ mouse atlas, tests batch integration |
| Human Cell Atlas (HCA) | varies | Reference atlases for annotation |
