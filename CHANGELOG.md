# Changelog

All notable changes to this project will be documented in this file.

Format follows [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

### Added
- cancer-multiomics skill: expression analysis with DESeq2 v1.50+, pathway analysis (GSEA/ORA/GSVA/ssGSEA), gene ID conversion, batch correction, visualization
- cancer-multiomics skill: mutation analysis with maftools v2.22+ (MAF handling, TMB, mutational signatures, driver detection), CNV analysis (segment processing, GISTIC2.0, gene-level mapping)
- cancer-multiomics skill: methylation analysis with minfi/ChAMP (450K/EPIC processing, Funnorm, probe filtering), DMPs (limma on M-values), DMRs (DMRcate), CIMP subtyping, methylation-expression integration, probe-bias-corrected pathway analysis (missMethyl)
- cancer-multiomics validation tests: expression (DEG benchmarks), mutation (driver frequencies, TMB), CNV (segment interpretation, gene mapping), methylation (DMP detection, beta-value QC) — all against TCGA-LUAD
- immune-deconvolution skill: unified immunedeconv interface for quanTIseq, EPIC, CIBERSORT, xCell, MCP-counter, TIMER, ESTIMATE; tumor purity correction; cross-method benchmarking; BayesPrism for scRNA-seq-reference-based deconvolution
- immune-deconvolution validation tests: deconvolution methods (quanTIseq, EPIC, MCP-counter output checks, cross-method CD8 correlation), ESTIMATE purity (anticorrelation with immune score, cross-validation with quanTIseq), subtype ordering (Basal vs Luminal A) — all against TCGA-BRCA
- survival-analysis skill: Kaplan-Meier with ggsurvfit, Cox PH with Schoenfeld diagnostics and time-varying coefficients, competing risks (cause-specific + Fine-Gray via tidycmprsk), RMST (survRM2), optimal cutpoint selection (maxstat with validation caveats), forest plots (forestmodel)
- Repository structure, contributing guidelines, security policy
