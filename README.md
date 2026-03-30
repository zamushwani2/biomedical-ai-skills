<div align="center">

# Biomedical Skills

SKILL.md files for cancer bioinformatics. Drop one into your project and your AI coding agent handles TCGA data, normalization, and statistics correctly.

[![GitHub Stars](https://img.shields.io/github/stars/zamushwani2/biomedical-ai-skills?style=for-the-badge&logo=github&logoColor=white&labelColor=1a1a2e&color=00d9ff)](https://github.com/zamushwani2/biomedical-ai-skills/stargazers)
[![License](https://img.shields.io/github/license/zamushwani2/biomedical-ai-skills?style=for-the-badge&labelColor=1a1a2e&color=4ecdc4)](LICENSE)
[![Last Commit](https://img.shields.io/github/last-commit/zamushwani2/biomedical-ai-skills?style=for-the-badge&logo=git&logoColor=white&labelColor=1a1a2e&color=ff6b6b)](https://github.com/zamushwani2/biomedical-ai-skills/commits/main)

[![R](https://img.shields.io/badge/R-≥_4.3-276DC3?style=flat-square&logo=r&logoColor=white&labelColor=1a1a2e)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.19-87b13f?style=flat-square&labelColor=1a1a2e)](https://bioconductor.org/)
[![TCGA](https://img.shields.io/badge/TCGA-GDC_Portal-e84d3c?style=flat-square&labelColor=1a1a2e)](https://portal.gdc.cancer.gov/)
[![Skills](https://img.shields.io/badge/Skills-3-00d9ff?style=flat-square&labelColor=1a1a2e)](skills/)

**Works with** Claude Code · Cursor · Codex CLI · Gemini CLI

</div>

---

```mermaid
graph LR
    A["Browse skills"] --> B["Copy SKILL.md<br>to your project"] --> C["Agent reads<br>domain protocols"] --> D["Correct code with<br>tested parameters"]
    style A fill:#1a1a2e,stroke:#00d9ff,color:#fff,stroke-width:2px
    style B fill:#1a1a2e,stroke:#4ecdc4,color:#fff,stroke-width:2px
    style C fill:#1a1a2e,stroke:#ff6b6b,color:#fff,stroke-width:2px
    style D fill:#1a1a2e,stroke:#87b13f,color:#fff,stroke-width:2px
```

## Skills

| Skill | Description | Tests |
|-------|-------------|-------|
| [`cancer-multiomics`](skills/cancer-multiomics/) | Multi-omics analysis for [TCGA](https://portal.gdc.cancer.gov/)/[GEO](https://www.ncbi.nlm.nih.gov/geo/) — expression ([DESeq2](https://bioconductor.org/packages/DESeq2/)), mutation ([maftools](https://bioconductor.org/packages/maftools/)), CNV ([GISTIC2](https://www.broadinstitute.org/cancer/cga/gistic)), methylation ([minfi](https://bioconductor.org/packages/minfi/), [DMRcate](https://bioconductor.org/packages/DMRcate/)) | [TCGA-LUAD](skills/cancer-multiomics/tests/) |
| [`immune-deconvolution`](skills/immune-deconvolution/) | Tumor microenvironment estimation via [immunedeconv](https://omnideconv.org/immunedeconv/) — [quanTIseq](https://icbi.i-med.ac.at/software/quantiseq/doc/), [EPIC](https://github.com/GfellerLab/EPIC), CIBERSORT, [xCell](https://xcell.ucsf.edu/), [MCP-counter](https://github.com/ebecht/MCPcounter), TIMER, ESTIMATE, tumor purity correction | [TCGA-BRCA](skills/immune-deconvolution/tests/) |
| [`survival-analysis`](skills/survival-analysis/) | Time-to-event analysis — Kaplan-Meier ([ggsurvfit](https://www.danieldsjoberg.com/ggsurvfit/)), Cox PH ([survival](https://cran.r-project.org/package=survival)), competing risks ([tidycmprsk](https://mskcc-epi-bio.github.io/tidycmprsk/)), RMST ([survRM2](https://cran.r-project.org/package=survRM2)), optimal cutpoints, forest plots | [TCGA-GBM](skills/survival-analysis/tests/) |

## Quick start

```bash
git clone https://github.com/zamushwani2/biomedical-ai-skills.git
```

Copy one or more skills into your project:

```bash
# Claude Code — single skill
cp skills/cancer-multiomics/SKILL.md your-project/.claude/skills/cancer-multiomics.md

# Claude Code — full pipeline (multiomics → deconvolution → survival)
cp skills/cancer-multiomics/SKILL.md your-project/.claude/skills/
cp skills/immune-deconvolution/SKILL.md your-project/.claude/skills/
cp skills/survival-analysis/SKILL.md your-project/.claude/skills/

# Cursor
cp skills/cancer-multiomics/SKILL.md your-project/.cursor/skills/

# Any agent that reads SKILL.md
cp skills/cancer-multiomics/SKILL.md your-project/SKILL.md
```

Skills cross-reference each other — the agent can chain multiomics data retrieval through deconvolution into survival modeling.

## What's a SKILL.md?

A file that gives AI coding agents domain knowledge for a specific field. The agent reads it before generating code and follows tested protocols instead of guessing at parameters.

**Without a skill:** agent runs DESeq2 without pre-filtering, skips `lfcShrink()`, uses wrong contrast syntax.
**With a skill:** agent pre-filters low-count genes, applies `apeglm` shrinkage, handles TCGA barcodes correctly.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) and [SECURITY.md](SECURITY.md).

## License

[MIT](LICENSE)
