<div align="center">

# Biomedical AI Skills

SKILL.md files for cancer bioinformatics. Drop one into your project and your AI coding agent handles TCGA data, normalization, and statistics correctly.

[![GitHub Stars](https://img.shields.io/github/stars/zamushwani2/biomedical-ai-skills?style=for-the-badge&logo=github&logoColor=white&labelColor=1a1a2e&color=00d9ff)](https://github.com/zamushwani2/biomedical-ai-skills/stargazers)
[![License](https://img.shields.io/github/license/zamushwani2/biomedical-ai-skills?style=for-the-badge&labelColor=1a1a2e&color=4ecdc4)](LICENSE)
[![Last Commit](https://img.shields.io/github/last-commit/zamushwani2/biomedical-ai-skills?style=for-the-badge&logo=git&logoColor=white&labelColor=1a1a2e&color=ff6b6b)](https://github.com/zamushwani2/biomedical-ai-skills/commits/main)

[![R](https://img.shields.io/badge/R-≥_4.3-276DC3?style=flat-square&logo=r&logoColor=white&labelColor=1a1a2e)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-3.19-87b13f?style=flat-square&labelColor=1a1a2e)](https://bioconductor.org/)
[![TCGA](https://img.shields.io/badge/TCGA-GDC_Portal-e84d3c?style=flat-square&labelColor=1a1a2e)](https://portal.gdc.cancer.gov/)
[![Skills](https://img.shields.io/badge/Skills-1-00d9ff?style=flat-square&labelColor=1a1a2e)](skills/)

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

### [`cancer-multiomics`](skills/cancer-multiomics/)

Multi-omics analysis for [TCGA](https://portal.gdc.cancer.gov/) and [GEO](https://www.ncbi.nlm.nih.gov/geo/) cancer datasets.

| Data type | Tools | Methods |
|-----------|-------|---------|
| **Expression** | [DESeq2](https://bioconductor.org/packages/DESeq2/), [clusterProfiler](https://bioconductor.org/packages/clusterProfiler/) | Differential expression, GSEA, ORA, GSVA, ssGSEA |
| **Mutation** | [maftools](https://bioconductor.org/packages/maftools/) | Driver detection, TMB calculation, mutational signatures |
| **Copy number** | [GISTIC2](https://www.broadinstitute.org/cancer/cga/gistic) | Segment analysis, gain/loss classification, gene-level mapping |
| **Methylation** | [minfi](https://bioconductor.org/packages/minfi/), [ChAMP](https://bioconductor.org/packages/ChAMP/), [DMRcate](https://bioconductor.org/packages/DMRcate/) | DMPs, DMRs, CIMP subtyping, methylation-expression integration |

[Validation tests](skills/cancer-multiomics/tests/) benchmarked against TCGA-LUAD.

## Quick start

```bash
git clone https://github.com/zamushwani2/biomedical-ai-skills.git
```

Copy a skill into your project:

```bash
# Claude Code
cp skills/cancer-multiomics/SKILL.md your-project/.claude/skills/

# Cursor
cp skills/cancer-multiomics/SKILL.md your-project/.cursor/skills/

# Any agent that reads SKILL.md
cp skills/cancer-multiomics/SKILL.md your-project/SKILL.md
```

## What's a SKILL.md?

A file that gives AI coding agents domain knowledge for a specific field. The agent reads it before generating code and follows tested protocols instead of guessing at parameters.

**Without a skill:** agent runs DESeq2 without pre-filtering, skips `lfcShrink()`, uses wrong contrast syntax.
**With a skill:** agent pre-filters low-count genes, applies `apeglm` shrinkage, handles TCGA barcodes correctly.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) and [SECURITY.md](SECURITY.md).

## License

[MIT](LICENSE)
