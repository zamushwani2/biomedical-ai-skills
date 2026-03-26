# Biomedical AI Skills

SKILL.md files that teach AI coding agents how to do cancer biology and multi-omics research correctly. Compatible with Claude Code, Cursor, Codex CLI, and Gemini CLI.

## Skills

| Skill | What it covers |
|-------|----------------|
| [`cancer-multiomics`](skills/cancer-multiomics/) | TCGA/GEO data retrieval, differential expression (DESeq2), pathway analysis (GSEA, ORA, GSVA), mutation analysis (maftools, TMB, signatures), copy number variation (GISTIC2, segment analysis), gene ID conversion, batch correction |

## Usage

```bash
git clone https://github.com/zamushwani2/biomedical-ai-skills.git

# Add a skill to your project
cp skills/cancer-multiomics/SKILL.md /path/to/your/project/.claude/skills/
```

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md).

## License

MIT
