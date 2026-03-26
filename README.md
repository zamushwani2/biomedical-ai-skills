# Biomedical AI Skills

SKILL.md files for cancer bioinformatics. Drop one into your project and your AI coding agent stops getting the normalization order wrong, using the wrong genome build, or running GSEA on filtered gene lists.

Compatible with Claude Code, Cursor, Codex CLI, and Gemini CLI.

## Skills

| Skill | What it covers |
|-------|----------------|
| [`cancer-multiomics`](skills/cancer-multiomics/) | TCGA/GEO data retrieval, differential expression (DESeq2), pathway analysis (GSEA, ORA, GSVA), gene ID conversion, batch correction |

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
