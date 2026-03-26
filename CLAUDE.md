# Biomedical AI Skills

This repository contains validated AI agent skills for cancer biology and multi-omics research.

## Repository Layout

- `skills/` contains individual SKILL.md files organized by domain

## When Working on This Repository

### Adding a New Skill

1. Create directory: `skills/skill-name/SKILL.md`
2. Create test directory: `skills/skill-name/tests/`
3. Use only publicly available datasets for examples
4. Test all code blocks before committing

### Modifying Existing Skills

1. Read the current SKILL.md completely before editing
2. Run existing validation tests after changes
3. Update the tests if behavior changes

### Security Rules

NEVER include in any file:
- Personal file paths or usernames
- Institutional names or server addresses
- API keys, tokens, or credentials
- Unpublished research data references
- Patient identifiers or clinical record numbers

Use generic paths: `/path/to/data/`, `~/analysis/`
Use public accessions: GSE12345, TCGA-LUAD, etc.

### Writing Quality

- No filler phrases ("It is important to note", "In this comprehensive guide")
- No field introductions ("Cancer is a complex disease")
- Direct, actionable content only
- Code must be tested and runnable
- Parameters must be justified, not arbitrary
