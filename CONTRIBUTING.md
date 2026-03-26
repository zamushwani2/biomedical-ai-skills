# Contributing

Pull requests welcome. Here's what matters.

## Accepted

- New skills for biomedical domains not yet covered
- Validation test cases using publicly available datasets
- Bug fixes (incorrect parameters, broken code, wrong defaults)
- Documentation improvements

## Quality Bar

Every SKILL.md must:

1. Be tested on at least one public dataset with reproducible results
2. Have code examples that run without modification in a clean environment
3. Pin package versions
4. Explain decision points (why method X over Y)

No placeholder content. No "TODO" sections. No filler.

## Submitting

1. Fork and create a branch: `git checkout -b skill/your-skill-name`
2. Add your skill in `skills/your-skill-name/SKILL.md`
3. Include tests in `skills/your-skill-name/tests/`
4. Open a pull request describing what the skill does and how you validated it

## Rules

- No personal file paths, credentials, or institutional details
- No references to unpublished or paywalled data
- All example datasets must be publicly accessible
- See [SECURITY.md](SECURITY.md)
