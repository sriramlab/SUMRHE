# Contributing

This is a code-only repository. Keep research inputs and outputs outside the
checkout. Do not commit participant identifiers, genotypes, phenotypes,
covariates, participant lists, derived datasets, analysis logs, or archives of
those files. Do not place them in issues, pull requests, releases, or build
artifacts either.

Before introducing any example, review its provenance. Simulated phenotypes,
renamed sample IDs, and subsets of real genotypes do not constitute fully
synthetic data. Prefer a small generator that uses no real participant data;
review the generator before adding it. Clear notebook outputs and attachments
before committing.

The ignore rules reduce accidental additions but do not inspect file contents or
block `git add -f`. Review `git diff --cached` and every new file before pushing.

## Clones made before the history cleanup

The repository history was rewritten to remove bundled datasets. Start from a
fresh clone. Do not merge, mirror-push, or force-push an older clone or restore
removed examples from an old commit. If you have unpublished code changes,
review and transfer only the needed code changes onto the new history without
bringing old ancestry, data, logs, or outputs with them. Keep any old local copy
restricted and disconnected from public push destinations until it is cleaned.
