#!/usr/bin/env bash
# Dry-run release check: build the sdist and wheel from the current tree
# and validate them, without uploading anything.
#
# Requires the dev extra: pip install -e ".[dev]"
#
# Note: the sdist is built from the working tree, so uncommitted or
# untracked files inside the package directory WILL be included. Run
# from a clean checkout to see what a real (CI-built) release contains.
set -euo pipefail

cd "$(dirname "$0")/.."

rm -rf dist/
python -m build
twine check dist/*

echo
echo "Contents of the sdist:"
tar -tzf dist/*.tar.gz | sort

echo
echo "OK: artifacts in dist/ pass twine check. Nothing was uploaded."
echo
echo "TO PUBLISH TO PyPI WITH GITHUB ACTIONS"
echo "--------------------------------------"
echo "# 1. bump the version in pyproject.toml and CITATION.cff, then:
    git add pyproject.toml stellar_geology/__init__.py CITATION.cff
    git commit -m \"Bump version to 0.1.0\"

# 2. create an annotated tag on that commit
    git tag -a v0.1.0 -m \"Release v0.1.0\"

# 3. push commit and tag together
    git push --follow-tags

# 4. GitHub Actions publishes to PyPI every release automatically"
echo