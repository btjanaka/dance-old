#!/bin/bash

#
# Testing
#

if ! command -v pytest >/dev/null; then
  echo 'pytest not on path; can not format. Please install pytest:'
  echo '    pip install pytest'
  exit 2
fi

pytest

#
# Formatting - YAPF - see
# https://github.com/google/yapf/blob/a4a48d89750cff7ac5239e314b0820d7f3bf408a/plugins/pre-commit.sh
#

# Find all staged Python files, and exit early if there aren't any.
PYTHON_FILES=$(git diff --name-only --cached --diff-filter=AM | grep --color=never '.py$')
if [ -z "${PYTHON_FILES}" ]; then
  exit 0
fi

# Verify that yapf is installed; if not, warn and exit.
if ! command -v yapf >/dev/null; then
  echo 'yapf not on path; can not format. Please install yapf:'
  echo '    pip install yapf'
  exit 2
fi

# Check for unstaged changes to files in the index.
CHANGED_FILES=$(git diff --name-only "${PYTHON_FILES[@]}")
if [ -n "${CHANGED_FILES}" ]; then
  echo 'You have unstaged changes to some files in your commit; skipping '
  echo 'auto-format. Please stage, stash, or revert these changes. You may '
  echo 'find `git stash -k` helpful here.'
  echo 'Files with unstaged changes:' "${CHANGED_FILES[@]}"
  exit 1
fi

# Format all staged files, then exit with an error code if any have uncommitted
# changes.
echo 'Formatting staged Python files . . .'

yapf -i -r "${PYTHON_FILES[@]}"

CHANGED_FILES=$(git diff --name-only "${PYTHON_FILES[@]}")
if [ -n "${CHANGED_FILES}" ]; then
  echo 'Reformatted staged files. Please review and stage the changes.'
  echo 'Files updated: ' "${CHANGED_FILES[@]}"
  exit 1
else
  exit 0
fi
