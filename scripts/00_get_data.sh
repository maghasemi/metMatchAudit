#!/usr/bin/env bash
# Downloads the MetLinkR repository at a fixed commit.
# We use the example input and output files that ship with it.
# Run from the project root: bash scripts/00_get_data.sh

set -euo pipefail

REPO_URL="https://github.com/ncats/MetLinkR"
COMMIT="a6aa145ef8cf76fb762172c6f6dd466c0b30c101"
DEST="data-raw/MetLinkR"

if [ -d "$DEST" ]; then
  echo "$DEST already exists. Delete it first if you want a fresh copy."
  exit 1
fi

git clone "$REPO_URL" "$DEST"
git -C "$DEST" checkout --quiet "$COMMIT"

# Print the commit we ended up on so it shows in the terminal log.
echo "MetLinkR checked out at: $(git -C "$DEST" rev-parse HEAD)"
