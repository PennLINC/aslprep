#!/bin/bash
#
# Regenerate .git-blame-ignore-revs with log entries for each ignored rev.
# Run from repository root. The file .git-blame-ignore-revs must exist
# and contain one SHA per line (comment lines starting with # are skipped).
#
# Usage: /bin/bash update_ignore_revs.sh
#

set -e

REVS_FILE=".git-blame-ignore-revs"
OUT_FILE="git-blame-ignore-revs"

if [[ ! -f "$REVS_FILE" ]]; then
    echo "No $REVS_FILE found. Create it with one SHA per line to ignore." >&2
    exit 1
fi

for SHA in $(grep -v "^#" "$REVS_FILE" || true); do
    git log --pretty=format:"# %ad - %ae - %s%n$SHA%n" -n 1 --date=short "$SHA"
done > "$OUT_FILE"

mv "$OUT_FILE" "$REVS_FILE"
echo "Updated $REVS_FILE"
