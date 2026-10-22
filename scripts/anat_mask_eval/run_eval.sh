#!/usr/bin/env bash
# Score the masks and draw the montages. Usage: run_eval.sh <output dir> [job names...]
set -euo pipefail
OUT=$1
shift
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# The worktree's aslprep (for aslprep.tests.truth_scoring), ahead of the env's editable install
export PYTHONPATH="$(cd "$HERE/../.." && pwd)"
ONLY=()
if [ "$#" -gt 0 ]; then ONLY=(--only "$@"); fi
micromamba run -n aslprep-pvc python "$HERE/evaluate.py" --jobs "$OUT/jobs.json" \
  --masks "$OUT/masks" --out "$OUT/results" "${ONLY[@]}"
micromamba run -n aslprep-pvc python "$HERE/montages.py" --jobs "$OUT/jobs.json" \
  --masks "$OUT/masks" --results "$OUT/results" --out "$OUT/figures"
