#!/usr/bin/env bash
# List the references and compute ASLPrep's SynthStrip and legacy masks for each, in the image.
# Usage: run_masks.sh <probe dir> <output dir> <image>
set -euo pipefail
PROBE=$1
OUT=$2
IMAGE=$3
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
mkdir -p "$OUT/masks" "$OUT/work"
micromamba run -n aslprep-pvc python "$HERE/jobs.py" --probe "$PROBE" --out "$OUT/jobs.json"
docker run --rm -u "$(id -u):$(id -g)" -e HOME=/tmp \
  -v "$PROBE:$PROBE:ro" -v "$OUT:$OUT" -v "$HERE:/scripts:ro" \
  --entrypoint python "$IMAGE" /scripts/aslmasks.py \
  --jobs "$OUT/jobs.json" --out "$OUT/masks" --work "$OUT/work" --nprocs 4
