#!/usr/bin/env bash
# Generate or update pixi.lock using the Pixi Docker image (linux-64).
# Run from repo root. Requires Docker. Commit the updated pixi.lock.
set -e
cd "$(dirname "$0")/.."
docker run --rm -v "$(pwd):/app" -w /app ghcr.io/prefix-dev/pixi:0.53.0 pixi lock
echo "pixi.lock updated. Review and commit."
