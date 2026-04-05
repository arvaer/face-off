#!/usr/bin/env bash
# Distribute GROMACS runs across nodes listed in nodes.txt using GNU Parallel.
# Run from sims/paper1/ with the project on the shared NAS mount.
#
# Usage:
#   bash launch_parallel.sh                        # all conditions
#   bash launch_parallel.sh artifact_aniso_pr ...  # specific conditions
set -euo pipefail

cd "$(dirname "$0")"

NODES_FILE="$(pwd)/nodes.txt"
WORK_DIR="$(pwd)"
OMP_THREADS="${OMP_THREADS:-16}"   # 16 physical cores on both nodes

if [ ! -f manifest.json ]; then
  echo "manifest.json not found; run: python3 scripts/prepare_runs.py" >&2
  exit 1
fi

mkdir -p analysis

# Build job list, skipping runs that already have output
python3 - "$@" <<'PY' > /tmp/gromacs_jobs.txt
import json, sys
from pathlib import Path

manifest = json.loads(Path("manifest.json").read_text())
requested = set(sys.argv[1:])

for condition, records in manifest["conditions"].items():
    if requested and condition not in requested:
        continue
    for record in records:
        if not (Path(record["run_dir"]) / "run.edr").exists():
            print(Path(record["run_dir"]) / "run.sh")
PY

TOTAL=$(wc -l < /tmp/gromacs_jobs.txt)
if [ "$TOTAL" -eq 0 ]; then
  echo "No pending runs found." >&2
  exit 0
fi
echo "Queuing $TOTAL runs | nodes: $(grep -v '^#' "$NODES_FILE" | grep -v '^$' | wc -l) | OMP_THREADS=$OMP_THREADS"

export OMP_THREADS

parallel \
  --sshloginfile "$NODES_FILE" \
  --workdir "$WORK_DIR" \
  --env OMP_THREADS \
  --env GMX \
  --joblog analysis/parallel_joblog.txt \
  --resume \
  --eta \
  bash {} \
  :::: /tmp/gromacs_jobs.txt
