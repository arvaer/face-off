#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

if [ ! -f manifest.json ]; then
  echo "manifest.json not found; run: python3 scripts/prepare_runs.py" >&2
  exit 1
fi

if [ "$#" -eq 0 ]; then
  mapfile -t RUNS < <(find runs -type f -name run.sh | sort)
else
  RUNS=()
  for condition in "$@"; do
    while IFS= read -r path; do
      RUNS+=("$path")
    done < <(find "runs/$condition" -type f -name run.sh | sort)
  done
fi

for script in "${RUNS[@]}"; do
  echo "==> $script"
  bash "$script"
done
