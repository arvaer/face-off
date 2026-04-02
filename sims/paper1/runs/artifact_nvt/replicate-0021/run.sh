#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

GMX="${GMX:-gmx}"
ROOT="$(cd ../../.. && pwd)"
OMP_THREADS="${OMP_THREADS:-8}"
export OMP_NUM_THREADS="$OMP_THREADS"

"$GMX" grompp   -f run.mdp   -c start.gro   -p "$ROOT/topology/water_antifreeze.top"   -o run.tpr   -maxwarn 1

"$GMX" mdrun -deffnm run -ntmpi 1 -ntomp "$OMP_THREADS"
