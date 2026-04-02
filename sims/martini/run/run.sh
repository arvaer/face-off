#!/usr/bin/env bash
# Run MARTINI DPPC sim: EM -> NVT -> NPT -> short production -> art export
# Usage: bash run.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"
DEFAULT_BLENDER_BIN="/Applications/Blender.app/Contents/MacOS/Blender"
GMX="${GMX:-/Users/mikeyalmeida/gromacs-install/bin/gmx_mpi}"
TOP="$SCRIPT_DIR/../topology"
MDP="$SCRIPT_DIR/../mdp"
ART_DIR="$SCRIPT_DIR/../art_prod"
BLEND_FILE="$ART_DIR/md_signal_art.blend"
if [[ -x "$DEFAULT_BLENDER_BIN" ]]; then
  BLENDER_BIN="${BLENDER_BIN:-$DEFAULT_BLENDER_BIN}"
else
  BLENDER_BIN="${BLENDER_BIN:-}"
fi
PROD_EXTEND_PS="${PROD_EXTEND_PS:-1000}"
ART_SOLVENT_STRIDE="${ART_SOLVENT_STRIDE:-6}"
ART_SOLVENT_SHELL_NM="${ART_SOLVENT_SHELL_NM:-0.75}"
ART_MAX_SOLVENT="${ART_MAX_SOLVENT:-24}"
UV_CACHE_DIR="${UV_CACHE_DIR:-/tmp/uv-cache}"
export OMPI_MCA_btl="${OMPI_MCA_btl:-self}"
export OMPI_MCA_oob="${OMPI_MCA_oob:-^tcp}"

cd "$SCRIPT_DIR"
mkdir -p "$UV_CACHE_DIR"

echo "=== Step 1: Energy Minimization ==="
$GMX grompp -f $MDP/em.mdp -c $TOP/system.gro -p $TOP/system.top -o em.tpr -maxwarn 2
$GMX mdrun -v -deffnm em

echo "=== Step 2: NVT Equilibration ==="
$GMX grompp -f $MDP/nvt.mdp -c em.gro -p $TOP/system.top -o nvt.tpr -maxwarn 2
$GMX mdrun -v -deffnm nvt

echo "=== Step 3: NPT Equilibration ==="
$GMX grompp -f $MDP/npt.mdp -c nvt.gro -t nvt.cpt -p $TOP/system.top -o npt.tpr -maxwarn 2
$GMX mdrun -v -deffnm npt

echo "=== Step 4: Short Production MD ==="
$GMX convert-tpr -s npt.tpr -extend "$PROD_EXTEND_PS" -o prod.tpr
$GMX mdrun -v -deffnm prod -s prod.tpr -cpi npt.cpt -noappend
PROD_TRAJECTORY="$(find "$SCRIPT_DIR" -maxdepth 1 -type f -name 'prod*.xtc' | sort | tail -n 1)"

if [[ -z "$PROD_TRAJECTORY" ]]; then
  echo "No production trajectory was written."
  exit 1
fi

echo "=== Step 5: Export production trajectory to art scene ==="
UV_CACHE_DIR="$UV_CACHE_DIR" uv run --python python3 "$ROOT_DIR/md_signal_art/export_gromacs_signal.py" \
  --trajectory "$PROD_TRAJECTORY" \
  --structure "$SCRIPT_DIR/prod.tpr" \
  --gmx-bin "$GMX" \
  --solvent-stride "$ART_SOLVENT_STRIDE" \
  --solvent-shell-nm "$ART_SOLVENT_SHELL_NM" \
  --max-solvent "$ART_MAX_SOLVENT" \
  --output "$ART_DIR"

UV_CACHE_DIR="$UV_CACHE_DIR" uv run --python python3 "$ROOT_DIR/md_signal_art/build_retro_viewer.py" \
  --scene "$ART_DIR/scene.json" \
  --output "$ART_DIR/retro_viewer.html"

if [[ -n "$BLENDER_BIN" && -x "$BLENDER_BIN" ]]; then
  echo "=== Step 6: Build Blender scene ==="
  "$BLENDER_BIN" --background --python "$ROOT_DIR/md_signal_art/blender_import_points.py" -- \
    "$ART_DIR/scene.json" "$BLEND_FILE"
else
  echo "=== Step 6: Blender handoff ==="
  echo "Set BLENDER_BIN to build the .blend automatically."
  echo "Command: blender --background --python $ROOT_DIR/md_signal_art/blender_import_points.py -- $ART_DIR/scene.json $BLEND_FILE"
fi

echo "=== Done. Production art assets are in $ART_DIR ==="
