#!/bin/bash
# Run both artifact and fixed NPT simulations side by side.
# Usage:
#   bash run.sh           # full run (50 ns)
#   bash run.sh --short   # quick test (5 ns)

GMX=gmx
EMGRO=../sim/em.gro   # reuse minimized water box from sim/

# Parse flags
NSTEPS=2500000   # 50 ns default
LABEL="50 ns"
for arg in "$@"; do
    if [ "$arg" = "--short" ]; then
        NSTEPS=250000   # 5 ns
        LABEL="5 ns (short)"
    fi
done

echo "=== Run length: $LABEL ($NSTEPS steps) ==="

# --- Artifact run (default VBT settings) ---
echo "=== Preparing artifact run ==="
$GMX grompp -f artifact_npt.mdp -c $EMGRO -p water.top -o artifact_npt.tpr \
    -maxwarn 1
echo "=== Running artifact ($LABEL, default settings) ==="
$GMX mdrun -v -deffnm artifact_npt -nsteps $NSTEPS

# --- Fixed run (generous rlist, no VBT) ---
echo "=== Preparing fixed run ==="
$GMX grompp -f fixed_npt.mdp -c $EMGRO -p water.top -o fixed_npt.tpr \
    -maxwarn 1
echo "=== Running fixed ($LABEL, generous rlist) ==="
$GMX mdrun -v -deffnm fixed_npt -nsteps $NSTEPS

# --- Extract box dimensions ---
echo "=== Extracting box dimensions ==="
echo -e "Box-X\nBox-Y\nBox-Z\n" | $GMX energy -f artifact_npt.edr -o box_artifact.xvg
echo -e "Box-X\nBox-Y\nBox-Z\n" | $GMX energy -f fixed_npt.edr -o box_fixed.xvg

echo "=== Done. Run: python3 compare_boxes.py ==="
