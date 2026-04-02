# Reproducing the Dramatic Box-Deformation Result

**What this shows:** A pure Martini water box in semi-isotropic NPT deforms its box geometry
over nanoseconds because the neighbor-list artifact creates a fake pressure anisotropy.
Box-Z expands while Box-XY contracts (artifact run), compared to a fixed run using a
generous manual buffer. Both runs start from the same 6×6×6 nm cubic box.

---

## System

- **Force field:** Martini v2.2
- **Molecules:** 1530 Martini water beads (W, P4-type), one bead per molecule
- **Box:** 6 × 6 × 6 nm, periodic boundary conditions
- **Topology file:** `water.top` (1530 W, `martini_v2.2.itp` + `w.itp`)

No membrane, no protein. Pure water was chosen because it isolates the neighbor-list
artifact from all other physics — the only thing that can deform the box is a fake
pressure anisotropy.

---

## Starting structure: `sim/em.gro`

The dramatic runs reuse the energy-minimized structure from `../sim/`. You must build it
once:

### 1. Create a single-bead template

**`single_w.gro`:**
```
Single Martini water
   1
    1W       W    1   0.000   0.000   0.000
   1.0   1.0   1.0
```

### 2. Fill the box

```bash
gmx insert-molecules -ci single_w.gro -box 6 6 6 -nmol 1530 -o water_box.gro
```

Random placement; some beads overlap — that's fine, EM fixes it.

### 3. Energy minimization

```bash
gmx grompp -f ../sim/em.mdp -c water_box.gro -p water.top -o em.tpr
gmx mdrun -v -deffnm em
```

Key `em.mdp` settings:
```ini
integrator  = steep
nsteps      = 5000
emtol       = 1000.0       ; kJ/mol/nm
rcoulomb    = 1.1          ; nm
rvdw        = 1.1
coulombtype = reaction-field
epsilon-rf  = 15
```

Converges in ~48 steps. Energy drops from +24,000,000 → −35,000 kJ/mol.
Output: `em.gro` (used by all dramatic runs as `-c $EMGRO`).

---

## The artifact run: `artifact_npt.mdp`

```ini
integrator      = md
dt              = 0.020          ; 20 fs — standard Martini timestep
nsteps          = 2500000        ; 50 ns (use --short for 5 ns, see below)

; Output — write frequently to see the sawtooth
nstenergy           = 25
nstxout-compressed  = 5000

; Neighbor list — DEFAULT settings, triggers the artifact
cutoff-scheme               = Verlet
nstlist                     = 20
verlet-buffer-tolerance     = 0.005   ; VBT ON → GROMACS auto-tunes

; Nonbonded
rcoulomb        = 1.1
rvdw            = 1.1
coulombtype     = reaction-field
epsilon-rf      = 15
vdwtype         = cutoff
vdw-modifier    = Potential-shift

; Thermostat
tcoupl   = v-rescale
tc-grps  = System
tau-t    = 1.0
ref-t    = 310

; Barostat — SEMI-ISOTROPIC: lets z scale independently from xy
pcoupl        = Parrinello-Rahman
pcoupltype    = semiisotropic
tau-p         = 4.0
ref-p         = 1.0 1.0
compressibility = 4.5e-5 4.5e-5
```

**What GROMACS does with these settings** (confirmed in `artifact_npt.log`):
```
Changing nstlist from 20 to 25, rlist from 1.266 to 1.357
  outer list: updated every 25 steps, buffer 0.257 nm, rlist 1.357 nm
  inner list: updated every  5 steps, buffer 0.031 nm, rlist 1.131 nm
```
VBT auto-tuned `nstlist` → 25 and set up the dual pair list. This is the regime that
produces the artifact.

**Why `semiisotropic` is critical:** with `isotropic` the fake ΔP averages out. With
`semiisotropic` the barostat treats z separately from xy, so it acts on the fake
Pzz ≠ Pxy and physically deforms the box.

---

## The fixed run: `fixed_npt.mdp`

Identical to `artifact_npt.mdp` except:

```ini
; Neighbor list — FIX: disable VBT, use generous manual buffer
nstlist                     = 20
verlet-buffer-tolerance     = -1      ; disables VBT and dual pair list
rlist                       = 1.45    ; generous buffer (0.35 nm over rcoulomb)
```

With `verlet-buffer-tolerance = -1`, GROMACS uses a single pair list with the explicit
`rlist`, no inner/outer split, no auto-tuning. The buffer is wide enough that essentially
no pairs are missed between rebuilds.

---

## Running (with nerfs for a small machine)

```bash
bash run.sh --short   # 5 ns instead of 50 ns
```

The `--short` flag overrides `nsteps` to 250000 (5 ns). This is enough to see
significant deformation. Actual run time on 1 MPI rank + 14 OpenMP threads: ~152 s wall
time (~2840 ns/day performance).

Full command sequence in `run.sh`:
```bash
GMX=~/gromacs-install/bin/gmx_mpi
EMGRO=../sim/em.gro

gmx grompp -f artifact_npt.mdp -c $EMGRO -p water.top -o artifact_npt.tpr -maxwarn 1
gmx mdrun -v -deffnm artifact_npt -nsteps $NSTEPS

gmx grompp -f fixed_npt.mdp -c $EMGRO -p water.top -o fixed_npt.tpr -maxwarn 1
gmx mdrun -v -deffnm fixed_npt -nsteps $NSTEPS

echo -e "Box-X\nBox-Y\nBox-Z\n" | gmx energy -f artifact_npt.edr -o box_artifact.xvg
echo -e "Box-X\nBox-Y\nBox-Z\n" | gmx energy -f fixed_npt.edr -o box_fixed.xvg
```

---

## Results (5 ns run)

| | Start | End Box-X/Y | End Box-Z |
|---|---|---|---|
| **Artifact** | 6.000 × 6.000 × 6.000 nm | 5.483 nm | 6.043 nm |
| **Fixed** | 6.000 × 6.000 × 6.000 nm | 6.130 nm | 4.901 nm |

The artifact run: Box-Z expands (+0.7%), Box-XY contracts (−8.6%) over just 5 ns.
The divergence is caused by the fake Pzz ≠ Pxy driving the semi-isotropic barostat.

Visualize:
```bash
python3 compare_boxes.py   # outputs box_comparison.png
```

---

## File inventory

| File | Role |
|---|---|
| `artifact_npt.mdp` | Artifact run parameters (VBT on, semi-isotropic NPT) |
| `fixed_npt.mdp` | Fixed run parameters (VBT off, generous rlist) |
| `water.top` | Topology: 1530 Martini W beads |
| `run.sh` | Full run script; `--short` flag for 5 ns |
| `compare_boxes.py` | Plots box-X/Y/Z vs time for both runs |
| `box_comparison.png` | Output figure |
| `box_artifact.xvg` / `box_fixed.xvg` | Raw box dimension data |
| `../sim/em.gro` | Minimized starting structure (shared input) |
| `../topology/martini_v2.2.itp` | Martini v2.2 force field |
| `../topology/w.itp` | Martini water bead definition |
