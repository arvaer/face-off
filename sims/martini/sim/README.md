# Reproducing the Neighbor List Pressure Artifact

Reproduction of the key result from Kim, Fabian & Hummer, "Neighbor List Artifacts in Molecular Dynamics Simulations," *J. Chem. Theory Comput.* 2023, 19, 8919-8929.

This guide builds a pure Martini water box from scratch and runs it with GROMACS's default neighbor list settings to observe anisotropic pressure oscillations.

---

## Background: What is the artifact?

GROMACS does not recompute all pairwise interactions every timestep. Instead, it builds a **neighbor list** — a catalogue of atom pairs close enough to interact — and reuses it for several steps. This is fast, but between rebuilds some pairs drift close enough to interact without being on the list. Their forces are silently missed.

GROMACS controls the error budget with the **Verlet Buffer Tolerance (VBT)**, which auto-tunes the list radius to keep total *energy* drift small. The problem: energy is a scalar, but pressure is a tensor. The missed interactions are not distributed uniformly across directions — GROMACS's SIMD clustering sorts atoms along *z*, making clusters elongated in that direction. The effective buffer is thinner along *z*, so more pairs are missed along *z* than along *x* or *y*. The result: $P_{zz}$ has a different error than $P_{xx}$ and $P_{yy}$.

On top of this directional bias, the pressure oscillates at specific frequencies:
- Every `nstlist` steps (default 25), the list is rebuilt and all pairs are suddenly accounted for, creating a discontinuity.
- GROMACS also maintains a shorter "inner" pair list that is pruned every ~5 steps, creating a second oscillation.

In an NPT simulation with a semi-isotropic barostat, the barostat sees the fake $P_{zz} \neq P_{xx}$ and "corrects" it by deforming the box. Over hundreds of nanoseconds, this crumples membranes (paper Figure 1A).

This simulation reproduces the artifact in its simplest form: pure Martini water in NVT, where you can observe the raw pressure oscillations without a barostat masking them.

---

## Prerequisites

- GROMACS (any recent version; we used 2027.0-dev)
- Martini v2.2 force field files in `../topology/`:
  - `martini_v2.2.itp` (the force field)
  - `w.itp` (Martini water bead definition: a single P4-type bead)

---

## Step 1: Create a single water bead template

GROMACS's `insert-molecules` tool stamps copies of a template molecule into a box. We need a `.gro` file containing exactly one Martini water bead.

**File: `single_w.gro`**
```
Single Martini water
   1
    1W       W    1   0.000   0.000   0.000
   1.0   1.0   1.0
```

This is one bead of type W at the origin. The box dimensions (1.0 nm) don't matter — `insert-molecules` will override them.

---

## Step 2: Fill a box with 1530 water beads

```bash
gmx insert-molecules -ci single_w.gro -box 6 6 6 -nmol 1530 -o water_box.gro
```

This creates a 6x6x6 nm periodic box and randomly places 1530 Martini water beads. The paper uses this system size. The beads are placed with random positions but no velocities — and some may overlap slightly.

**Output**: `water_box.gro` (1530 atoms in a 6 nm cubic box).

---

## Step 3: Write the topology

The topology tells GROMACS what force field to use and how many of each molecule type are in the system. It must match the coordinate file exactly.

**File: `water.top`**
```
; Pure Martini water topology

#include "../topology/martini_v2.2.itp"
#include "../topology/w.itp"

[ system ]
Martini water

[ molecules ]
; name  nmol
W       1530
```

The `#include` lines pull in the Martini force field parameters and the W molecule definition (one P4 bead per molecule). The `[ molecules ]` section declares 1530 copies.

---

## Step 4: Energy minimization

Random insertion creates overlapping beads with enormous repulsive forces (initial potential energy ~ +24,000,000 kJ/mol). If you try to run MD directly, the forces are so large that the simulation crashes immediately (segfault from numerical overflow).

Energy minimization (EM) iteratively moves atoms downhill on the potential energy surface until forces are reasonable.

**File: `em.mdp`**
```ini
integrator  = steep        ; steepest descent — simple, robust
nsteps      = 5000         ; max iterations (usually converges much faster)
emtol       = 1000.0       ; stop when max force < 1000 kJ/mol/nm
emstep      = 0.01         ; initial step size (nm)

; Same nonbonded settings as production
cutoff-scheme = Verlet
nstlist     = 10
rcoulomb    = 1.1
rvdw        = 1.1
coulombtype = reaction-field
epsilon-rf  = 15
vdwtype     = cutoff
vdw-modifier = Potential-shift
```

```bash
gmx grompp -f em.mdp -c water_box.gro -p water.top -o em.tpr
gmx mdrun -v -deffnm em
```

EM converges in ~48 steps. The potential energy drops from +24,000,000 to -35,000 kJ/mol — the overlaps are resolved and the system sits in a physically reasonable configuration.

**Output**: `em.gro` (minimized coordinates, ready for MD).

---

## Step 5: Run the artifact simulation

**File: `artifact.mdp`**
```ini
; Integrator
integrator      = md
dt              = 0.020          ; 20 fs — standard Martini timestep
nsteps          = 500000         ; 10 ns

; Output — CRITICAL: write every step
nstenergy       = 1              ; default (100) aliases the oscillation
nstxout-compressed = 5000

; Neighbor list — these defaults cause the artifact
cutoff-scheme   = Verlet
nstlist         = 20             ; VBT auto-tunes this to 25
verlet-buffer-tolerance = 0.005  ; activates VBT → rlist ≈ 1.27 nm

; Nonbonded
rcoulomb        = 1.1
rvdw            = 1.1
coulombtype     = reaction-field
epsilon-rf      = 15
vdwtype         = cutoff
vdw-modifier    = Potential-shift

; Thermostat
tcoupl          = v-rescale
tc-grps         = System
tau-t           = 1.0
ref-t           = 310

; No barostat — pure NVT to see raw pressure errors
pcoupl          = no
```

Key settings and why they matter:

| Parameter | Value | Role in the artifact |
|---|---|---|
| `verlet-buffer-tolerance` | 0.005 | Activates VBT. GROMACS auto-tunes `nstlist` to 25 and `rlist` to ~1.27 nm. This buffer is wide enough for energy conservation but too narrow for pressure accuracy. |
| `nstlist` | 20 (tuned to 25) | The list is rebuilt every 25 steps = 0.5 ps. Between rebuilds, pairs can drift into range unnoticed. |
| `nstenergy` | 1 | Writes pressure every step. The default (100) would alias the 25-step oscillation, hiding the artifact entirely. |
| `pcoupl` | no | NVT: no barostat to mask or amplify the pressure error. We observe the raw signal. |

```bash
gmx grompp -f artifact.mdp -c em.gro -p water.top -o artifact.tpr
gmx mdrun -v -deffnm artifact
```

GROMACS will report:
```
Changing nstlist from 20 to 25, rlist from 1.268 to 1.359
```
This confirms VBT is active and has auto-tuned the neighbor list parameters to the regime identified by the paper.

---

## Step 6: Extract and analyze pressure

After the run completes:

```bash
echo -e "Pres-XX\nPres-YY\nPres-ZZ\n" | gmx energy -f artifact.edr -o pres.xvg
```

### Power spectral density

```python
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('pres.xvg', comments=['#', '@'])
t    = data[:, 0]   # time in ps
pxx  = data[:, 1]   # Pres-XX in bar
pyy  = data[:, 2]   # Pres-YY
pzz  = data[:, 3]   # Pres-ZZ

dt   = t[1] - t[0]  # 0.02 ps
freq = np.fft.rfftfreq(len(pxx), d=dt)  # in ps^-1

psd_xx = np.abs(np.fft.rfft(pxx - pxx.mean()))**2
psd_zz = np.abs(np.fft.rfft(pzz - pzz.mean()))**2

plt.figure(figsize=(10, 5))
plt.semilogy(freq, psd_xx, alpha=0.7, label='Pres-XX')
plt.semilogy(freq, psd_zz, alpha=0.7, label='Pres-ZZ')
plt.xlabel('Frequency (ps$^{-1}$)')
plt.ylabel('PSD (bar$^2$)')
plt.title('Pressure Power Spectrum — Neighbor List Artifact')
plt.axvline(x=2.0,  color='red',   ls='--', label='f = 1/(25 × 0.02 ps) = 2 ps⁻¹ (list rebuild)')
plt.axvline(x=10.0, color='green', ls='--', label='f ≈ 10 ps⁻¹ (dynamic pruning)')
plt.legend()
plt.xlim(0, 25)
plt.tight_layout()
plt.savefig('pressure_psd.png', dpi=150)
plt.show()
```

### What you should see

1. **Spikes at 2 ps^-1 and its harmonics** (4, 6, 8, ...): these come from the outer neighbor list rebuilding every 25 steps = 0.5 ps.

2. **Spikes near 10 ps^-1**: these come from dynamic pruning of the inner pair list every ~5 steps = 0.1 ps.

3. **Pres-ZZ has different amplitude than Pres-XX/YY**: the z-axis clustering bias makes more interactions missed along z.

This matches the paper's Figure 5.

---

## What would fix the artifact?

For reference (not run here), three approaches:

**Disable VBT and widen the buffer:**
```ini
verlet-buffer-tolerance = -1     ; disables VBT and dual pair list
rlist                   = 1.422  ; generous manual buffer
```

**Rebuild every step (expensive):**
```ini
nstlist                 = 1
verlet-buffer-tolerance = -1
rlist                   = 1.1    ; no buffer needed
```

**Synchronize barostat with list rebuilds (practical for NPT):**
```ini
nsttcouple = 25   ; = nstlist
nstpcouple = 25   ; = nstlist
```

---

## File inventory

| File | What it is |
|---|---|
| `single_w.gro` | One Martini water bead — template for `insert-molecules` |
| `water_box.gro` | 1530 beads in 6x6x6 nm box (randomly placed, not minimized) |
| `water.top` | Topology: Martini v2.2 FF + 1530 W molecules |
| `em.mdp` | Energy minimization parameters |
| `em.gro` | Minimized coordinates (input for MD) |
| `artifact.mdp` | Production MD parameters (artifact-triggering defaults) |
| `artifact.tpr` | Run input file |
| `artifact.edr` | Energy/pressure output (one frame per step) |
| `artifact.xtc` | Trajectory (compressed, every 5000 steps) |
| `artifact.log` | Run log — check for `nstlist` and `rlist` auto-tuning |
