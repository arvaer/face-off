#!/usr/bin/env python3
"""Compare box dimensions between artifact and fixed NPT runs.

If the artifact is present, the semi-isotropic barostat will cause
Box-Z to drift away from Box-X/Y over time. The fixed run should
keep all three dimensions equal (cubic box stays cubic).
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

art = np.loadtxt('box_artifact.xvg', comments=['#', '@'])
fix = np.loadtxt('box_fixed.xvg', comments=['#', '@'])

fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)

# --- Artifact ---
ax = axes[0]
t = art[:, 0] / 1000  # ps → ns
ax.plot(t, art[:, 1], label='Box-X', alpha=0.8)
ax.plot(t, art[:, 2], label='Box-Y', alpha=0.8)
ax.plot(t, art[:, 3], label='Box-Z', alpha=0.8, color='red')
ax.set_xlabel('Time (ns)')
ax.set_ylabel('Box dimension (nm)')
ax.set_title('Default settings (artifact)')
ax.legend()
ax.axhline(y=art[0, 1], color='gray', ls=':', alpha=0.5)

# --- Fixed ---
ax = axes[1]
t = fix[:, 0] / 1000
ax.plot(t, fix[:, 1], label='Box-X', alpha=0.8)
ax.plot(t, fix[:, 2], label='Box-Y', alpha=0.8)
ax.plot(t, fix[:, 3], label='Box-Z', alpha=0.8, color='red')
ax.set_xlabel('Time (ns)')
ax.set_title('Fixed settings (generous rlist, no VBT)')
ax.legend()
ax.axhline(y=fix[0, 1], color='gray', ls=':', alpha=0.5)

plt.suptitle('Box deformation from neighbor list artifact\n'
             '(semi-isotropic Parrinello-Rahman barostat, 1530 Martini waters)',
             fontsize=13)
plt.tight_layout()
plt.savefig('box_comparison.png', dpi=150)
print("Saved box_comparison.png")

# Print summary
for name, d in [('Artifact', art), ('Fixed', fix)]:
    bx, by, bz = d[-1, 1], d[-1, 2], d[-1, 3]
    bx0, by0, bz0 = d[0, 1], d[0, 2], d[0, 3]
    print(f"\n{name}:")
    print(f"  Start: {bx0:.3f} x {by0:.3f} x {bz0:.3f} nm")
    print(f"  End:   {bx:.3f} x {by:.3f} x {bz:.3f} nm")
    print(f"  Drift: ΔX={bx-bx0:+.3f}  ΔY={by-by0:+.3f}  ΔZ={bz-bz0:+.3f} nm")
