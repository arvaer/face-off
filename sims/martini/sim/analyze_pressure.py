#!/usr/bin/env python3
"""Analyze pressure tensor from the neighbor list artifact simulation.

Usage:
    1. Extract pressure:
       echo -e "Pres-XX\nPres-YY\nPres-ZZ\n" | gmx energy -f artifact.edr -o pres.xvg
    2. Run this script:
       python3 analyze_pressure.py
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

data = np.loadtxt('pres.xvg', comments=['#', '@'])
t   = data[:, 0]
pxx = data[:, 1]
pyy = data[:, 2]
pzz = data[:, 3]

# Basic stats
print(f"Mean pressures:  Pxx={pxx.mean():.1f}  Pyy={pyy.mean():.1f}  Pzz={pzz.mean():.1f} bar")
print(f"Std devs:        Pxx={pxx.std():.1f}  Pyy={pyy.std():.1f}  Pzz={pzz.std():.1f} bar")

# PSD
dt   = t[1] - t[0]
freq = np.fft.rfftfreq(len(pxx), d=dt)
psd_xx = np.abs(np.fft.rfft(pxx - pxx.mean()))**2
psd_yy = np.abs(np.fft.rfft(pyy - pyy.mean()))**2
psd_zz = np.abs(np.fft.rfft(pzz - pzz.mean()))**2

# Plot
fig, axes = plt.subplots(2, 1, figsize=(10, 8))

# Time series — first 2 ps (100 steps)
ax = axes[0]
mask = t < 2.0
ax.plot(t[mask], pxx[mask], alpha=0.7, label='Pres-XX')
ax.plot(t[mask], pyy[mask], alpha=0.7, label='Pres-YY')
ax.plot(t[mask], pzz[mask], alpha=0.7, label='Pres-ZZ')
ax.set_xlabel('Time (ps)')
ax.set_ylabel('Pressure (bar)')
ax.set_title('Pressure time series (first 2 ps)')
ax.legend()

# Power spectrum
ax = axes[1]
ax.semilogy(freq, psd_xx, alpha=0.5, label='Pres-XX')
ax.semilogy(freq, psd_yy, alpha=0.5, label='Pres-YY')
ax.semilogy(freq, psd_zz, alpha=0.5, label='Pres-ZZ')
ax.axvline(x=2.0,  color='red',   ls='--', alpha=0.8, label=r'2 ps$^{-1}$ (list rebuild, nstlist=25)')
ax.axvline(x=10.0, color='green', ls='--', alpha=0.8, label=r'~10 ps$^{-1}$ (dynamic pruning)')
ax.set_xlabel(r'Frequency (ps$^{-1}$)')
ax.set_ylabel(r'PSD (bar$^2$)')
ax.set_title('Pressure Power Spectrum')
ax.set_xlim(0, 25)
ax.legend()

plt.tight_layout()
plt.savefig('pressure_analysis.png', dpi=150)
print("Saved pressure_analysis.png")
