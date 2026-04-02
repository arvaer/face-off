#!/usr/bin/env python3
"""
analyze_artifact.py — Run harness at multiple nstlist values and plot
the pressure-tensor discontinuity artifact (Kim et al. 2023, Figure 2).

Usage:
    python analyze_artifact.py --harness ./harness --nsteps 5000 --N 500
"""

import argparse
import subprocess
import numpy as np
import csv
import io
import os
import sys

try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("matplotlib not found — will print stats only", file=sys.stderr)


def run_harness(exe, nstlist, rlist, nsteps, N):
    """Run harness binary and return parsed CSV as dict of arrays."""
    cmd = [exe, str(nstlist), str(rlist), str(nsteps), str(N)]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    reader = csv.DictReader(io.StringIO(result.stdout))
    data = {k: [] for k in reader.fieldnames}
    for row in reader:
        for k, v in row.items():
            data[k].append(float(v))
    return {k: np.array(v) for k, v in data.items()}


def pressure_jump_stats(data, nstlist):
    """
    At each rebuild step, compute |ΔPzz - ΔPxx|.
    Returns mean and std of the jump magnitude.
    """
    steps = data['step'].astype(int)
    Pxx = data['Pxx_bar']
    Pzz = data['Pzz_bar']
    rebuilt = data['rebuilt'].astype(int)

    rebuild_indices = np.where(rebuilt == 1)[0]
    jumps_xx, jumps_zz = [], []
    for idx in rebuild_indices[1:]:   # skip first (no prev data)
        before = idx - 1
        after  = idx
        if before >= 0:
            jumps_xx.append(abs(Pxx[after] - Pxx[before]))
            jumps_zz.append(abs(Pzz[after] - Pzz[before]))

    anisotropy = np.array(jumps_zz) - np.array(jumps_xx)
    return {
        'mean_jump_Pzz': np.mean(jumps_zz),
        'mean_jump_Pxx': np.mean(jumps_xx),
        'mean_anisotropy': np.mean(np.abs(anisotropy)),
        'rms_Pzz': np.std(Pzz),
        'rms_Pxx': np.std(Pxx),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--harness', default='./harness')
    parser.add_argument('--nsteps', type=int, default=5000)
    parser.add_argument('--N',      type=int, default=500)
    parser.add_argument('--rlist',  type=float, default=1.2)
    parser.add_argument('--nstlist-values', nargs='+', type=int,
                        default=[1, 5, 10, 20, 40],
                        dest='nstlist_values')
    args = parser.parse_args()

    if not os.path.isfile(args.harness):
        print(f"ERROR: harness binary not found at {args.harness}")
        print("Compile with:  c++ -O2 -std=c++17 virial_artifact_harness.cpp -o harness")
        sys.exit(1)

    print(f"{'nstlist':>8}  {'<|ΔPzz|> bar':>14}  {'<|ΔPxx|> bar':>14}  "
          f"{'anisotropy':>12}  {'RMS Pzz':>10}")
    print("-" * 65)

    results = {}
    for nstlist in args.nstlist_values:
        data = run_harness(args.harness, nstlist, args.rlist, args.nsteps, args.N)
        stats = pressure_jump_stats(data, nstlist)
        results[nstlist] = (data, stats)
        print(f"{nstlist:>8}  {stats['mean_jump_Pzz']:>14.2f}  "
              f"{stats['mean_jump_Pxx']:>14.2f}  "
              f"{stats['mean_anisotropy']:>12.2f}  "
              f"{stats['rms_Pzz']:>10.2f}")

    if not HAS_MPL:
        return

    # --- Figure 1: pressure vs time for two nstlist values ---
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=False)
    for ax, nstlist in zip(axes, [1, args.nstlist_values[-1]]):
        data, _ = results[nstlist]
        ax.plot(data['step'], data['Pzz_bar'], label='Pzz', lw=0.8, color='steelblue')
        ax.plot(data['step'], data['Pxx_bar'], label='Pxx', lw=0.8, color='tomato', alpha=0.7)
        # Mark rebuild steps
        rebuild_steps = data['step'][data['rebuilt'] == 1]
        ax.vlines(rebuild_steps, *ax.get_ylim(), color='gray', alpha=0.3, lw=0.5)
        ax.set_ylabel('Pressure (bar)')
        ax.set_title(f'nstlist = {nstlist}')
        ax.legend(fontsize=8)
    axes[-1].set_xlabel('MD step')
    fig.suptitle('Pressure tensor vs step — neighbor-list artifact (LJ fluid)', y=1.01)
    plt.tight_layout()
    plt.savefig('pressure_vs_step.png', dpi=150, bbox_inches='tight')
    print("\nSaved: pressure_vs_step.png")

    # --- Figure 2: jump magnitude vs nstlist ---
    nstlist_vals = sorted(results.keys())
    jump_zz = [results[n][1]['mean_jump_Pzz'] for n in nstlist_vals]
    jump_xx = [results[n][1]['mean_jump_Pxx'] for n in nstlist_vals]

    fig2, ax2 = plt.subplots(figsize=(6, 4))
    ax2.plot(nstlist_vals, jump_zz, 'o-', label='|ΔPzz|', color='steelblue')
    ax2.plot(nstlist_vals, jump_xx, 's--', label='|ΔPxx|', color='tomato')
    ax2.set_xlabel('nstlist (steps)')
    ax2.set_ylabel('Mean pressure jump at rebuild (bar)')
    ax2.set_title('Pressure discontinuity vs nstlist')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('jump_vs_nstlist.png', dpi=150)
    print("Saved: jump_vs_nstlist.png")

    # --- Figure 3: PSD to see spike at f = 1/(nstlist*dt) ---
    fig3, ax3 = plt.subplots(figsize=(8, 4))
    dt = 0.002  # 2 fs
    for nstlist in args.nstlist_values:
        data, _ = results[nstlist]
        Pzz = data['Pzz_bar'] - data['Pzz_bar'].mean()
        freqs = np.fft.rfftfreq(len(Pzz), d=dt)  # ps^-1
        psd = np.abs(np.fft.rfft(Pzz))**2
        ax3.semilogy(freqs, psd, lw=0.8, label=f'nstlist={nstlist}',
                     alpha=0.8)
        f_artifact = 1.0 / (nstlist * dt)
        ax3.axvline(f_artifact, color='gray', lw=0.5, ls='--')
    ax3.set_xlim(0, 5)
    ax3.set_xlabel('Frequency (ps⁻¹)')
    ax3.set_ylabel('PSD |Pzz|²')
    ax3.set_title('Power spectral density — artifact spike at f=1/(nstlist·Δt)')
    ax3.legend(fontsize=7)
    ax3.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('psd_artifact.png', dpi=150)
    print("Saved: psd_artifact.png")


if __name__ == '__main__':
    main()
