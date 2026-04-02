# CS2050 HPC Final: CUDA Neighbor List Benchmark

**Course:** CS2050 — High Performance Computing (Harvard)
**Approach:** Write a CUDA cell-list neighbor search kernel from scratch, then benchmark it against GROMACS Verlet and MACE-OFF on the same water system.

## The Question

Which neighbor-search strategy wins on modern GPU hardware — classical cell-list, GROMACS's optimized Verlet with MxN SIMD clustering, or MACE-OFF's per-step graph construction?

## What's Being Built

| Component | File | What |
|-----------|------|------|
| Cell-list kernel | `src/cell_list.cu` | `compute_cell_index`, `find_neighbors`, `compute_lj_forces` |
| Validation | `src/validate.py` | Compare vs GROMACS reference |
| Benchmark harness | `src/benchmark.py` | 3-way × 4 scales × 3 runs |

## Benchmarks

**3 approaches:**
- Your CUDA cell-list kernel (this repo)
- GROMACS Verlet (VBT=0.005, nstlist=25, rlist=1.269nm)
- MACE-OFF medium (Kovács, Moore, **Witt** et al., arXiv:2312.15211)

**4 system sizes:** N = 1530, 6120, 24480, 97920 (SPC/E water, same density + rcut)

**Metrics:** neighbor search time, force eval time, total step time, GPU utilization, pressure tensor accuracy

## Compute Resources

- Local: RTX 5090, A400
- Cluster: 160× T4 (20 nodes, 8 T4 each)
- Profiling: Nsight Compute

## Phases

| Phase | What | Target |
|-------|------|--------|
| 2 | CUDA kernel + validation | Apr 10 |
| 3 | Profiling + optimization | Apr 19 |
| 4 | Scaling study + figures | Apr 26 |
| 5 | Report + presentation | May 3 |

## Background Reading

- Kim, Fábián, Hummer (2023) — Neighbor list artifacts in MD
- Kovács, Moore, Witt et al. (2023) — MACE-OFF (arXiv:2312.15211)
