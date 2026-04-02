/*
 * virial_artifact_harness.cpp
 *
 * Standalone harness for measuring the neighbor-list virial artifact
 * described in Kim, Fábián, Hummer JCTC 2023.
 *
 * What it does:
 *   1. Builds a random box of Lennard-Jones particles.
 *   2. Constructs a Verlet pair list with configurable nstlist / rlist.
 *   3. Evaluates nonbonded forces + virial using the plain-C 1x1 kernel.
 *   4. Rebuilds the list every nstlist steps and records:
 *        - ΔP_xx, ΔP_zz at each rebuild (the "jump")
 *        - RMS pressure over the run
 *   5. Prints a CSV to stdout: step, Pxx, Pyy, Pzz, rebuilt
 *
 * Compile (outside GROMACS build tree — pure stdlib/Eigen):
 *   c++ -O2 -std=c++17 virial_artifact_harness.cpp -o harness
 *
 * Run:
 *   ./harness <nstlist> <rlist_nm> <nsteps> <N_particles>
 *   ./harness 20 1.2 5000 500
 *
 * Physics notes:
 *   - Pressure = (2K + W) / (3V), but here we track W_αβ only (virial tensor).
 *   - W_αβ = Σ_{i<j} r_{ij,α} * F_{ij,β}
 *   - Lennard-Jones:  U = 4ε[(σ/r)^12 - (σ/r)^6],  F = -dU/dr * r̂
 *   - We use a hard cutoff at rlist (no shift/switch) to maximise the
 *     discontinuity signal — just like GROMACS Verlet without vdw-modifier.
 *   - The artifact is the JUMP in W_αβ at rebuild steps due to pairs that
 *     were just inside rlist last step suddenly contributing at the new step.
 */

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

struct Vec3 { double x, y, z; };

struct Particle {
    Vec3 pos;
    Vec3 vel;
    Vec3 force;
};

// Pair list entry
struct Pair { int i, j; };

// ---------------------------------------------------------------------------
// Minimal periodic box — cubic, side L
// ---------------------------------------------------------------------------

static inline Vec3 pbc(Vec3 dr, double L) {
    auto wrap = [&](double d) {
        d -= L * std::round(d / L);
        return d;
    };
    return {wrap(dr.x), wrap(dr.y), wrap(dr.z)};
}

// ---------------------------------------------------------------------------
// Build pair list (O(N^2) — fine for small N used in harness)
// ---------------------------------------------------------------------------

static std::vector<Pair> build_pairlist(
        const std::vector<Particle>& p, double rlist, double L)
{
    std::vector<Pair> pairs;
    const double rlist2 = rlist * rlist;
    const int N = static_cast<int>(p.size());
    for (int i = 0; i < N - 1; ++i) {
        for (int j = i + 1; j < N; ++j) {
            Vec3 dr = pbc({p[j].pos.x - p[i].pos.x,
                           p[j].pos.y - p[i].pos.y,
                           p[j].pos.z - p[i].pos.z}, L);
            double r2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
            if (r2 < rlist2) pairs.push_back({i, j});
        }
    }
    return pairs;
}

// ---------------------------------------------------------------------------
// Evaluate forces + virial from current pair list
// LJ: ε=1, σ=0.34 nm (argon-like), hard cutoff at rvdw = rlist - buffer
// ---------------------------------------------------------------------------

static void eval_forces(
        std::vector<Particle>& p,
        const std::vector<Pair>& pairs,
        double rvdw,              // actual force cutoff (< rlist)
        double L,
        std::array<double,9>& W)  // virial tensor [xx,xy,xz,yx,yy,yz,zx,zy,zz]
{
    const double eps = 1.0, sig = 0.34;
    const double sig2 = sig*sig;
    const double rvdw2 = rvdw*rvdw;

    // Zero forces
    for (auto& pi : p) pi.force = {0,0,0};
    W.fill(0.0);

    for (const auto& pr : pairs) {
        int i = pr.i, j = pr.j;
        Vec3 dr = pbc({p[j].pos.x - p[i].pos.x,
                       p[j].pos.y - p[i].pos.y,
                       p[j].pos.z - p[i].pos.z}, L);
        double r2 = dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
        if (r2 >= rvdw2) continue;  // outside force cutoff (buffer region)

        // Soft-core floor: if r < 0.7σ, clamp to the force at 0.7σ.
        // This prevents numerical explosion when particles briefly overlap
        // (can happen between list rebuilds in a gas).  The artifact signal
        // is in pairs near rlist, not the core, so this doesn't suppress it.
        const double r_floor2 = 0.7*0.7*sig2;
        double r2_eff = std::max(r2, r_floor2);
        double ir2  = sig2 / r2_eff;
        double ir6  = ir2*ir2*ir2;
        double ir12 = ir6*ir6;
        // dU/dr^2 * 2 = 24ε/r^2 * (2(σ/r)^12 - (σ/r)^6)
        double fscal = 24.0 * eps / r2_eff * (2.0*ir12 - ir6);

        Vec3 fvec = {fscal*dr.x, fscal*dr.y, fscal*dr.z};

        p[i].force.x += fvec.x;
        p[i].force.y += fvec.y;
        p[i].force.z += fvec.z;
        p[j].force.x -= fvec.x;
        p[j].force.y -= fvec.y;
        p[j].force.z -= fvec.z;

        // Virial: W_αβ = Σ r_{ij,α} F_{ij,β}  (Newton's 3rd law symmetric)
        // dr points from i→j,  F on j from i = -fvec
        // Standard convention: W_αβ = r_{α} * F_{β}
        double* w = W.data();
        w[0] += dr.x * fvec.x;   // Wxx
        w[1] += dr.x * fvec.y;   // Wxy
        w[2] += dr.x * fvec.z;   // Wxz
        w[3] += dr.y * fvec.x;   // Wyx
        w[4] += dr.y * fvec.y;   // Wyy
        w[5] += dr.y * fvec.z;   // Wyz
        w[6] += dr.z * fvec.x;   // Wzx
        w[7] += dr.z * fvec.y;   // Wzy
        w[8] += dr.z * fvec.z;   // Wzz
    }
}

// ---------------------------------------------------------------------------
// Leapfrog integrator step
// ---------------------------------------------------------------------------
static void integrate(std::vector<Particle>& p, double dt, double mass) {
    for (auto& pi : p) {
        pi.vel.x += dt * pi.force.x / mass;
        pi.vel.y += dt * pi.force.y / mass;
        pi.vel.z += dt * pi.force.z / mass;
        pi.pos.x += dt * pi.vel.x;
        pi.pos.y += dt * pi.vel.y;
        pi.pos.z += dt * pi.vel.z;
    }
}

// ---------------------------------------------------------------------------
// Kinetic energy
// ---------------------------------------------------------------------------
static double kinetic(const std::vector<Particle>& p, double mass) {
    double ke = 0;
    for (const auto& pi : p)
        ke += 0.5 * mass * (pi.vel.x*pi.vel.x + pi.vel.y*pi.vel.y + pi.vel.z*pi.vel.z);
    return ke;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main(int argc, char** argv) {
    if (argc < 5) {
        std::fprintf(stderr,
            "Usage: %s <nstlist> <rlist_nm> <nsteps> <N_particles>\n"
            "  nstlist   : pair list update frequency (steps)\n"
            "  rlist_nm  : neighbor list cutoff (nm), e.g. 1.2\n"
            "  nsteps    : total MD steps\n"
            "  N_particles: number of LJ atoms\n"
            "Example: %s 20 1.2 5000 500\n",
            argv[0], argv[0]);
        return 1;
    }

    const int    nstlist = std::atoi(argv[1]);
    const double rlist   = std::atof(argv[2]);
    const int    nsteps  = std::atoi(argv[3]);
    const int    N       = std::atoi(argv[4]);

    // Force cutoff is rlist - 0.1 nm buffer (mirrors GROMACS default)
    const double rvdw    = rlist - 0.1;
    const double dt      = 0.001;   // 1 fs — halved for stability with hard cutoff
    const double mass    = 40.0;    // argon, g/mol  (amu in MD units)
    const double kB      = 0.008314; // kJ/(mol·K)
    const double T0      = 300.0;   // K target
    const double sig     = 0.34;    // LJ sigma (nm) — used for spacing

    // Place particles on a cubic lattice with spacing = 2.5σ — gas-phase density
    // (ρ* = Nσ³/V ≈ 0.05).  At 1.5σ particles quickly fall into the deep LJ well
    // and an integration timestep of 1 fs can't follow the core repulsion cleanly.
    const int    nx      = static_cast<int>(std::ceil(std::cbrt(N)));
    const double spacing = 2.5 * sig;         // 0.85 nm — well outside the attractive well
    double L = nx * spacing;
    if (L < rlist * 2.5) L = rlist * 2.5;    // ensure box > 2*rlist for PBC

    std::fprintf(stderr,
        "Box L=%.3f nm, spacing=%.3f nm (%.2fσ), rlist=%.3f nm, rvdw=%.3f nm, nstlist=%d\n",
        L, spacing, spacing/sig, rlist, rvdw, nstlist);

    // Initialize particles
    std::vector<Particle> particles(N);
    {
        int idx = 0;
        std::mt19937 rng(42);
        std::normal_distribution<double> vdist(0.0, std::sqrt(kB * T0 / mass));
        for (int ix = 0; ix < nx && idx < N; ++ix)
        for (int iy = 0; iy < nx && idx < N; ++iy)
        for (int iz = 0; iz < nx && idx < N; ++iz, ++idx) {
            particles[idx].pos = {(ix+0.5)*spacing, (iy+0.5)*spacing, (iz+0.5)*spacing};
            particles[idx].vel = {vdist(rng), vdist(rng), vdist(rng)};
            particles[idx].force = {0,0,0};
        }
        // Zero net momentum
        Vec3 vcm = {0,0,0};
        for (const auto& p : particles) { vcm.x += p.vel.x; vcm.y += p.vel.y; vcm.z += p.vel.z; }
        vcm.x /= N; vcm.y /= N; vcm.z /= N;
        for (auto& p : particles) { p.vel.x -= vcm.x; p.vel.y -= vcm.y; p.vel.z -= vcm.z; }
    }

    // CSV header
    std::printf("step,Pxx_bar,Pyy_bar,Pzz_bar,Wxx,Wyy,Wzz,KE,rebuilt\n");

    std::vector<Pair> pairlist;
    std::array<double,9> W;
    W.fill(0.0);
    const double V = L*L*L;
    // Conversion: virial in kJ/mol → pressure in bar
    // P = (2KE + W_trace/3) / (3V)  in kJ/mol/nm^3
    // 1 kJ/mol/nm^3 = 16.6054 bar
    const double kJpermol_per_nm3_to_bar = 16.6054;
    // Degrees of freedom: 3N - 3 (subtract COM)
    const int    ndof  = 3 * N - 3;
    // v-rescale every step — isokinetic thermostat, keeps KE exactly at T0.
    // The virial artifact is entirely in W_αβ (positions+forces), not KE,
    // so rescaling each step doesn't suppress the artifact; it isolates it.
    const int    nstrescale = 1;

    // Pre-equilibrate with small dt for 1000 steps to relax lattice forces
    {
        const double dt_eq = 0.0002;
        std::vector<Pair> pl_eq = build_pairlist(particles, rlist, L);
        std::array<double,9> W_eq; W_eq.fill(0.0);
        for (int s = 0; s < 1000; ++s) {
            eval_forces(particles, pl_eq, rvdw, L, W_eq);
            integrate(particles, dt_eq, mass);
            for (auto& pi : particles) {
                auto wrap = [&](double x) { return x - L * std::floor(x / L); };
                pi.pos.x = wrap(pi.pos.x); pi.pos.y = wrap(pi.pos.y); pi.pos.z = wrap(pi.pos.z);
            }
        }
        // Rescale to T0 after pre-eq
        double ke_eq = kinetic(particles, mass);
        double T_eq  = 2.0 * ke_eq / (ndof * kB);
        if (T_eq > 0) {
            double scale = std::sqrt(T0 / T_eq);
            for (auto& pi : particles) { pi.vel.x *= scale; pi.vel.y *= scale; pi.vel.z *= scale; }
        }
        std::fprintf(stderr, "Pre-equilibration done, T_final=%.1f K\n", T_eq);
    }

    for (int step = 0; step <= nsteps; ++step) {
        bool rebuilt = (step % nstlist == 0);
        if (rebuilt) {
            pairlist = build_pairlist(particles, rlist, L);
        }

        eval_forces(particles, pairlist, rvdw, L, W);

        double ke = kinetic(particles, mass);
        // Instantaneous pressure tensor (diagonal only for output):
        // P_αα = (2KE/3 + W_αα) / V  — simplified isotropic kinetic part
        double ke_contribution = 2.0 * ke / 3.0;
        double Pxx = (ke_contribution + W[0]) / V * kJpermol_per_nm3_to_bar;
        double Pyy = (ke_contribution + W[4]) / V * kJpermol_per_nm3_to_bar;
        double Pzz = (ke_contribution + W[8]) / V * kJpermol_per_nm3_to_bar;

        std::printf("%d,%.4f,%.4f,%.4f,%.6f,%.6f,%.6f,%.4f,%d\n",
                    step, Pxx, Pyy, Pzz, W[0], W[4], W[8], ke,
                    rebuilt ? 1 : 0);

        // Integrate (leapfrog)
        integrate(particles, dt, mass);

        // v-rescale thermostat: rescale velocities to T0 every nstrescale steps
        if (step % nstrescale == 0) {
            double ke_now = kinetic(particles, mass);
            double T_now  = 2.0 * ke_now / (ndof * kB);
            if (T_now > 0) {
                double scale = std::sqrt(T0 / T_now);
                for (auto& pi : particles) { pi.vel.x *= scale; pi.vel.y *= scale; pi.vel.z *= scale; }
            }
        }

        // Wrap positions into box
        for (auto& pi : particles) {
            auto wrap = [&](double x) { return x - L * std::floor(x / L); };
            pi.pos.x = wrap(pi.pos.x);
            pi.pos.y = wrap(pi.pos.y);
            pi.pos.z = wrap(pi.pos.z);
        }
    }

    return 0;
}
