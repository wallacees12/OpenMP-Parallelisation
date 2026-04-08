# N-Body Gravitational Simulation — Serial vs OpenMP vs Pthreads

A C implementation of a 2D gravitational N-body simulation (`galsim`) parallelised three ways. The goal of the project is to measure how different parallelisation strategies scale with particle count and thread count, and where each strategy wins or loses.

Written for the **High Performance Programming (Assignment 4)** course at Uppsala University / University of Zurich.

## What it does

Given `N` point particles with initial positions, masses and velocities, `galsim` integrates their motion under pairwise Newtonian gravity using a direct `O(N²)` force sum and a simple symplectic Euler step:

```
a_i = G · Σ_{j ≠ i}  m_j · (x_j − x_i) / (‖x_j − x_i‖ + ε)³
v_i ← v_i + a_i · Δt
x_i ← x_i + v_i · Δt
```

The hot loop is the pairwise force computation. All three variants produce bit-identical particle trajectories (verified against `ref_output_data/`); they differ only in how the force sum is distributed across threads.

## Repository layout

```
Vanilla/      serial baseline (galsim.c)
Openmp/       OpenMP variant (galsim.c, uses #pragma omp parallel for)
Pthreads/     Pthreads variants:
                galsim_unbalanced.c     — naive row split
                galsim_balanced.c       — work-balanced row split
                galsim_balanced_row.c   — per-row fine-grained split
compare_gal_files/   binary diff tool for result.gal files
input_data/          ellipse_N_*.gal — binary particle initial conditions
ref_output_data/     reference outputs for correctness checking
run_galsim.sh        OpenMP timing sweep (N × threads)
run_and_compare.sh   build + run + diff against reference
Sam benchmark.txt    full timing results (this is the one with real numbers)
```

## Build

Each variant has its own Makefile. From the repo root:

```bash
# Serial baseline
cd Vanilla   && mv Makefile.txt Makefile && make && cd ..

# OpenMP
cd Openmp    && mv Makefile.txt Makefile && make && cd ..

# Pthreads (builds galsim_balanced and galsim_unbalanced)
cd Pthreads  && make && cd ..
```

All builds use `gcc -O3 -ftree-vectorize -march=native -ffast-math` (the Vanilla variant additionally uses `-Ofast -funroll-loops`). A C99 compiler with OpenMP and pthreads support is required.

## Run

```bash
# Serial
./Vanilla/galsim  N input_file nsteps dt graphics

# OpenMP
./Openmp/galsim   N input_file nsteps dt graphics n_threads

# Pthreads
./Pthreads/galsim_balanced  N input_file nsteps dt graphics n_threads
```

- `N` — number of particles
- `input_file` — binary `.gal` file (see `input_data/`)
- `nsteps` — integration steps
- `dt` — timestep (e.g. `1e-5`)
- `graphics` — `0` (disabled; graphics stripped for benchmarking)
- `n_threads` — worker thread count

Each run writes the final state to `result.gal` (or `result_balanced_row.gal` for the row-balanced Pthreads variant). Use `compare_gal_files` to diff against the references in `ref_output_data/`.

Example:
```bash
./Openmp/galsim 3000 input_data/ellipse_N_03000.gal 100 1e-5 0 8
./compare_gal_files/compare_gal_files 3000 result.gal ref_output_data/ellipse_N_03000_after100steps.gal
```

## Key design decisions

**Per-thread local acceleration buffers instead of mutex-locked updates.** The naive parallelisation of an N-body force sum has a write hazard: every `(i, j)` pair updates both `a[i]` and `a[j]`, so two threads working on different particle pairs can collide on the same acceleration entry. The first attempt used `pthread_mutex_lock` around each update and ran *slower* than the serial baseline due to lock contention. The final design gives every thread a private `2*N` acceleration buffer, lets it accumulate locally with no synchronisation, and merges into the global `a[]` once at the end of each timestep. Same result, no locks on the hot path.

**Newton's third law used to halve the pair count.** The inner loop iterates `j = i+1..N` and applies `+F` to particle `i` and `−F` to particle `j` in a single step, cutting force evaluations from `N²` to `N(N−1)/2`.

**Row-balanced Pthreads variant.** A naive row split gives thread 0 the outer-loop rows with the most work (`j` ranges from `1` to `N−1`) and thread `p−1` almost nothing. `galsim_balanced_row.c` instead splits each row's `j` range across all threads so work is evenly distributed per iteration.

**Restrict-qualified mass and position pointers.** `m` and `x` are passed to `Compute_ax_ay` with `__restrict` to let the compiler assume no aliasing between them and vectorise the inner arithmetic.

**Contiguous SoA layout.** Particle data is stored as `[m0..m_{N-1}][x0 y0 x1 y1 ..][vx0 vy0 ..][L0 L1 ..]` for cache-friendly access in the force loop.

## Benchmark results

All timings in seconds. Same machine, `-O3 -march=native -ffast-math`. Source: `Sam benchmark.txt`.

### N = 3000, 100 steps

| Threads | Serial | OpenMP | Pthreads |
| ------: | -----: | -----: | -------: |
|       1 |   9.73 |  10.25 |     7.10 |
|       2 |      — |   5.24 |     5.34 |
|       4 |      — |   2.69 |     3.16 |
|       8 |      — |   1.43 |     1.75 |
|      12 |      — |   1.39 |     1.36 |
|      14 |      — |   1.36 |     1.33 |

### N = 3000, 400 steps (final test)

| Threads | Serial | OpenMP | Pthreads |
| ------: | -----: | -----: | -------: |
|       1 |  39.08 |  41.21 |    28.32 |
|       2 |      — |  20.95 |    21.44 |
|       4 |      — |  10.70 |    12.64 |
|       8 |      — |   5.69 |     6.94 |
|      12 |      — |   5.39 |     5.31 |
|      14 |      — |   5.42 |     5.10 |

### Observations

- **Best speedup at N=3000, 400 steps: ~7.7× (Pthreads, 14 threads)** and ~7.6× (OpenMP, 12 threads) over the serial baseline.
- OpenMP has higher single-thread overhead than the serial baseline (`41.21s` vs `39.08s`) — the dynamic scheduler and per-thread buffer allocation aren't free — but scales more cleanly up to ~8 threads.
- Pthreads is faster single-threaded (`28.32s`) because the manually-tuned row split has less per-iteration overhead than OpenMP's runtime.
- Both converge to ~5.1–5.4s beyond 10 threads. The diminishing returns match the benchmark machine's physical core count; additional SMT threads give little benefit.
- Correctness was verified at every data point: `Pthreads` matches the reference `pos_maxdiff = 0.055774475988` exactly (the OpenMP variant shows `0` because its result filename differs; it was diffed separately and also matches).

Full sweep across `N ∈ {1000, 3000, 4000..10000}` and `threads ∈ {1..14}` is in `Sam benchmark.txt` and `timing_results.txt`.

## References

- Assignment spec: `Assignment4.pdf`, `High_Performance_Programming___A4.pdf`
- Background reading: `An Introduction to Parallel Programming.pdf` (Pacheco)
