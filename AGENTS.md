# Merzbild.jl AGENTS.md file

Merzbild.jl is a variable-weight **DSMC + linear Fokker–Planck** code for rarefied gas dynamics,
written in Julia. It is a **library**: it provides the building blocks (particle sampling/indexing,
collisions, merging, property computation, I/O) and the user writes their own time loop. Supports
0D/1D, serial and multithreaded operation.

## Repository layout

* Source is in `src/`; entry point `src/Merzbild.jl` (module, `include`s, exports). Feature
  subdirectories, each fronted by a same-named glue include file:
  * `collisions/` — NTC DSMC, SWPM, Fokker–Planck, cross-sections, scattering
  * `grids/` — 1D uniform grid, particle sorting, grid I/O
  * `merging/` — grid / octree / roulette / NNLS particle merging
  * `properties/` — physical, moment, flux, collisional, surface properties
  * `convection/` — boundary conditions, 1D convection
  * `pic/` — constant-field acceleration
  Top-level: `particles.jl`, `distributions_and_sampling.jl`, `io.jl`, `io_moments.jl`,
  `parallel.jl`, `constants.jl`, `utils.jl`, `abstract_types.jl`.
* Project documentation apart from docstrings is in `docs/src` (Documenter.jl). Particle indexing is
  documented in `docs/src/overview_blocks.md` and `docs/src/contiguous_indexing.md`.
* Example simulations are in `simulations/0D` and `simulations/1D`.
* **Coding conventions** (naming, argument order, contribution checklist) are in `CONTRIBUTING.md` —
  follow these when writing or changing code.

## General instructions and tips

* Ignore all directories and files listed in `.gitignore` (e.g. `scratch/`, `vibe/`, `plots/`).
* Memory allocations are to be avoided at all costs (unless new instances of structs are
  instantiated or completely new particles added); prefer explicit loops and re-using data to
  cleaner-looking operations using vectorized syntax.
* Assume vector indices are never out-of-bounds unless writing debugging functions, use `@inbounds`
  where appropriate.
* NEVER modify `Project.toml` and `Manifest.toml` by yourself.
* Only include comments when absolutely necessary. When the function name or implementation clearly
  indicates its purpose or behavior, redundant comments are unnecessary.
* Naming: a `ParticleIndexer` variable is called `particle_indexer` (since `pi` could mean a
  particle index), but a `ParticleIndexerArray` variable is called `pia`. See `CONTRIBUTING.md` for
  the full argument-order rules.

## Testing instructions

* Run all project tests with `julia --project=. -e 'using Pkg; Pkg.test()'`.
* The suite is large and organized into grouped runners `test/runtests_*.jl` (basics, computes,
  collisions, 1D, threading, merging, FP, io, swpm, reference solutions, malloc) included by
  `test/runtests.jl`. Individual `test/test_*.jl` files are not standalone — they assume the `using`
  imports from `runtests.jl` and are `include`d into that scope.
* Unexpected-allocation regressions are guarded by `test/runtests_malloc.jl` /
  `test/test_malloc_*.jl`; keep these green when touching hot paths.
* New features must have corresponding tests and documentation added (see the checklist in
  `CONTRIBUTING.md`).
* `run_examples.py` runs all examples in `simulations/` for 10 timesteps/ensembles and writes errors to directory
  given by the `logdir` command line argument (defaults to `scratch/logs`). Use to verify simulation examples work after
  changes to code.