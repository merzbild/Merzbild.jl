# Merging algorithms

Particle merging reduces the number of particles in a cell while (approximately) preserving the
distribution function and its moments. This is essential for variable-weight simulations, where the
particle count in a cell can grow without bound. All merging routines operate on the particles of a
single species in a single grid cell, addressed via the [`ParticleIndexerArray`](@ref).

Merging deletes particles and therefore introduces "holes" in the particle indexing; see
[Particle buffers and contiguous indexing](@ref "Particle buffers and contiguous indexing") for how
to restore contiguity with [`squash_pia!`](@ref) afterwards. Several of the routines come in two
variants: a plain one, and one that additionally takes a [`Grid1DUniform`](@ref) and clamps any
post-merge particle positions that fall outside the simulation domain back into it.

## Overview

| Algorithm | Setup struct | Function |
| --- | --- | --- |
| Velocity grid-based (N:2) | [`GridN2Merge`](@ref) | [`merge_grid_based!`](@ref) |
| Octree (N:1 and N:2) | [`OctreeMerge`](@ref) | [`merge_octree!`](@ref) |
| Roulette (+ conservative variant) | — | [`merge_roulette!`](@ref) |
| NNLS, moment-preserving | [`NNLSMerge`](@ref) | [`merge_nnls_based!`](@ref) |
| NNLS, approximate electron–neutral rate-preserving | [`NNLSMerge`](@ref) | [`merge_nnls_based_rate_preserving!`](@ref) |
| NNLS, exact electron–neutral rate-preserving | [`NNLSMerge`](@ref) | [`merge_nnls_based_rate_preserving!`](@ref) |

## Grid-based merging

[`merge_grid_based!`](@ref) groups particles using a Cartesian grid in velocity space (particles
outside the grid are grouped by velocity octant) and merges the particles in each cell/octant down
to 2 particles, following the algorithm of [Oblapenko et al. (2020)](https://doi.org/10.1016/j.jcp.2020.109302)
(see also [Vranic et al. (2015)](https://doi.org/10.1016/j.cpc.2015.01.020)). The velocity-space grid is described by a
[`GridN2Merge`](@ref) instance. The grid extent can either be derived from the species temperature
(via a `PhysProps` instance) or supplied explicitly, and each of these has a
[`Grid1DUniform`](@ref) overload that keeps post-merge particles inside the domain.

## Octree merging

[`merge_octree!`](@ref) performs octree-based N:M merging following
[Martin and Cambier (2016)](https://doi.org/10.1016/j.jcp.2016.01.020): velocity
space is recursively subdivided into octants, and the particles in each leaf bin are merged down to
`M` particles. The merge target `M` is a type parameter of the [`OctreeMerge`](@ref) struct, with
`M = 1` (N:1) and `M = 2` (N:2) supported. The construction of the octree is controlled by three
enums:

* [`Merzbild.OctreeBinSplit`](@ref) — how the splitting velocity is chosen
  ([`OctreeBinMidSplit`](@ref Merzbild.OctreeBinSplit), [`OctreeBinMeanSplit`](@ref Merzbild.OctreeBinSplit), [`OctreeBinMedianSplit`](@ref Merzbild.OctreeBinSplit)).
* [`Merzbild.OctreeInitBin`](@ref) — how the bounds of the initial bin are computed
  ([`OctreeInitBinMinMaxVel`](@ref Merzbild.OctreeInitBin), [`OctreeInitBinMinMaxVelSym`](@ref Merzbild.OctreeInitBin), [`OctreeInitBinC`](@ref Merzbild.OctreeInitBin)).
* [`Merzbild.OctreeBinBounds`](@ref) — how the bounds of a split sub-octant bin are computed
  ([`OctreeBinBoundsInherit`](@ref Merzbild.OctreeBinBounds), [`OctreeBinBoundsRecompute`](@ref Merzbild.OctreeBinBounds)).

A [`Grid1DUniform`](@ref) overload of [`merge_octree!`](@ref) is available that keeps post-merge
particles inside the domain.

## Roulette merging

[`merge_roulette!`](@ref) deletes random particles until the target particle count is reached and
re-weights the remaining particles to conserve number density, following
[Watrous, Seidel et al. (2023)](https://www.osti.gov/servlets/purl/2431184). Passing
the `conservative=true` keyword additionally corrects the post-merge velocities to conserve
momentum and energy.

## NNLS merging

The NNLS (non-negative least squares) merging routines cast moment preservation as a non-negative
least squares problem, following [Oblapenko and Torrilhon (2026)](https://doi.org/10.48550/arXiv.2604.00668).
The set of preserved moments and the
workspace are held in an [`NNLSMerge`](@ref) instance; [`compute_multi_index_moments`](@ref) is a
helper for building the multi-index moment set. All variants return `-1` on failure (residual too
large, or no reduction in particle count achieved).

* [`merge_nnls_based!`](@ref) — general moment-preserving merging.
* [`merge_nnls_based_rate_preserving!`](@ref) — merging of electrons that additionally preserves
  electron–neutral elastic scattering and electron-impact ionization rates. Two methods are
  provided: one preserving **approximate** rates, and one (taking an additional neutral collision
  partner `ParticleVector`) preserving **exact** rates for a specific neutral species.
