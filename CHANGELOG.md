# Changelog

## v0.8.3
* `sort_particles!` and `sort_particles_after_exchange!` record the first and last cell holding particles
in the `occ_lo`/`occ_hi` fields of the
`GridSortInPlace` instance, for use by `update_occupancy_bounds!`
* Added electrostatic Particle-in-Cell on a 1-D uniform grid: charge deposition (`deposit_charge!`,
`normalize_charge_density!`), a Poisson solve (`PoissonSolver1DUniform`, `solve_poisson!`), and the gather + push
(`accelerate_electric_field_x!`); `ElectrostaticFieldProps` struct to hold nodal charge density, potential, and electric field
* Field boundary conditions `DirichletFieldBC1D`, `NeumannFieldBC1D`, and `PeriodicFieldBC1D`
* NetCDF output of the electrostatic field quantities stored in an `ElectrostaticFieldProps` instance
via `NCDataHolderField` and `write_netcdf`
* Added `debye_length` and `plasma_frequency` functions, and the vacuum permittivity `eps_0`
* `reduce_field_props!` for reduction after thread-safe charge deposition and 
* Elastic scattering models are no longer hard-coded to VHS: the VSS model (`VSS`) has been added alongside
VHS (`VHS`), and both are supported by `ntc!`, `ntc_equal_weight!`, and `swpm!` (see new bundled interaction data file `vss.toml`)
* The elastic scattering model is fixed per species pair and is stored in the `Interaction` instance as a
`ScatteringModel` enum value (`ScatteringVHS`, `ScatteringVSS`), set by the optional `model` key
of the species pair's entry in the interaction data TOML file (defaults to `"VHS"`); the VSS model additionally
requires a `vss_alpha` key
* `collide_2particles_vhs!` and `collide_2particles_vhs_equal_weight!` renamed to `collide_2particles!` and
`collide_2particles_equal_weight!`, and now take the scattering model tag as their second argument
* Improved test coverage

## v0.8.2
* NNLS merging speed-ups
* Batched particle deletion for more efficient deletion of particles in merging routines
* Use unrolled function for norm of 3-vectors
* Faster exponentiation in computation of VHS cross-sections
* Load re-balancing via `LoadBalancerCellQ`
* Speed-ups in particle exchange and re-sorting in multi-threaded simulations
* `ChunkExchanger` now stores its indexing as a structure of arrays (`start1`, `n_group1`, `start2`, `n_group2`)
instead of an array of `ParticleIndexer` instances; the group ends are no longer stored, as they are given by
`start + n_group - 1`
* Added `update_occupancy_bounds!`, which records the first and last cell in which a chunk holds particles,
so that `exchange_particles!` can reject a pair of chunks with nothing to exchange without scanning any cells
* Better documentation through use of `DocumenterCodeBlocks.jl`
* Documentation on load balancing
* Statistical tests for isotropic (VHS) scattering

## v0.8.1
* Clean-up of roulette merging code and minor optimizations
* Documentation improvements
* Minor optimizations in computation of fluxes on grid
* Minor optimizations in computation of collision cross-sections (`1 - 2 vhs_o` factor used to compute VHS cross-sections now pre-computed
and stored for each interaction)
* Added simplified `Interaction` constructor
* NNLS merging speed-ups and code clean-up

## v0.8.0

### Breaking changes
* Usage of `Vector{Particle}` completely deprecated; tests have been updated
* `Particle` type now has a `D` type parameter describing the dimension of the position vector. `D` defaults to 3 unless
specified. `ParticleVector` now also has a `D` type parameter.
* Merging routines now also have a `D` type parameter describing the dimension of the position vector of the particles to merge.
`D` defaults to 3 unless specified.
* Octree merging struct renamed to `OctreeMerge` and also has a second type parameter `M`: the number of post-merge particles in a bin (refactoring
by Yanliang Zhu)
* Removed unused `update_particle_indexer_new_lower_count!` function
* NNLS merging functions now do not take `v_multipliers`, `n_rand_pairs`, and `centered_at_mean` as arguments, and do not use any fictitious particles
* The function `check_speed_bounds` has been removed
* NetCDF output now does not write species' names to a variable, but rather to a global attribute as a single comma-separated string (implementation by Yanliang Zhu)
* NetCDF PhysProps output now does not write total moments; an `NCDataHolderMoments` has been added for that purpose (implementation by Yanliang Zhu)
* PhysProps no longer stores any moment data and `compute_props_with_total_moments!` has been removed (implementation by Yanliang Zhu)
* Particle netCDF I/O simplified, skips output of position and cell data for 0-D particles
* `load_electron_neutral_interactions` replaced by an `ElectronNeutralInteractions` constructor
* `MaxwellWalls1D` removed
* Boundary conditions now need to be passed to `convect_particles!` or `convect_particles_and_compute_cell!` as a `Tuple` of boundary conditions, one for each surface.
* Order of parameters passed to `fp_linear!` and `mean_collision_frequency` changed to be consistent with guidelines in `CONTRIBUTING.md`.
* `merge_octree_N2_based!` renamed to `merge_octree!`

### New functionality
* `squash_pia!` can be called less frequently, i.e. only before particle convection, as better tracking of last index
in particle arrays has been implemented. See documentaton on contiguous indexing.
* NetCDF output now allows to choose exact format and defalts to `NC_64BIT_OFFSET`, significantly speeding-up I/O (implementation by Yanliang Zhu)
* `compute_moment_scaling!`, `compute_moments!` functions to compute scaled total velocity moments (implementation by Yanliang Zhu)
* `NCDataHolderMoments` struct type added for I/O of computed moment data (implementation by Yanliang Zhu)
* More specialized boundary conditions added for 1D simulations: `FullyDiffuseBC1D`, `MaxwellWallBC1D`, `SpecularWallBC1D`
* The `ParticleIndexerArray` type now directly stores number of cells and species it is tracking (`.n_cells`, `.n_species`)
* `MERZBILD_DATA_PATH` now exported for easier loading of particle and interaction data bundled with Merzbild.jl
* `MERZBILD_SIMULATIONS_PATH` now exported to expose path to example simulations bundled with Merzbild.jl
* `MERZBILD_SCRIPTS_PATH` now exported to expose path to Python scripts bundled with Merzbild.jl
* Octree-based merging now supports N:M merging in each bin (currently, only conservative N:1 and N:2 merging implemented; implementation by Yanliang Zhu)

### Misc
* Added tests that check that for unexpected memory allocations
* Various performance optimizations
* Minor bug fixes
* Fixed linear Fokker-Planck for variable-weight particles
* Documentation improvements

## v0.7.10
* Documentation improvements
* Improved test coverage
* `convect_particles_and_compute_cell!` function added that computes new cell indices during the convection, a new version of `sort_particles!`
that takes advantage of these pre-computed cell indices has also been added
* `squash_pia!` now also squashes particle cell indices

## v0.7.9
* `restore_particle_ordering!` added, this restores optimal indexing of particles and can lead to simulation speed-ups due to improved cache usage
* Minor optimizations in octree merging
* Simplified elastic VHS collision code via a unified `collide_2particles_vhs!` function, minor speed-ups
* New keyword parameter `dw_tol` in collision routines that sets tolerance in weight difference under which particle
collisions are treated as equal weight collisions and no additional particle splitting is performed
* Minor improvements in readability and speed of electron-neutral collision routines
* `@inbounds` added to acceleration routine
* `ntc_equal_weight` added for slightly faster collisions in equal-weight simulations
* Documentation improvements
* Improved test coverage

## v0.7.8
* Fixed ionization simulations without event splitting
* Tests for simulations with elastic and ionizing electron-neutral collisions (without event splitting)
* Threaded swap of particles for multi-threaded simulations now possible
* New `sample_particles_phase_box_weighted!` sampling function added for variable-weight particle sampling
from an equilibrium distribution
* Variable-weight simulation examples now perform merging before the time loop is started, and
collision factor estimates use mean post-merge particle weight
* `run_examples.py` now runs each example only for 10 timesteps
* Documentation of simulation files in `simulations` directory added (`simulations/README.md`)
* Exact rate-conserving version of NNLS added
* Simulation examples added for Fourier flow
* Simulation example added that simply samples and merges particles and computes various statistics
* Reproducibility setups added for "Moment-preserving particle merging via non-negative least squares" (simulation parameters and
post-processing scripts)

## v0.7.7
* Improved scattering modelling for ionization reactions based on Nanbu's paper
* Added tests for ionization reaction scattering
* Documentation improvements (new sections on ionization modelling and on debugging)
* `check_unique_buffer` utility function for debugging

## v0.7.6
* Speed-up of variable weight SWPM collisions
* Rate-preserving merging now makes use of reference process cross-sections
* Fix in `simulations/OD/BKW/bkw.jl`
* Clean-up of some `@inbounds` macros in sampling routines
* Option to set merging grid extent explicitly when calling `merge_grid_based!`
* Improved test coverage

## v0.7.5
* Octree and grid-based merging can now handle bins/cells with only zero-weight particles inside
* Grid-based merging 1D version added

## v0.7.4
* Fixed rate-preserving NNLS merging, it can be considered more or less stable now, but only for the 0D case,
as preservation of spatial moments is not implemented
* NNLS merging improvements (store original, non-scaled velocities and positions separately)
* Improved test coverage
* Slight re-ordering of `pretty_print_pia` output
* `conservative` keyword added to `merge_roulette!` merging algorithm to ensure conservation of momentum and energy

## v0.7.3
* Fixes in particle exchange for multi-threaded computations
* Python plotting script for 1-D flows now allows for multiple files + plot labels
* `physical_props.jl`, `surface_props.jl`, `flux_props.jl` files moved to new `src/properties` directory
* Added computation of mean collision frequency and mean free path for single-species VHS gases: `mean_free_path`,
`mean_collision_frequency`
* Speed-up of computation of `FluxProps` (unnecessary allocations removed)
* Reduced allocation in computation of `PhysProps` for sorted particles

## v0.7.2
* Added `w_threshold` keyword argument to NNLS merging to discard particles with very small weights
* NNLS merging now also accepts a list of spatial moments to preserve (requires explicitly conserving 1-st order spatial moments, otherwise
the corresponding coordinates are set to 0); not implemented for rate-preserving merging
* 1D example added for NNLS merging (`simulations/1D/couette_varweight_nnls.jl`)
* NNLS merging speed-ups and stability improvements
* New weighted median with interpolation function for octree merging
* `FluxProps` structure added, along with `compute_flux_props!`, `compute_flux_props_sorted!` functions to compute fluxes in grid cells
* I/O for `FluxProps` added
* Documentation improvements
* Improved test coverage
* Single-species Stochastic Weighted Particle Method (SWPM) implemented (`swpm!` function)
* Tabulated cross-section computations and associated NTC routines for electron-neutral interactions now correctly deal with
how the out-of-table values should be handled
* Speed-up of variable weight collisions
* Creation of a `ParticleVector(n_particles)` now leads to it being filled with particles with weight, velocity, position 0
* Resizing of a `ParticleVector` also leads to new particles being instantiated immediately and filled with weight, velocity, position 0

## v0.7.1
* More use of Muladd in a few places
* Roulette merge added (`merge_roulette!`)
* Documentation improvements
* CI improvements
* `update_particle_indexer_new_lower_count!`, `delete_particle_end_group2!` now set `end2` to -1 if last particle from group2 is deleted
* Now possible to write particle data via `write_netcdf`

## v0.7.0
* Added use of MuladdMacro (via `@muladd`), updated tolerances in tests and reference solutions
* Removed unused functions `sample_maxwellian_single!`, `sample_maxwellian!(rng, particles, nparticles, m, T, v0)`
* Removed `create_vdf`, `create_unit_dvgrid`, `create_noiseless_dvgrid`, replaced with constructors
* `update_particle_indexer_new_particle` renamed to `update_particle_indexer_new_particle!`
* `update_particle_indexer_new_lower_count` renamed to `update_particle_indexer_new_lower_count!`
* `update_particle_buffer_new_particle` renamed to `update_particle_buffer_new_particle!`
* `write_netcdf_surf_props`, `write_netcdf_phys_props` replaced with `write_netcdf`
* Fokker-Planck speed-up via use of pre-allocated arrays to store sampled velocities
* `fp!` renamed to `fp_linear!`
* `count_disordered_particles` added
* `check_pia_is_correct` and `check_unique_index` functions added for diagnostics
* Multithreaded simulations now possible (currently not verified for variable-weight simulations)
* NNLS merging now accepts additional keyword arguments `centered_at_mean`, `v_multipliers`, `iteration_mult`, `scaling`
    and defaults to a new scaling algorithm
* Documentation improvements
* Improved test coverage

## v0.6.6
* Utility function `Merzbild.add_particle!` added
* Improved test coverage
* Internal code improvements (use of `inbounds`, minor speed-up of `compute_props_sorted!`)
* `Merzbild.compute_octree!` function added to reduce code duplication
* `Merzbild.bin_bounds_recompute!` function fixed, as it could give erroneous bounds
* `squash_pia!` now correctly deals with empty cells
* Grid sorting logic simplified and now simply calls `squash_pia`! if non-contiguous indices
* No more memory allocations in `avg_props!`

## v0.6.5
* Possible to compute number density in cell when computing calling `compute_props_sorted!` and passing in a grid
* Documentation improvements
* Improved test coverage

## v0.6.4
* Now possible to compute surface properties due to particle-surface interactions.
* Unified `update_particle_indexer_new_particle` and `update_particle_buffer_new_particle`
    into single routine `update_buffer_index_new_particle!`.
* Documentation improvements
* Fixed bug where `compute_props_sorted!` did not compute the number of particles
* `pretty_print_pia` added
* `AbstractNCDataHolder` abstract type added, `NCDataHolder` and `NCDataHolderSurf` are subtypes thereof
* Example for 1-D variable weight DSMC simulation added

## v0.6.3
* Documentation improvements
* Function to restore continuity to indices in a `ParticleVector`/`ParticleIndexerArray` pair

## v0.6.2
* Fixes in resizing of `ParticleVector` instances.
* Tests for buffers for `ParticleVector` instances.
* `delete_particle!`, `delete_particle_end!`, `delete_particle_end_group1!`, and `delete_particle_end_group2!`
    internal functions added for particle deletion.
* Merging routines now update particle buffers upon particle count reduction.
* `contiguous` field added to `ParticleIndexerArray` structure to keep track of lack of gaps in indexing. Convection
    and sorting routines now can deal with gaps in indexing; sorting removes these gaps.
* Documentation improvements.
* Performance improvements in convection and sorting for contiguous `ParticleIndexerArray` instances.

## v0.6.1
* Utility function `load_species_and_interaction_data` added.

## v0.6.0
* The Octree and Grid-based merging routines now require the random number generator to be passed.
* Tests now rely on the `StableRNGs` package to avoid changes in RNGs between julia versions. BKW sampling fixed and
now uses the RNG. CI Github workflow added. Test values and tolerances updated. Examples in the `simulations` directory updated.

## v0.5.0
* First draft of documentation.
* Removed one old version of particle sampling routine used in some of the tests and examples.
* `resize!` now works with `ParticleVector` instances (preliminary work for variable-weight DSMC in non-0D settings).
* `buffer` field added to `ParticleVector` instances to keep track of unused particles. Additional tests.
* Fixed simple estimator of (sigma_w_g)_max for the multi-species case.
* `MaxwellWalls` renamed to `MaxwellWalls1D`.

## v0.4.1
* Implemented linear Fokker-Planck model for a single species gas without internal degrees of freedom (implementation by Leo Basov).

## v0.4.0
* Proper constructors added for a lot of the structs used in the code; old initialization functions removed.
* Other additions include code coverage reports via Coverage.jl and more tests.
* `compute_props!` renamed to `compute_props_with_total_moments!`, `compute_props!` now does not compute the moments
    and computes energy faster. 
* `compute_props_sorted_without_moments!` now renamed to `compute_props_sorted!`.
* Can now pass a list of physical property (i.e. `["T", "v"]`) names to I/O to skip in output.

## v0.3.0
* 1-D simulations on uniform grids now possible.
* Package name in Project.toml is now fixed (was "merzbild", now is "Merzbild"), and tests don't need to include the src files anymore.
* Aqua.jl testing added.

## v0.2.0
* Consistent argument order in all functions. Clean-up of simulation examples.

## v0.1.0
* This is the first development version, with support for 0-D simulations and some basic plasma processes implemented.