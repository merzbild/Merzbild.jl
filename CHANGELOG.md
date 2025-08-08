# Changelog

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
* NNLS merging now accepts additional keyword arguments `centered_at_mean`, `v_multipliers`, `iteration_mult`
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
* Implemented linear Fokker-Planck model for a single species gas without internal degrees of freedom.

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