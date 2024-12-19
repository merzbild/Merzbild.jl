# Changelog

## v0.6.0
The Octree and Grid-based merging routines now require the random number generator to be passed.
Tests now rely on the `StableRNGs` package to avoid changes in RNGs between julia versions. BKW sampling fixed and
now uses the RNG. CI Github workflow added. Test values and tolerances updated. Examples in the `simulations` directory updated.

## v0.5.0
First draft of documentation. Removed one old version of particle sampling routine used in some of the tests and examples.
`resize!` now works with `ParticleVector` instances (preliminary work for variable-weight DSMC in non-0D settings).
`buffer` field added to `ParticleVector` instances to keep track of unused particles. Additional tests.
Fixed simple estimator of (sigma_w_g)_max for the multi-species case.
`MaxwellWalls` renamed to `MaxwellWalls1D`.

## v0.4.1
Implemented linear Fokker-Planck model for a single species gas without internal degrees of freedom.

## v0.4.0
Proper constructors added for a lot of the structs used in the code; old initialization functions removed.
Other additions include code coverage reports via Coverage.jl and more tests.
`compute_props!` renamed to `compute_props_with_total_moments!`, `compute_props!` now does not compute the moments
and computes energy faster. 
`compute_props_sorted_without_moments!` now renamed to `compute_props_sorted!`.
Can now pass a list of physical property (i.e. `["T", "v"]`) names to I/O to skip in output.

## v0.3.0
1-D simulations on uniform grids now possible. Package name in Project.toml is now fixed (was "merzbild", now is "Merzbild"), and tests don't need to include the src files anymore.
Aqua.jl testing added.

## v0.2.0
Consistent argument order in all functions. Clean-up of simulation examples.

## v0.1.0
This is the first development version, with support for 0-D simulations and some basic plasma processes implemented.