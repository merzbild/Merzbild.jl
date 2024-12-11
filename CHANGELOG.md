# Changelog

## v0.5.0
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