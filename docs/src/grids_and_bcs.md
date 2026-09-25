# Grids and particle-surface interaction models

For non-0D simulations, Merzbild.jl provides a spatial grid, particle sorting onto that grid,
boundary conditions describing particle–surface interactions, and convection routines that move
particles and apply the boundary conditions.

## Grids

Currently a single 1-D grid type is implemented.

| Structure | Description | Type/constructor |
| --- | --- | --- |
| 1-D uniform grid | Uniform grid on ``[0, L]`` with `nx` cells | [`Grid1DUniform`](@ref) |

A [`Grid1DUniform`](@ref) discretizes the domain ``[0, L]`` into `nx` equal cells; it subtypes
[`AbstractGrid`](@ref). The
grid can be written to disk with [`write_grid`](@ref).

## Particle sorting


| Structure | Description | Type/constructor |
| --- | --- | --- |
| In-place grid sorter | Sorts particles into their grid cells | [`GridSortInPlace`](@ref) |

Particles are sorted into their cells with [`sort_particles!`](@ref) using a
[`GridSortInPlace`](@ref) instance, which also restores contiguity of the particle indexing (see
[Particle buffers and contiguous indexing](@ref "Particle buffers and contiguous indexing")). 

## Boundary conditions: 1D

Boundary conditions describe how particles interact with the walls at the ends of the 1-D domain.
They are collected into a tuple (left wall, right wall) that is passed to the convection routines.
All wall types subtype `AbstractBC`.

| Boundary condition | Description | Type |
| --- | --- | --- |
| Fully specular wall | Reflects the particle by flipping its normal velocity component | [`FullySpecularBC1D`](@ref) |
| Fully diffuse wall | Re-emits the particle from a Maxwellian at the wall temperature | [`FullyDiffuseBC1D`](@ref) |
| Maxwell wall | Mixture of specular and diffuse reflection set by an accommodation coefficient | [`MaxwellWallBC1D`](@ref) |
| Outflow | Deletes the particle: it leaves the domain | [`OutflowBC`](@ref) |
| Inflow | Reservoir boundary: emits particles into the domain, deletes those reaching it | [`InflowBC1D`](@ref) |

[`OutflowBC`](@ref) is neither species-specific nor dimension-specific and simply deletes particles
reaching the surface. Deletion makes `pia` non-contiguous so a sort is required before any per-cell
operation runs. Outflowing particles are recorded in [`SurfProps`](@ref) as incident only.

[`InflowBC1D`](@ref) describes an equilibrium reservoir of a given number density, temperature and
streaming velocity adjacent to the surface. It is species-specific, and particles of the domain
reaching it are deleted exactly as for an [`OutflowBC`](@ref). The velocity of an entering particle
is sampled with [`sample_inflow_velocity`](@ref) using the method described in [Garcia and Wagner (2006)](https://doi.org/10.1016/j.jcp.2006.01.025), and the number flux across the surface is given by
[`compute_inflow_number_flux`](@ref).

TODO: describe the emission routine that injects the particles into the domain once it is
implemented, including the per-particle residual flight time within a timestep.

Both [`FullyDiffuseBC1D`](@ref) and [`MaxwellWallBC1D`](@ref) are species-specific, since the
thermal reflection velocity depends on the species mass; [`MaxwellWallBC1D`](@ref) reduces to
specular reflection at an accommodation coefficient of 0 and to fully diffuse reflection at a value
of 1. These specular/diffuse gas–surface interaction models are the standard DSMC wall models found in
[Bird (1994)](https://doi.org/10.1093/oso/9780198561958.001.0001).

Periodic BCs for 1D grids are implemented via convection functions that wrap the particle positions,
see below.

## Convection

Convection moves particles according to their velocities and applies the boundary conditions when
particles reach a wall (periodic versions do not take any
lists of boundary conditions as input parameters and simply wrap
particle coordinates on 1D uniform grids). Six routines are available, distinguished by whether they accumulate
surface properties into a [`SurfProps`](@ref) instance, and by whether they also compute and store
each particle's post-convection cell index (into `particles.cell`, which avoids required an additional
pass over all particles during sorting), plus periodic convection versions:

| Function | Computes surface properties | Computes post-convection cell index |
| --- | --- | --- |
| [`convect_particles!`](@ref) | no | no |
| [`convect_particles_periodic!`](@ref) (1D uniform grid only) | no | no |
| [`convect_particles!`](@ref) (with `SurfProps`) | yes | no |
| [`convect_particles_and_compute_cell!`](@ref) | no | yes |
| [`convect_particles_and_compute_cell_periodic!`](@ref) (1D uniform grid only) | no | yes |
| [`convect_particles_and_compute_cell!`](@ref) (with `SurfProps`) | yes | yes |

The surface-property-computing variants take an additional [`SurfProps`](@ref) argument, into which
per-wall incident and reflected fluxes of mass, momentum, and energy are accumulated as particles
strike the walls.

The convection routines return the number of particles deleted at
[`OutflowBC`](@ref) and [`InflowBC1D`](@ref) boundaries (particles reaching an inflow/reservoir boundary are deleted);
the periodic convection routine always returns 0 to be consistent with non-periodic convection functions.
