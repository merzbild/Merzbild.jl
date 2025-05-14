# 1D DSMC simulations

In this section, setting up DSMC simulations on a uniform 1D grid will be discussed. Some of the
concepts and algorithms used here are also applicable to other grids (upcoming), but
some things are specific to the uniform 1D grids.

## Creating a grid
Creating a 1D uniform grid is very easy, as one needs to specify only the domain length `L` and the number
of cells `nx`:
```julia
grid = Grid1DUniform(L, nx)
```
The Grid1DUniform structure stores additional properties such as cell volume, required for collisions.

We can also immediately create the `ParticleIndexerArray` instance (assuming we have a single-species flow):
```julia
pia = ParticleIndexerArray(grid.n_cells, 1)
```

The grid information can be written to a NetCDF file by calling `write_grid`:
```julia
write_grid("grid_info.nc", grid)
```

## Sampling particles in each cell
To sample equal-weight particles in each cell, we need to compute the required `Fnum`.
If we want to have `ppc` particles and a number density of `ndens`, then we can compute `Fnum` as
(since the grid is uniform, we can use the volume of any cell for the computation)
```julia
Fnum = grid.cells[1].V * ndens / ppc
```

We then initialize an `Vector` of `ParticleVector`'s to store our particles for each species:
```julia
particles = [ParticleVector(n_particles)]
```
and perform the sampling:
```julia
sample_particles_equal_weight!(rng, grid, particles[1], pia, 1,
                               species_data, ndens, T, Fnum)
```

## Particle sorting
We will need to sort particles on the grid in case we will be convecting them; therefore
a structure for particle sorting needs to be created as well. An in-place bin sorting algorithm is used,
it requires an estimate of the number of particles in the simulation to pre-allocate arrays.
For a fixed-weight simulation where `ppc` (particles per cell) are sampled in each cell at the start
and no particles are created during the course of the simulation, we can compute the estimate
as `n_particles = ppc * nx`. So we can initialize the `GridSortInPlace` instance like this:
```julia
n_particles = ppc * nx

gridsorter = GridSortInPlace(grid, n_particles)
```

If we want to sort our particles (held in a `ParticleVector` instance), we can
simply call

```julia
species_id = 1

# we assume that particles has type Vector{ParticleVector} (a ParticleVector per species)
sort_particles!(gridsorter, grid, particles[species_id], pia, species_id)
```

## Creating boundary conditions
Next, we need to create boundary conditions for the left and right walls.
Currently, a diffusely reflecting wall is implemented with a user-defined accommodation coefficient
(if it is equal to 0, the reflection is fully specular; if it is equal to 1, the reflection is fully diffuse).

This type of boundary condition (which stores the wall temperature, wall velocity, and accommodation coefficient)
is described by the `MaxwellWallBC` structure, which should be defined for each wall.

However, it is not intended to be defined or used directly; instead, for a 1-D simulation, a higher-level
`MaxwellWalls1D` structure is used, which holds not only the two `MaxwellWallBC` instances (for the left and right walls),
but also some species-wise precomputed quantities for the diffuse reflection. The `MaxwellWalls1D` struct assumes that
the wall velocity in the `x` and `z` directions is 0, but a non-zero `y` velocity may be specified.

A `MaxwellWalls1D` instance can be initialized like this (this will create two walls with equal temperatures
and `y` velocities in opposite directions):
```julia

boundaries = MaxwellWalls1D(species_data, T_wall, T_wall, -v_wall, v_wall, 1.0, 1.0)
```

## Performing convection
Having set up the grid and boundary conditions, we can convect particles.
This is done by calling the `convect_particles!` function.
The convection should be followed by particle sorting before any computations of physical properties are done.

```julia
convect_particles!(rng, grid, boundaries, particles[species_id], pia, species_id, species_data, Δt)
```

## Bringing it all together
Now we can combine all the pieces to set up a simulation of a single-species Couette flow in a channel
with a width of 0.5 mm, discretized with 50 cells. The y-velocity of the left wall is assumed to be -500 m/s,
and that of the right wall 500 m/s; the temperature of both walls is 300 K. The solution is initialized with 100
particles per cell and a number density of 5e22 1/m^3. A timestep of 2.59 ns is used.
The simulation runs for 50K steps and the solution time-averaged is averaged after the first 14K steps.

```julia
using Merzbild
using Random

# set our random seed for reproducibility
seed = 1
Random.seed!(seed)
rng = Xoshiro(seed)

# set physical and discretization parameters
T_wall = 300.0
v_wall = 500.0
L = 5e-4
ndens = 5e22
nx = 50
ppc = 100
Δt = 2.59e-9
n_timesteps = 50000
avg_start = 14000

# load particle and interaction data
particles_data_path = joinpath("data", "particles.toml")
species_data = load_species_data(particles_data_path, "Ar")
interaction_data_path = joinpath("data", "vhs.toml")
interaction_data::Array{Interaction, 2} = load_interaction_data(interaction_data_path, species_data)

# create our grid and BCs
grid = Grid1DUniform(L, nx)
boundaries = MaxwellWalls1D(species_data, T_wall, T_wall, -v_wall, v_wall, 1.0, 1.0)

# init particle vector, particle indexer, grid particle sorter
# we will not be creating or destroying any particles, so we can compute the exact number
# of particles we will have in the simulation
n_particles = ppc * nx
particles = [ParticleVector(n_particles)]
pia = ParticleIndexerArray(grid.n_cells, 1)
gridsorter = GridSortInPlace(grid, n_particles)

# sample particles
# Fnum * ppc = Np in cell = ndens * V_cell
Fnum = grid.cells[1].V * ndens / ppc
sample_particles_equal_weight!(rng, grid, particles[1], pia, 1,
                                species_data, ndens, T_wall, Fnum)

# create collision structs
collision_data = CollisionData()

# create struct for computation of physical properties
phys_props = PhysProps(pia)

# create second struct for averaging of physical properties
phys_props_avg = PhysProps(pia)

# create struct for time-averaged output netCDF
ds_avg = NCDataHolder("couette_example.nc",
                        species_data, phys_props)

# init collision factors
collision_factors = create_collision_factors_array(pia, interaction_data, species_data, T_wall, Fnum)

# write out grid data
write_grid("scratch/data/couette_$(L)_$(nx)_grid.nc", grid)

# number of timesteps we are averaging for
n_avg = n_timesteps - avg_start + 1

for t in 1:n_timesteps

    # output timestep every 1000 timesteps
    if t % 1000 == 0
        println(t)
    end

    # collide particles
    for cell in 1:grid.n_cells
        ntc!(rng, collision_factors[1, 1, cell],
                collision_data, interaction_data, particles[1], pia, cell, 1, Δt, grid.cells[cell].V)
    end

    # convect particles
    convect_particles!(rng, grid, boundaries, particles[1], pia, 1, species_data, Δt)

    # sort particles
    sort_particles!(gridsorter, grid, particles[1], pia, 1)

    # compute props and do averaging
    if (t >= avg_start)
        compute_props_sorted!(particles, pia, species_data, phys_props)
        avg_props!(phys_props_avg, phys_props, n_avg)
    end
end

write_netcdf_phys_props(ds_avg, phys_props_avg, n_timesteps)

close_netcdf(ds_avg)
```