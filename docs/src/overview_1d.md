# 1D DSMC simulations

In this section, setting up fixed-weight DSMC simulations on a uniform 1D grid will be discussed,
along with computation of surface properties due to particle-surface interactions. Some of the
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

We then initialize an `Vector` of `ParticleVector`'s to store our particles for each species. We now use 1-dimensional position vectors:
```julia
particles = [ParticleVector{1}(n_particles)]
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

Since sorting indices only can lead to increase fragmentation of the particle layout in memory,
a function [`count_disordered_particles`](@ref) is available that counts the number of 
non-continuously laid out particles; this can serve as a metric as to whether the underlying
particles (and not just their indices) need to be re-sorted.
The fragmentation of particles indices can be fixed by calling [`restore_particle_ordering!`](@ref) after the sorting routine;
calling it every 10 timesteps or so gives a good balance between cost of re-indexing and speed-up due to improved memory access.

The re-indexing can be called as follows:
```julia
restore_particle_ordering!(particles[species_id], index_inv_map)
```

Here `index_inv_map` is a pre-allocated array of integers of length equal to the number of particles in the simulation (if it is smaller,
it will be resized in by the `restore_particle_ordering!` function).

## Creating boundary conditions
Next, we need to create boundary conditions for the left and right walls.
Currently, fully diffusely reflecting walls, specularly reflecting walls, and a Maxwell
boundary with a user-defined accommodation coefficient are available, with specific version optimized for 1D
simulations. Due to dynamical dispatch, more generic conditions can be used in place of the specialized 1D ones
even in 1D simulations, but this will be less efficient.

We instantiate two fully diffuse boundary conditions and pack them into a `Tuple` (the order being left and right wall,
`1` being the species index):

```julia
bc_list = (FullyDiffuseBC1D(1, species_data, T_wall, [0.0, -v_wall, 0.0]),
           FullyDiffuseBC1D(1, species_data, T_wall, [0.0, v_wall, 0.0]))
```

## Calculation of surface properties
To compute surface properties due to particle-surface interactions, one needs to first set up the corresponding struct
that will hold the computed values. This is done by.
```julia
surf_props = SurfProps(pia, grid)
```
If one wants to compute the properties on a given timestep, the structure needs to be passed to the convection
routine --- otherwise they will not be computed, as one needs to know the particle properties before and after
its interaction with a surface.

## I/O of surface properties
To set up NetCDF output of computed surface properties, one has code similar to the one used for the output of
grid quantities:
```julia
ds_surf = NCDataHolderSurf("scratch/data/couette_example_surf.nc", species_data, surf_props)

for t in 1:n_timesteps
# simulation loop here
    write_netcdf(ds_surf, surf_props, t)  # write computed surface properties to file
end
```

## Performing convection
Having set up the grid and boundary conditions, we can convect particles.
This is done by calling the [`convect_particles!`](@ref) function.
The convection should be followed by particle sorting before any computations of physical properties are done.

```julia
convect_particles!(rng, grid, bc_list, particles[species_id], pia, species_id, species_data, Δt)
```

The function `convect_particles` as called above will **not** compute surface properties. To do that,
a `SurfProps` instance needs to be passed:
```julia
convect_particles!(rng, grid, bc_list, particles[species_id], pia, species_id, species_data, surf_props, Δt)
```

## Convection and sorting with precomputed particle/cell indices
The approach described above assumes that during the convection process, the indices of the cells the particles find themselves in
are not computed; therefore the call to `sort_particles!` requires passing in the `grid` instance, and the `sort_particles!`
calls a `get_cell` function internally. For 1-D uniform grids, this is an efficient operation, but for other grid types,
the computation of the cell index given only the particle position can be more expensive than keeping track of the cell index inside the convection routine.

To this purpose, one can call [`convect_particles_and_compute_cell!`](@ref), which will also set the values of the `cell` array of the `ParticleVector`
instance. One can then call another version of `sort_particles!` that does **not** take the grid as a parameter and instead uses the pre-computed `cell` values
to sort the particles:

```
sort_particles!(gridsorter, particles[species_id], pia, species_id)
```

## Bringing it all together
Now we can combine all the pieces to set up a simulation of a single-species Couette flow in a channel
with a width of 0.5 mm, discretized with 50 cells. The y-velocity of the left wall is assumed to be -500 m/s,
and that of the right wall 500 m/s; the temperature of both walls is 300 K. The solution is initialized with 100
particles per cell and a number density of 5e22 1/m^3. A timestep of 2.59 ns is used.
The simulation runs for 50K steps and the solution is time-averaged after the first 14K steps.

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
particles_data_path = joinpath(MERZBILD_DATA_PATH, "particles.toml")
species_data = load_species_data(particles_data_path, "Ar")
interaction_data_path = joinpath(MERZBILD_DATA_PATH, "vhs.toml")
interaction_data::Array{Interaction, 2} = load_interaction_data(interaction_data_path, species_data)

# create our grid and BCs
grid = Grid1DUniform(L, nx)
bc_list = (FullyDiffuseBC1D(1, species_data, T_wall, [0.0, -v_wall, 0.0]),
           FullyDiffuseBC1D(1, species_data, T_wall, [0.0, v_wall, 0.0]))

# init particle vector, particle indexer, grid particle sorter
# we will not be creating or destroying any particles, so we can compute the exact number
# of particles we will have in the simulation
n_particles = ppc * nx
particles = [ParticleVector{1}(n_particles)]
pia = ParticleIndexerArray(grid.n_cells, 1)
gridsorter = GridSortInPlace(grid, n_particles)
index_inv_map = zeros(Int64, n_particles)

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

# create struct for computation of surface properties
surf_props = SurfProps(pia, grid)

# create second struct for averaging of surface properties
surf_props_avg = SurfProps(pia, grid)

# create struct for time-averaged output netCDF for grid properties
ds_avg = NCDataHolder("scratch/data/couette_example.nc",
                      species_data, phys_props)

# create struct for time-averaged output netCDF for surface properties
ds_surf_avg = NCDataHolderSurf("scratch/data/couette_example_surf.nc",
                               species_data, surf_props)

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
    convect_particles!(rng, grid, bc_list, particles[1], pia, 1, species_data, surf_props, Δt)

    # sort particles
    sort_particles!(gridsorter, grid, particles[1], pia, 1)

    # restore indexing/ordering every 10 timesteps
    if t%10 == 0
        restore_particle_ordering!(particles[1], index_inv_map)
    end

    # compute props and do averaging
    if (t >= avg_start)
        compute_props_sorted!(particles, pia, species_data, phys_props)
        avg_props!(phys_props_avg, phys_props, n_avg)
        avg_props!(surf_props_avg, surf_props, n_avg)
    end
end

write_netcdf(ds_avg, phys_props_avg, n_timesteps)
write_netcdf(ds_surf_avg, surf_props_avg, n_timesteps)

close_netcdf(ds_avg)
close_netcdf(ds_surf_avg)
```