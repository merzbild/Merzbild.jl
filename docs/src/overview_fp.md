# Fokker-Planck simulations

Instead of resolving individual binary collisions, as done in DSMC (see [Fixed-weight DSMC simulations](@ref)),
the Fokker-Planck (FP) approach models the effect of many collisions as a continuous stochastic (drift-diffusion)
process acting on the particle velocities. Each particle velocity is updated by integrating a Langevin equation
over the timestep, rather than by pairing particles up and colliding them.

The main advantage is cost: the FP update is applied to each particle once per timestep, so its cost per cell scales
linearly with the number of particles in the cell, independent of the local collision rate. In near-continuum
(small mean free path) regions, where a DSMC simulation would need very many collisions --- and hence a very small
timestep and many particles --- the FP model can use a coarser timestep, making it much cheaper. It is thus well
suited to dense or near-continuum flows, and to hybrid schemes.

Currently, the linear Fokker-Planck model of
[Gorji, Torrilhon and Jenny (2011)](https://doi.org/10.1017/jfm.2011.188)
is implemented for single-species elastic collisions, for both fixed- and variable-weight particles. Since
no pair-wise particle collisions are performed, no particle merging is required even if variable-weight particles are used.

## Loading interaction data
As for DSMC, the FP model needs interaction data for the species being collided: the VHS parameters are used to
compute the viscosity, which sets the relaxation time ``\tau = 2\mu/p`` of the drift-diffusion process. This data is
loaded exactly as in the DSMC case, using [`load_interaction_data`](@ref), and stored in an `n_species x n_species`
matrix of `Interaction` instances (see [Fixed-weight DSMC simulations](@ref) for details).

## FP-specific data: CollisionDataFP
The Langevin update requires, per cell, a set of normally distributed random numbers (one per velocity component per
particle) plus some temporary quantities (the cell-average velocity, and the mean/standard deviation used to rescale
the random numbers so that momentum and energy are conserved exactly). These are stored in a [`CollisionDataFP`](@ref)
instance. Unlike the NTC algorithm, no per-cell collision factors are needed --- a single `CollisionDataFP` instance
is created and reused for all cells.

The internal arrays are resized automatically if a cell holds more particles than expected, but it is cheaper to
pre-allocate them for the maximum expected number of particles in a cell by passing that number to the constructor:
```julia
# pre-allocate for up to ppc*2 particles in a cell
collision_data_fp = CollisionDataFP(ppc * 2)
```

## Performing FP collisions
Single-species elastic collisions using the linear FP model are performed by calling [`fp_linear!`](@ref) once per
cell:
```julia
fp_linear!(rng, collision_data_fp, interaction_data[1, 1], particles[1], pia, cell, 1, species_data, Δt, grid.cells[cell].V)
```
Here `interaction_data[1, 1]` is the (single) self-interaction, `cell` is the index of the cell being collided,
`Δt` is the timestep, and the last argument is the volume of the physical cell. The function updates the particle
velocities in place; it computes the relaxation time from the local temperature and density
(see [`Merzbild.compute_relaxation_time`](@ref)), advances the velocities via the Langevin step, and finally rescales
them to conserve energy exactly.

Note that the FP model requires at least 7 particles in a cell to produce a meaningful update; if fewer are present,
[`fp_linear!`](@ref) returns without modifying the particles.

## Example: bringing it all together
The FP collision routine slots into the same grid/boundary-condition/convection machinery used for 1D DSMC
simulations (see [1D DSMC simulations](@ref)): the only difference is that the NTC collision step is replaced by a
call to [`fp_linear!`](@ref), and no collision factors need to be created. The example below simulates a single-species
Couette flow in a channel of width 0.5 mm, discretized with 50 cells. The y-velocity of the left wall is -500 m/s and
that of the right wall 500 m/s; both walls are at 300 K. The solution is initialized with 100 particles per cell and
a number density of 5e22 1/m^3, and runs for 50K steps with time-averaging after the first 14K steps.

```julia
using Merzbild
using Random

# set our random seed for reproducibility
seed = 1234
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
species_data = load_species_data(joinpath(MERZBILD_DATA_PATH, "particles.toml"), "Ar")
interaction_data = load_interaction_data(joinpath(MERZBILD_DATA_PATH, "vhs.toml"), species_data)

# create our grid and BCs
grid = Grid1DUniform(L, nx)
bc_list = (FullyDiffuseBC1D(species_data, 1, T_wall, [0.0, -v_wall, 0.0]),
           FullyDiffuseBC1D(species_data, 1, T_wall, [0.0, v_wall, 0.0]))

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

# create FP collision struct, pre-allocating for the expected number of particles per cell
collision_data_fp = CollisionDataFP(ppc * 2)

# create struct for computation of physical properties
phys_props = PhysProps(pia)

# create second struct for averaging of physical properties
phys_props_avg = PhysProps(pia)

# create struct for time-averaged output netCDF
ds_avg = NCDataHolder("scratch/data/couette_fp_example.nc", species_data, phys_props)

# write out grid data
write_grid("scratch/data/couette_fp_$(L)_$(nx)_grid.nc", grid)

# number of timesteps we are averaging for
n_avg = n_timesteps - avg_start + 1

for t in 1:n_timesteps

    # output timestep every 1000 timesteps
    if t % 1000 == 0
        println(t)
    end

    # collide particles: one FP update per cell, no collision factors needed
    for cell in 1:grid.n_cells
        fp_linear!(rng, collision_data_fp, interaction_data[1, 1], 
                   particles[1], pia, cell, 1, species_data, Δt, grid.cells[cell].V)
    end

    # convect particles
    convect_particles!(rng, grid, bc_list, particles[1], pia, 1, species_data, Δt)

    # sort particles
    sort_particles!(gridsorter, grid, particles[1], pia, 1)

    # restore indexing/ordering every 10 timesteps
    if t % 10 == 0
        restore_particle_ordering!(particles[1], index_inv_map)
    end

    # compute props and do averaging
    if (t >= avg_start)
        compute_props_sorted!(particles, pia, species_data, phys_props)
        avg_props!(phys_props_avg, phys_props, n_avg)
    end
end

write_netcdf(ds_avg, phys_props_avg, n_timesteps)
close_netcdf(ds_avg)
```

## Summary
Now we have an overview of how to

1. Set up the FP-specific `CollisionDataFP` data structure
2. Perform single-species elastic collisions using the linear Fokker-Planck model
3. Drop the FP collision step into the existing 1D grid/convection machinery in place of DSMC collisions
