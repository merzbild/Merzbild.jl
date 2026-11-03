# Inflow and outflow BCs

## Outflow BC
To allow gas particles to leave the domain, one can use an [`OutflowBC`](@ref)
boundary condition. Any particle touching a surface with such a BC gets deleted,
i.e. this corresponds to a vacuum boundary condition (equivalently, this can be used
as fully adsorbing BC), which is a good approximation of supersonic outflow at high Mach numbers.
The convection routines return the number of deleted particles.

An associated utility type [`Merzbild.ParticleDeletingBC`](@ref) is present, which is used
in conjunction with the [`Merzbild.has_particle_deleting_bc`](@ref) function to
determine at compile-time whether any of the boundaries of the domain have a particle-deleting BC
(since potentially other BCs may also delete particles, such as as the inflow BCs described below).

## Inflow BC
Inflow BCs are described by the [`InflowBC1D`](@ref) type.
They allow for sampling from a Maxwellian distribution with an optional drift/streaming
velocity pointing **into** the domain (a streaming velocity pointing **out** of the domain is
currently not supported).
An [`InflowBC1D`](@ref) also deletes any particles that hit it, so it is a particle-deleting BC.
The sampling is performed according to direct sampling in case of a zero streaming velocity (where
this is understood to mean the velocity multiplied with the inward surface normal of the boundary),
or via "Envelope 4" in case of a non-zero streaming velocity, as described in
[A.L. Garcia, W. Wagner, Generation of the Maxwellian inflow distribution, J. Comput. Phys., 2006](https://doi.org/10.1016/j.jcp.2006.01.025).

### Creating inflow BCs
To create an [`InflowBC1D`](@ref) that can
sample **uniformly weighted** particles at inflow (with particles
sampled from a Poisson process), one needs to specify the inflow number density, temperature, streaming velocity, and surface normal;
the sampling routine is then chosen automatically based on the streaming velocity.

For example, for inflow with a number
density `n_inflow`, temperature `T_inflow`, and x-velocity `vx_inflow`
from the left boundary (with a in-domain normal of `(1.0, 0.0, 0.0)`, therefore `vx_inflow` should be >= 0), one can write
```julia
# 1 is the species' index
# vx_inflow should be >= 0
inflow_bc = InflowBC1D(1, species_data, T_inflow, [vx_inflow, 0.0, 0.0], [1.0, 0.0, 0.0], n_inflow)
```
**Important note 1**: a specific instance of an `InflowBC1D` is species-specific, i.e. each species needs its own `bc_list` for emission.

**Important note 2**: the sampling method/envelope is chosen at run time based on the value of the streaming velocity
and the direction of the surface normal,
and is a type parameter of `InflowBC1D` (i.e. one has `InflowBC1D{S<:AbstractInflowSampler}`).
As such, in order to avoid on-the-fly memory allocations and possible type instabilities,
one can either use a constructor for `InflowBC1D` with an explicitly specified sampler,
or pass the instance (or the `Tuple` of boundary conditions containing it) through a function barrier.
So one has to create a list of InflowBCs based on the prescribed velocities, so that the sampling
algorithms are chosen, and then pass this list to another function that will actually run the simulation.
The latter is shown in the example below.

**Important note 3**: no run-time checks are performed
to see if the inflow BCs are used on the correct boundaries
(i.e. the first inflow boundary condition in the list of BCs
should correspond to the left wall), so incorrect order
will lead to erroneous results.

### Sampling particles from inflow BCs
The function that emits and immediately performs convection of the emitted particles is
[`emit_particles!`](@ref). It should be called **after** the convection of the existing particles (otherwise
the emitted+convected particles will be convected again),
and **before** particle sorting. The emitted particles can also update the computed properties of the surface they were
emitted from, and are counted as reflected particles for that purpose (with no corresponding incident particle).

## Example: outflow into vacuum.
This example describes a gas flowing into a 1-D domain from the left side, and flowing out into vacuum on the right side.
The full example can be found in `simulations/1D/inflow_into_vacuum_outflow.jl`.
Note that the function `run` only sets up the boundary conditions (for type stability), which are then passed to
`run_with_inflow_bc`, which performs the actual simulation.

```julia
using Merzbild
using Random
using TimerOutputs

"""
    run(seed, T_R, n_R, L, n_cells, n_emit_per_timestep,
        Δt, output_freq, n_timesteps, avg_start,
        name; vx_inflow=0.0, reorder_freq=10)

Run a simulation of a gas flowing from an equilibrium reservoir into vacuum. The domain is initially empty,
particles enter it from the left boundary (an `InflowBC1D`) and leave it from the right boundary (an `OutflowBC`),
as well as from the left boundary (`InflowBC1D` acts as supersonic outflow
for incoming particles). For a reservoir at rest, this corresponds to evaporation from a surface
into vacuum, with particles coming back to the evaporating surface condensing on it.
Fixed-weight particles are used, with the computational weight computed from the number flux of the reservoir
and the mean number of particles to be emitted per timestep.

Positional arguments:
* `seed`: random seed value
* `T_R`: temperature of the reservoir
* `n_R`: number density of the reservoir
* `L`: length of domain in meters
* `n_cells`: number of grid cells
* `n_emit_per_timestep`: mean number of particles emitted across the inflow boundary per timestep
* `Δt`: timestep size
* `output_freq`: how often instantaneous properties are written to disk before averaging starts
* `n_timesteps`: total number of timesteps to run for
* `avg_start`: timestep after which averaging begins
* `name`: name to append to start of output files

Keyword arguments:
* `vx_inflow`: x-velocity of the reservoir, cannot be negative (default is 0.0, i.e. the reservoir is at rest)
* `reorder_freq`: how frequently the ordering of particle indices is restored
"""
function run(seed, T_R, n_R, L, n_cells, n_emit_per_timestep,
             Δt, output_freq, n_timesteps, avg_start,
             name; vx_inflow=0.0, reorder_freq=10)
    particles_data_path = joinpath(MERZBILD_DATA_PATH, "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")

    # the type of the inflow BC depends on the value of vx_inflow and is not known at compile time,
    # so the simulation is run in a separate function (a function barrier) to avoid dynamic dispatch
    # in the time loop
    inflow_bc = InflowBC1D(1, species_data, T_R, [vx_inflow, 0.0, 0.0], [1.0, 0.0, 0.0], n_R)
    run_with_inflow_bc(seed, species_data, inflow_bc, T_R, n_R, L, n_cells, n_emit_per_timestep,
                       Δt, output_freq, n_timesteps, avg_start, name, vx_inflow, reorder_freq)
end

function run_with_inflow_bc(seed, species_data, inflow_bc, T_R, n_R, L, n_cells, n_emit_per_timestep,
                            Δt, output_freq, n_timesteps, avg_start, name, vx_inflow, reorder_freq)
    reset_timer!()

    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    # load interaction data
    interaction_data_path = joinpath(MERZBILD_DATA_PATH, "vhs.toml")
    interaction_data::Array{Interaction, 2} = load_interaction_data(interaction_data_path, species_data)

    # create our grid and BCs
    grid = Grid1DUniform(L, n_cells)
    bc_list = (inflow_bc, OutflowBC())

    # the boundaries of a 1-D grid have unit area, so Fnum * n_emit_per_timestep = ndens_flux * Δt
    ndens_flux = inflow_bc.ndens_flux
    Fnum = ndens_flux * Δt / n_emit_per_timestep

    λ = mean_free_path(interaction_data, 1, n_R, T_R)
    cfreq = mean_collision_frequency(interaction_data, 1, species_data, n_R, T_R)
    v_thermal = sqrt(2 * k_B * T_R / species_data[1].mass)

    println("Reservoir properties:")
    println("Γ (flux) = $ndens_flux; Fnum = $Fnum")
    println("λ = $λ; Kn(reservoir) = $(λ / L); λ/dx = $(λ / grid.Δx)")
    println("1/ν_coll = $(1/cfreq); Δt * ν_coll = $(Δt * cfreq)")
    println("CFL time: $(grid.Δx / v_thermal); L / v_thermal in timesteps: $(L / v_thermal / Δt)")

    # init particle vector, particle indexer, grid particle sorter
    # the domain is initially empty, the particle vector and sorter are resized as needed
    n_particles = ceil(Int64, 2 * n_R * L / Fnum)
    particles = [ParticleVector{1}(n_particles)]
    pia = ParticleIndexerArray(grid.n_cells, 1)
    gridsorter = GridSortInPlace(grid, n_particles)

    # create collision structs
    collision_data = CollisionData()
    collision_factors = create_collision_factors_array(pia, interaction_data, species_data, T_R, Fnum)

    # create structs for computation and averaging of physical and surface properties
    phys_props = PhysProps(pia; ndens_not_Np=true)
    phys_props_avg = PhysProps(pia; ndens_not_Np=true)
    surf_props = SurfProps(pia, grid)
    surf_props_avg = SurfProps(pia, grid)

    # create structs for netCDF output
    fname_base = "$(name)_$(L)_$(n_cells)_$(T_R)_$(n_R)_$(vx_inflow)"
    ds = NCDataHolder("scratch/data/$(fname_base).nc", species_data, phys_props)
    ds_avg = NCDataHolder("scratch/data/avg_$(fname_base)_after$(avg_start).nc", species_data, phys_props)
    ds_surf_avg = NCDataHolderSurf("scratch/data/avg_$(fname_base)_surf_after$(avg_start).nc",
                                   species_data, surf_props_avg)

    write_grid("scratch/data/$(name)_$(L)_$(n_cells)_grid.nc", grid)

    n_avg = n_timesteps - avg_start + 1

    index_inv_map = zeros(Int64, n_particles)

    for t in 1:n_timesteps
        if t % 1000 == 0
            println("$t, n_particles = $(pia.n_total[1])")
        end

        for cell in 1:grid.n_cells
            @timeit "collide" ntc_equal_weight!(rng, collision_factors[1, 1, cell],
                                                collision_data, interaction_data, particles[1],
                                                pia, cell, 1, Δt, grid.cells[cell].V)
        end

        # convect particles, then emit new ones; these are convected for the part of the timestep
        # remaining after they enter the domain
        if (t < avg_start)
            @timeit "convect" convect_particles_and_compute_cell!(rng, grid, bc_list, particles[1], pia, 1,
                                                                  species_data, Δt)
            @timeit "emit" emit_particles!(rng, grid, bc_list, particles[1], pia, 1, species_data, Fnum, Δt)
        else
            @timeit "convect + surface compute" convect_particles_and_compute_cell!(rng, grid, bc_list, particles[1], pia, 1,
                                                                                    species_data, surf_props, Δt)
            @timeit "emit + surface compute" emit_particles!(rng, grid, bc_list, particles[1], pia, 1, species_data,
                                                             surf_props, Fnum, Δt)
            @timeit "avg surfprops" avg_props!(surf_props_avg, surf_props, n_avg)
        end

        # sort particles
        @timeit "sort" sort_particles!(gridsorter, particles[1], pia, 1)

        if t % reorder_freq == 0
            @timeit "restore ordering" restore_particle_ordering!(particles[1], pia, 1, index_inv_map)
        end

        # compute props and do I/O
        if (t < avg_start)
            if (t % output_freq == 0)
                @timeit "props compute" compute_props_sorted!(particles, pia, species_data, phys_props, grid)
                @timeit "I/O" write_netcdf(ds, phys_props, t)
            end
        else
            @timeit "props compute" compute_props_sorted!(particles, pia, species_data, phys_props, grid)
            @timeit "avg physprops" avg_props!(phys_props_avg, phys_props, n_avg)
        end
    end

    @timeit "I/O final" write_netcdf(ds_avg, phys_props_avg, n_timesteps)
    @timeit "I/O final" write_netcdf(ds_surf_avg, surf_props_avg, n_timesteps)

    close_netcdf(ds)
    close_netcdf(ds_avg)
    close_netcdf(ds_surf_avg)

    flux_emitted = -surf_props_avg.flux_reflected[1, 1]
    flux_returned = surf_props_avg.flux_incident[1, 1]
    if avg_start > n_timesteps
        println("Averaging starts at timestep $avg_start > $n_timesteps timesteps the simulation has; ", 
                "so no statistics can be computed")
    else
        println("Mass fluxes: emitted = $flux_emitted (expected $(ndens_flux * species_data[1].mass)), ",
                "returned to inflow boundary = $flux_returned, outflow = $(surf_props_avg.flux_incident[2, 1])")
        println("Fraction of emitted mass flux returning to inflow boundary: $(flux_returned / flux_emitted)")
    end

    print_timer()
end

# run 5000 timesteps (with dt = 5e-7 s) with
# following inflow parameters: T=300.0, ndens of 1e21, vx=0 (by default),
# 100 particles emitted per timestep
# domain length of length 0.05 m discretized with 100 cells,
# writing output every 500 timesteps and starting averaging after 2000 timesteps
const n_t = 5000
run(1234, 300.0, 1e21, 0.05, 100, 100, 5e-7, 500, n_t, 2000, "inflow_into_vacuum")
```