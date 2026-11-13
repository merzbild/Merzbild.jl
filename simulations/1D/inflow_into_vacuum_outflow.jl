using Merzbild
using Random
using TimerOutputs

isdir("scratch/data") || error("Output directory `scratch/data` does not exist (the path is relative to the working directory). Create it, e.g. with `mkdir -p scratch/data`.")

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

Instantaneous output is written every `output_freq` timesteps before averaging starts to a file with name
`[name]_[L]_[n_cells]_[T_R]_[n_R]_[vx_inflow].nc`; time-averaged output (physical and surface properties)
follows a similar naming scheme. All files are written to `scratch/data/`.

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

const n_t = 5000
run(1234, 300.0, 1e21, 0.05, 100, 100, 5e-7, 500, n_t, 2000, "inflow_into_vacuum")
