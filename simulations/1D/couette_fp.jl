include("../../src/Merzbild.jl")

using ..Merzbild
using Random
using TimerOutputs

function run(seed, T_wall, v_wall, L, ndens, nx, ppc, Δt, output_freq, n_timesteps, avg_start)
    reset_timer!()

    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    # load particle and interaction data
    particles_data_path = joinpath("data", "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")
    interaction_data_path = joinpath("data", "vhs.toml")
    interaction_data::Array{Interaction, 2} = load_interaction_data(interaction_data_path, species_data)

    # create our grid and BCs
    grid = Grid1DUniform(L, nx)
    boundaries = MaxwellWalls(species_data, T_wall, T_wall, -v_wall, v_wall, 1.0, 1.0)

    # init particle vector, particle indexer, grid particle sorter
    n_particles = ppc * nx
    particles = particles = [ParticleVector(n_particles)]
    pia = ParticleIndexerArray(grid.n_cells, 1)
    gridsorter = GridSortInPlace(grid, n_particles)

    # sample particles
    # Fnum * ppc = Np in cell = ndens * V_cell
    Fnum = grid.cells[1].V * ndens / ppc
    @timeit "sampling" sample_particles_equal_weight!(rng, grid, particles[1], pia, 1,
                                                      species_data, ndens, T_wall, Fnum)

    # create collision structs
    collision_factors = create_collision_factors_array(1, grid.n_cells)
    collision_data_fp = CollisionDataFP()
    
    # create struct for computation of physical properties
    phys_props = PhysProps(pia)

    # create second struct for averaging of physical properties
    phys_props_avg = PhysProps(pia)

    # create struct for netCDF output
    ds = NCDataHolder("couette_$(L)_$(nx)_$(v_wall)_$(T_wall)_$(ppc).nc", species_data, phys_props)

    # create struct for second netCDF, this one is for time-averaged 
    ds_avg = NCDataHolder("avg_couette_$(L)_$(nx)_$(v_wall)_$(T_wall)_$(ppc)_after$(avg_start).nc",
                          species_data, phys_props)

    # init collision structs
    for cell in 1:grid.n_cells
        collision_factors[1, 1, cell].sigma_g_w_max = estimate_sigma_g_w_max(interaction_data[1,1],
                                                                             species_data[1], T_wall, Fnum)
    end

    # compute and write data at t=0
    compute_props!(particles, pia, species_data, phys_props)
    write_netcdf_phys_props(ds, phys_props, 0)
    write_grid("couette_$(L)_$(nx)_grid.nc", grid)

    n_avg = n_timesteps - avg_start + 1

    for t in 1:n_timesteps
        if t % 1000 == 0
            println(t)
        end

        # collide particles
        for cell in 1:grid.n_cells
            @timeit "collide" fp!(rng, collision_data_fp, interaction_data[1, 1], species_data, particles[1], pia, cell, 1, Δt, grid.cells[cell].V)
        end

        # convect particles
        @timeit "convect" convect_particles!(rng, grid, boundaries, particles[1], pia, 1, species_data, Δt)

        # sort particles
        @timeit "sort" sort_particles!(gridsorter, grid, particles[1], pia, 1)

        # compute props and do I/O

        if (t < avg_start)
            if (t % output_freq == 0)
                @timeit "props compute" compute_props_sorted!(particles, pia, species_data, phys_props)
            end
        else
            @timeit "props compute" compute_props_sorted!(particles, pia, species_data, phys_props)
            avg_props!(phys_props_avg, phys_props, n_avg)
        end


        if (t % output_freq == 0)
            @timeit "I/O" write_netcdf_phys_props(ds, phys_props, t)
        end
    end

    @timeit "I/O" write_netcdf_phys_props(ds_avg, phys_props_avg, n_timesteps)

    close_netcdf(ds)
    close_netcdf(ds_avg)

    print_timer()
end

N_steps = 50000

run(1234, 300.0, 500.0, 5e-4, 5e22, 50, 1000, 2.59e-9, 1000, N_steps, 14000)