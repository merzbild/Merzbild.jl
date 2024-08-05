include("../../src/Merzbild.jl")

using ..Merzbild
using Random
using TimerOutputs

function run(seed, T_wall, v_wall, L, ndens, nx, ppc, Δt, output_freq, n_timesteps, avg_start)
    reset_timer!()

    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    # create our grid and BCs
    grid = create_grid1D_uniform(L, nx)
    boundaries = create_1D_boundaries(T_wall, T_wall, -v_wall, v_wall, 1.0, 1.0)

    # load particle and interaction data
    particles_data_path = joinpath("data", "particles.toml")
    species_data = load_species_list(particles_data_path, "Ar")
    interaction_data_path = joinpath("data", "vhs.toml")
    interaction_data::Array{Interaction, 2} = load_interaction_data(interaction_data_path, species_data)

    # init particle vector, particle indexer, grid particle sorter
    n_particles = ppc * nx
    particles = particles = [create_particle_vector(n_particles)]
    pia = create_particle_indexer_array(grid.n_cells, 1)
    gridsorter = create_grid_sort_inplace(grid, n_particles)

    # sample particles
    # Fnum * ppc = Np in cell = ndens * V_cell
    Fnum = grid.cells[1].V * ndens / ppc
    @timeit "sampling" sample_particles_equal_weight!(rng, grid, particles[1], pia, 1,
                                                      species_data, ndens, T_wall, Fnum)

    # create collision structs
    collision_factors = create_collision_factors(1, grid.n_cells)
    collision_data = create_collision_data()
    
    # create struct for computation of physical properties
    phys_props = create_props(pia)

    # create second struct for averaging of physical properties
    phys_props_avg = create_props(pia)

    # create struct for netCDF output
    ds = create_netcdf_phys_props("scratch/data/couette_$(L)_$(nx)_$(v_wall)_$(T_wall)_$(ppc).nc", species_data, phys_props)

    # create struct for second netCDF, this one is for time-averaged 
    ds_avg = create_netcdf_phys_props("scratch/data/avg_couette_$(L)_$(nx)_$(v_wall)_$(T_wall)_$(ppc)_after$(avg_start).nc",
                                      species_data, phys_props)

    # init collision structs
    for cell in 1:grid.n_cells
        collision_factors[1, 1, cell].sigma_g_w_max = estimate_sigma_g_w_max(interaction_data[1,1],
                                                                             species_data[1], T_wall, Fnum)
    end

    # compute and write data at t=0
    compute_props!(particles, pia, species_data, phys_props)
    write_netcdf_phys_props(ds, phys_props, 0)

    n_avg = n_timesteps - avg_start + 1

    for t in 1:n_timesteps
        if t % 1000 == 0
            println(t)
        end

        # collide particles
        for cell in 1:grid.n_cells
            @timeit "collide" ntc!(rng, collision_factors[1, 1, cell],
                                   collision_data, interaction_data, particles[1], pia, cell, 1, Δt, grid.cells[cell].V)
        end

        # convect particles
        @timeit "convect" convect_particles!(rng, grid, boundaries, particles[1], pia, 1, species_data, Δt)

        # sort particles
        @timeit "sort" sort_particles!(gridsorter, grid, particles[1], pia, 1)

        # compute props and do I/O

        if (t < avg_start)
            if (t % output_freq == 0)
                @timeit "props compute" compute_props_sorted_without_moments!(particles, pia, species_data, phys_props)
            end
        else
            @timeit "props compute" compute_props_sorted_without_moments!(particles, pia, species_data, phys_props)
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

run(1234, 300.0, 500.0, 5e-4, 5e22, 50, 1000, 2.59e-9, 1000, 50000, 14000)