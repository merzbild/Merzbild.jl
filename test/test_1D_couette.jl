
@testset "couette test" begin

    seed = 1234
    T_wall = 300.0
    v_wall = 500.0
    L = 5e-4
    ndens = 5e22
    nx = 50
    ppc = 1000
    Δt = 2.59e-9
    output_freq = 1000
    n_timesteps = 6000

    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    # load particle and interaction data
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")
    interaction_data_path = joinpath(@__DIR__, "..", "data", "vhs.toml")
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

    sample_particles_equal_weight!(rng, grid, particles[1], pia, 1,
                                   species_data, ndens, T_wall, Fnum)

    # create collision structs
    collision_factors = create_collision_factors_array(1, grid.n_cells)
    collision_data = CollisionData()
    
    # create struct for computation of physical properties
    phys_props = PhysProps(pia)


    # create struct for netCDF output
    sol_path = joinpath(@__DIR__, "data", "tmp_couette.nc")
    ds = NCDataHolder(sol_path, species_data, phys_props)

    # # create struct for second netCDF, this one is for time-averaged 
    # ds_avg = NCDataHolder("scratch/data/avg_couette_$(L)_$(nx)_$(v_wall)_$(T_wall)_$(ppc)_after$(avg_start).nc",
    #                                   species_data, phys_props)
    # create second struct for averaging of physical properties
    # phys_props_avg = PhysProps(pia)

    # init collision structs
    for cell in 1:grid.n_cells
        collision_factors[1, 1, cell].sigma_g_w_max = estimate_sigma_g_w_max(interaction_data[1,1],
                                                                             species_data[1], T_wall, Fnum)
    end

    # compute and write data at t=0
    compute_props!(particles, pia, species_data, phys_props)
    write_netcdf_phys_props(ds, phys_props, 0)

    # n_avg = n_timesteps - avg_start + 1

    for t in 1:n_timesteps
        # collide particles
        for cell in 1:grid.n_cells
            ntc!(rng, collision_factors[1, 1, cell],
                 collision_data, interaction_data, particles[1], pia, cell, 1, Δt, grid.cells[cell].V)
        end

        # convect particles
        convect_particles!(rng, grid, boundaries, particles[1], pia, 1, species_data, Δt)

        # sort particles
        sort_particles!(gridsorter, grid, particles[1], pia, 1)

        # # compute props and do I/O
        # if (t < avg_start)
        #     if (t % output_freq == 0)
        #         compute_props_sorted!(particles, pia, species_data, phys_props)
        #     end
        # else
        #     compute_props_sorted!(particles, pia, species_data, phys_props)
        #     avg_props!(phys_props_avg, phys_props, n_avg)
        # end

        # if (t < avg_start)
        if (t % output_freq == 0)
            compute_props_sorted!(particles, pia, species_data, phys_props)
        end

        if (t % output_freq == 0)
            write_netcdf_phys_props(ds, phys_props, t)
        end
    end

    # write_netcdf_phys_props(ds_avg, phys_props_avg, n_timesteps)

    close_netcdf(ds)
    # close_netcdf(ds_avg)


    ref_sol_path = joinpath(@__DIR__, "data", "couette_0.0005_50_500.0_300.0_1000.nc")
    ref_sol = NCDataset(ref_sol_path, "r")
    sol = NCDataset(sol_path, "r")
 
    @test maximum(ref_sol["ndens"][:, 1, 1:5] .- sol["ndens"][:, 1, 1:5]) < 2 * eps()
    @test maximum(ref_sol["T"][:, 1, 1:5] .- sol["T"][:, 1, 1:5]) < 1.2e-13

    close(sol)
    rm(sol_path)

    # the reference solution was also time averaged betwen t=14001 and t=50000 and compared to a SPARTA solution
    # SPARTA version from 7 March 2024
    # see test/data/external
end
