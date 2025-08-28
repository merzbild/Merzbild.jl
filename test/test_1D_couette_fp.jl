
@testset "couette test linear Fokker-Planck" begin

    T_wall = 300.0
    v_wall = 500.0
    L = 5e-4
    ndens = 5e22
    nx = 50
    ppc = 200
    Δt = 2.59e-9
    output_freq = 1000
    n_timesteps = 3000

    seed = 1234
    rng = StableRNG(seed)

    # load particle and interaction data
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")
    interaction_data_path = joinpath(@__DIR__, "..", "data", "vhs.toml")
    interaction_data::Array{Interaction, 2} = load_interaction_data(interaction_data_path, species_data)

    # create our grid and BCs
    grid = Grid1DUniform(L, nx)
    boundaries = MaxwellWalls1D(species_data, T_wall, T_wall, -v_wall, v_wall, 1.0, 1.0)

    # init particle vector, particle indexer, grid particle sorter
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
    collision_data_fp = CollisionDataFP(ppc * 2)
    
    # create struct for computation of physical properties
    phys_props = PhysProps(pia)


    # create struct for netCDF output
    sol_path = joinpath(@__DIR__, "data", "tmp_couette.nc")
    ds = NCDataHolder(sol_path, species_data, phys_props)

    # compute and write data at t=0
    compute_props!(particles, pia, species_data, phys_props)
    write_netcdf(ds, phys_props, 0)


    for t in 1:n_timesteps
        # collide particles
        for cell in 1:grid.n_cells
            fp_linear!(rng, collision_data_fp, interaction_data[1, 1],
                       species_data, particles[1], pia, cell, 1, Δt, grid.cells[cell].V)
        end

        # convect particles
        convect_particles!(rng, grid, boundaries, particles[1], pia, 1, species_data, Δt)

        # sort particles
        sort_particles!(gridsorter, grid, particles[1], pia, 1)

        if (t % output_freq == 0)
            compute_props_sorted!(particles, pia, species_data, phys_props)
        end

        if (t % output_freq == 0)
            write_netcdf(ds, phys_props, t)
        end
    end


    close_netcdf(ds)


    ref_sol_path = joinpath(@__DIR__, "data", "couette_0.0005_50_500.0_300.0_100_fp_linear.nc")
    ref_sol = NCDataset(ref_sol_path, "r")
    sol = NCDataset(sol_path, "r")

    for nt in 1:4
        @test sum(sol["np"][:, 1, nt]) == 50 * ppc
    end
 
    @test maximum(ref_sol["np"][:, 1, 1] .- sol["np"][:, 1, 1]) < eps()
    @test maximum(ref_sol["ndens"][:, 1, 1:4] .- sol["ndens"][:, 1, 1:4]) < 2 * eps()
    @test maximum(ref_sol["T"][:, 1, 1:4] .- sol["T"][:, 1, 1:4]) < 2.4e-13

    close(sol)
    rm(sol_path)
end
