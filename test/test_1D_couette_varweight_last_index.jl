
@testset "couette test variable weight + no intermediate squashing" begin

    T_wall = 300.0
    v_wall = 500.0
    L = 5e-4
    ndens = 5e22
    nx = 50
    ppc = 1000
    Δt = 2.59e-9
    output_freq = 100
    n_timesteps = 2000
    merge_threshold = 165
    merge_target = 150

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
    collision_data = CollisionData()
    
    # create struct for computation of physical properties
    phys_props = PhysProps(pia)


    # create struct for netCDF output
    sol_path = joinpath(@__DIR__, "data", "tmp_couette.nc")
    ds = NCDataHolder(sol_path, species_data, phys_props)

    # init collision factors
    collision_factors = create_collision_factors_array(pia, interaction_data, species_data, T_wall, Fnum)


    oc = OctreeN2Merge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel, max_Nbins=6000)

    for cell in 1:grid.n_cells
        if pia.indexer[cell,1].n_local > merge_threshold
            merge_octree_N2_based!(rng, oc, particles[1], pia, cell, 1, merge_target, grid)
        end
    end
    squash_pia!(particles, pia)

    # compute and write data at t=0
    compute_props!(particles, pia, species_data, phys_props)
    write_netcdf(ds, phys_props, 0)

    # n_avg = n_timesteps - avg_start + 1

    for t in 1:n_timesteps
        # collide particles
        merged = false

        for cell in 1:grid.n_cells
            ntc!(rng, collision_factors[1, 1, cell],
                 collision_data, interaction_data, particles[1], pia, cell, 1, Δt, grid.cells[cell].V)

            if pia.indexer[cell,1].n_local > merge_threshold
                merge_octree_N2_based!(rng, oc, particles[1], pia, cell, 1, merge_target, grid)
                merged = true
            end
        end

        if t % 10 == 0
            @test check_pia_is_correct(pia, 1) == (true, 0)
        end

        if merged
            squash_pia!(particles, pia)
        end
        # convect particles
        convect_particles!(rng, grid, boundaries, particles[1], pia, 1, species_data, Δt)

        # sort particles
        sort_particles!(gridsorter, grid, particles[1], pia, 1)

        # if (t < avg_start)
        if (t % output_freq == 0)
            compute_props_sorted!(particles, pia, species_data, phys_props)
        end

        if (t % output_freq == 0)
            write_netcdf(ds, phys_props, t)
        end
    end

    # write_netcdf(ds_avg, phys_props_avg, n_timesteps)

    close_netcdf(ds)
    # close_netcdf(ds_avg)

    @test check_pia_is_correct(pia, 1) == (true, 0)
    @test check_unique_index(particles[1], pia, 1) == (true, 0)
    @test check_unique_buffer(particles[1]) == (true, 0)

    ref_sol_path = joinpath(@__DIR__, "data", "couette_0.0005_50_500.0_300.0_1000_vw200to150.nc")
    ref_sol = NCDataset(ref_sol_path, "r")
    sol = NCDataset(sol_path, "r")
 
    ndens_conservation = true
    ref_ndens = 2.5e19

    for t in 1:5
        if abs(sum(sol["ndens"][:, 1, t]) - ref_ndens) / ref_ndens > 2*eps()
            ndens_conservation = false
        end
    end

    @test ndens_conservation == true

    close(sol)
    rm(sol_path)

    # the reference solution was also time averaged betwen t=14001 and t=50000 and compared to a SPARTA solution
    # SPARTA version from 7 March 2024
    # see test/data/external
end
