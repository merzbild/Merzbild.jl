@testset "test couette VW serial with chunking and particle exchange + dynamic rebalancing" begin

    # this tests that when we run a VW couette sim with merging and chunking (even though it's serial)
    # and also does rebalancing based on number of collisions
    # we don't lose/overwrite particles, etc.
    n_chunks = 4

    T_wall = 300.0
    v_wall = 500.0
    L = 5e-4

    # set up a linear profile for the density so that number of collisions varies
    # and rebalacing occrus
    ndens_1 = 5e21
    ndens_2 = 5e22

    nx = 50
    ppc_sampled = 500
    merge_threshold = 180
    merge_target = 150

    Δt = 2.59e-9
    preallocation_margin_multiplier = 2.0

    n_timesteps = 50

    seed = 1234
    rng_chunks = [StableRNG(seed + i) for i in 0:n_chunks-1]

    lbq = LoadBalancerCellQ(nx, n_chunks)

    # load particle and interaction data
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")
    interaction_data_path = joinpath(@__DIR__, "..", "data", "vhs.toml")
    interaction_data::Array{Interaction, 2} = load_interaction_data(interaction_data_path, species_data)

    # create our grid and BCs
    grid = Grid1DUniform(L, nx)
    bc_list = (MaxwellWallBC1D(1, species_data, T_wall, [0.0, -v_wall, 0.0], 1.0),
               MaxwellWallBC1D(1, species_data, T_wall, [0.0, v_wall, 0.0], 1.0))

    # split cell indices into chunks
    cell_indices = Vector(1:nx)
    cell_chunks::Vector{UnitRange{Int64}} = [(chunk_index) for chunk_index in index_chunks(cell_indices; n=n_chunks)]

    # alternatively,
    # cell_chunks = [UnitRange{Int64}(chunk_index) for chunk_index in index_chunks(cell_indices; n=n_chunks)]

    # init per-chunk particle vectors, particle indexers, grid particle sorters
    n_particles_chunks = [floor(Int64, ppc_sampled * length(cell_chunk) * preallocation_margin_multiplier) for cell_chunk in cell_chunks]
    particles_chunks = [[ParticleVector(n_particles)] for n_particles in n_particles_chunks]
    pia_chunks = [ParticleIndexerArray(grid.n_cells, 1) for cell_chunk in cell_chunks]
    gridsorter_chunks = [GridSortInPlace(grid, n_particles) for n_particles in n_particles_chunks]

    # this is used for moving particles between chunks after they have been sorted into grid cells
    chunk_exchanger = ChunkExchanger(cell_chunks, nx)

    # sample particles
    # Fnum * ppc_sampled = Np in cell = ndens * V_cell
    Fnum = grid.cells[1].V * 0.5 * (ndens_1 + ndens_2) / ppc_sampled

    # sample particles per-chunk
    for (chunk_id, cell_chunk) in enumerate(cell_chunks)

        # sample explicitly varying ndens
        for cell in cell_chunk
            ndens = ndens_1 + (ndens_2 - ndens_1) * cell / nx
            @inbounds n_in_cell = ndens * grid.cells[cell].V

            ppc = n_in_cell / Fnum
            ppc_int = floor(Int64, ppc)
            remainder = ppc - ppc_int

            R = rand(rng_chunks[chunk_id])
            if R < remainder
                ppc_int += 1
            end

            @inbounds sample_particles_equal_weight!(rng_chunks[chunk_id], particles_chunks[chunk_id][1],
                                        pia_chunks[chunk_id], cell, 1,
                                        ppc_int, species_data[1].mass, T_wall, Fnum,
                                        grid.cells[cell].xlo, grid.cells[cell].xhi,
                                        0.0, 1.0,
                                        0.0, 1.0;
                                        distribution=:Maxwellian, vx0=0.0, vy0=0.0, vz0=0.0)
        end
    end
     
    # create collision structs
    collision_data = [CollisionData() for cell_chunk in cell_chunks]
    
    # create struct for computation of physical properties, sizes of pia are the same
    phys_props = PhysProps(pia_chunks[1])

    # create and estimate collision factors
    collision_factors = [create_collision_factors_array(pia, interaction_data, species_data, T_wall, Fnum * ppc_sampled / merge_target)
                         for pia in pia_chunks]

    # create merging structs
    oc_chunks = [OctreeMerge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel, max_Nbins=6000) for cell_chunks in cell_chunks]

    # merge and compute data at t=0
    for chunk_id in 1:n_chunks
        for cell in cell_chunks[chunk_id]
            merge_octree!(rng_chunks[chunk_id], oc_chunks[chunk_id], particles_chunks[chunk_id][1], pia_chunks[chunk_id], cell, 1, merge_target, grid)
        end
        squash_pia!(particles_chunks[chunk_id], pia_chunks[chunk_id])
        compute_props_sorted!(particles_chunks[chunk_id], pia_chunks[chunk_id], species_data, phys_props, cell_chunks[chunk_id])
    end

    for cell in 2:nx
        @test phys_props.n[cell] > phys_props.n[cell-1]
    end

    ndens_t0 = sum(phys_props.n)
    old_chunks = copy(cell_chunks)

    for t in 1:n_timesteps
        # collide, convect, sort particles
        for chunk_id in 1:n_chunks
            for cell in cell_chunks[chunk_id]
                ntc!(rng_chunks[chunk_id], collision_factors[chunk_id][1, 1, cell],
                               collision_data[chunk_id], interaction_data, particles_chunks[chunk_id][1],
                               pia_chunks[chunk_id], cell, 1, Δt, grid.cells[cell].V)
                update_lb_cellq!(lbq, chunk_id, cell, collision_factors[chunk_id][1, 1, cell].n_coll_performed; averaging_window=1.0)
                if pia_chunks[chunk_id].indexer[cell,1].n_local > merge_threshold
                    merge_octree!(rng_chunks[chunk_id], oc_chunks[chunk_id], particles_chunks[chunk_id][1], pia_chunks[chunk_id], cell, 1, merge_target, grid)
                    squash_pia!(particles_chunks[chunk_id], pia_chunks[chunk_id])
                end
            end

            convect_particles!(rng_chunks[chunk_id], grid, bc_list,
                                particles_chunks[chunk_id][1], pia_chunks[chunk_id],
                                1, species_data, Δt)
        

            # sort particles
            sort_particles!(gridsorter_chunks[chunk_id], grid, particles_chunks[chunk_id][1], pia_chunks[chunk_id], 1)
            update_occupancy_bounds!(chunk_exchanger, gridsorter_chunks[chunk_id], pia_chunks[chunk_id], chunk_id, 1)
        end

        # perform re-balancing and test that we find all cells, chunks are consistent
        if t%10 == 0
            rebalance_lb!(lbq)
            reset_lb!(lbq)
            cell_chunks = lbq.chunked_indices
            for i in 1:n_chunks-1
                @test cell_chunks[i][end] + 1 == cell_chunks[i+1][1]
            end
            for cell in 1:nx
                found = false
                for i in 1:n_chunks
                    if cell in cell_chunks[i]
                        found = true
                    end
                end
                @test found == true
            end
        end

        if t == 10
            # test that first chunk now has more cells than before
            # and last chunk has fewer cells than before
            @test length(cell_chunks[1]) > length(old_chunks[1])
            @test length(cell_chunks[end]) < length(old_chunks[end])
        end

        # move particles between chunks
        exchange_particles!(chunk_exchanger, particles_chunks, pia_chunks, cell_chunks, 1)

        # reset indexing, compute physical properties if needed
        for chunk_id in 1:n_chunks
            sort_particles_after_exchange!(chunk_exchanger, gridsorter_chunks[chunk_id],
                                           particles_chunks[chunk_id][1], pia_chunks[chunk_id],
                                           cell_chunks[chunk_id], 1)
            compute_props_sorted!(particles_chunks[chunk_id], pia_chunks[chunk_id],
                                  species_data, phys_props, cell_chunks[chunk_id])
        end
        
        # check indexing correctness
        pia_correct = [check_pia_is_correct(pia_chunks[chunk_id], 1) for chunk_id in 1:n_chunks]
        for chunk_id in 1:n_chunks
            @test pia_correct[chunk_id] == (1,0)
        end

        index_correct = [check_unique_index(particles_chunks[chunk_id][1], pia_chunks[chunk_id], 1) for chunk_id in 1:n_chunks]
        for chunk_id in 1:n_chunks
            @test index_correct[chunk_id] == (1,0)
        end

        # check that total number density is not lost
        @test abs(sum(phys_props.n) - ndens_t0) / ndens_t0< 4*eps() 
    end
end