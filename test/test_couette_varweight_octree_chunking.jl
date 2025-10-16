@testset "test couette VW serial with chunking and particle exchange" begin

    # this tests that when we run a VW couette sim with merging and chunking (even though it's serial)
    # we don't lose/overwrite particles, etc.
    n_threads = 4
    chunk_count_multiplier = 1
    n_chunks = n_threads * chunk_count_multiplier

    T_wall = 300.0
    v_wall = 500.0
    L = 5e-4
    ndens = 5e22
    nx = 50
    ppc_sampled = 400
    merge_threshold = 180
    merge_target = 150

    Δt = 2.59e-9
    preallocation_margin_multiplier = 1.5

    n_timesteps = 50

    seed = 1234
    rng_chunks = [StableRNG(seed + i) for i in 0:n_chunks-1]

    # load particle and interaction data
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")
    interaction_data_path = joinpath(@__DIR__, "..", "data", "vhs.toml")
    interaction_data::Array{Interaction, 2} = load_interaction_data(interaction_data_path, species_data)

    # create our grid and BCs
    grid = Grid1DUniform(L, nx)
    boundaries = MaxwellWalls1D(species_data, T_wall, T_wall, -v_wall, v_wall, 1.0, 1.0)

    # split cell indices into chunks
    cell_indices = Vector(1:nx)
    cell_chunks = chunks(cell_indices; n=n_chunks)

    # init per-chunk particle vectors, particle indexers, grid particle sorters
    n_particles_chunks = [floor(Int64, ppc_sampled * length(cell_chunk) * preallocation_margin_multiplier) for cell_chunk in cell_chunks]
    particles_chunks = [[ParticleVector(n_particles)] for n_particles in n_particles_chunks]
    pia_chunks = [ParticleIndexerArray(grid.n_cells, 1) for cell_chunk in cell_chunks]
    gridsorter_chunks = [GridSortInPlace(grid, n_particles) for n_particles in n_particles_chunks]

    # this is used for moving particles between chunks after they have been sorted into grid cells
    chunk_exchanger = ChunkExchanger(cell_chunks, nx)

    # sample particles
    # Fnum * ppc_sampled = Np in cell = ndens * V_cell
    Fnum = grid.cells[1].V * ndens / ppc_sampled

    # sample particles per-chunk
    for (chunk_id, cell_chunk) in enumerate(cell_chunks)
        sample_particles_equal_weight!(rng_chunks[chunk_id], grid, particles_chunks[chunk_id][1],
                                                    pia_chunks[chunk_id],
                                                    1, species_data, ndens, T_wall, Fnum, cell_chunk)
    end
     
    # create collision structs
    collision_data = [CollisionData() for cell_chunk in cell_chunks]
    
    # create struct for computation of physical properties, sizes of pia are the same
    phys_props = PhysProps(pia_chunks[1])

    # create and estimate collision factors
    collision_factors = [create_collision_factors_array(pia, interaction_data, species_data, T_wall, Fnum)
                         for pia in pia_chunks]

    # create merging structs
    oc_chunks = [OctreeN2Merge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel, max_Nbins=6000) for cell_chunks in cell_chunks]

    # merge and compute data at t=0
    for chunk_id in 1:n_chunks
        for cell in cell_chunks[chunk_id]
            merge_octree_N2_based!(rng_chunks[chunk_id], oc_chunks[chunk_id], particles_chunks[chunk_id][1], pia_chunks[chunk_id], cell, 1, merge_target, grid)
        end
        squash_pia!(particles_chunks[chunk_id], pia_chunks[chunk_id])
        compute_props_sorted!(particles_chunks[chunk_id], pia_chunks[chunk_id], species_data, phys_props, cell_chunks[chunk_id])
    end

    # compute data at t=0
    for (chunk_id, cell_chunk) in enumerate(cell_chunks)
        compute_props_sorted!(particles_chunks[chunk_id], pia_chunks[chunk_id], species_data, phys_props, cell_chunk)
    end

    for t in 1:n_timesteps
        
        # check indexing correctness
        pia_correct = [check_pia_is_correct(pia_chunks[chunk_id], 1) for chunk_id in 1:n_chunks]
        for chunk_id in 1:n_chunks
            @test pia_correct[chunk_id] == (1,0)
        end

        index_correct = [check_unique_index(particles_chunks[chunk_id][1], pia_chunks[chunk_id], 1) for chunk_id in 1:n_chunks]
        for chunk_id in 1:n_chunks
            @test index_correct[chunk_id] == (1,0)
        end
        # collide, convect, sort particles
        for chunk_id in 1:n_chunks
            for cell in cell_chunks[chunk_id]
                ntc!(rng_chunks[chunk_id], collision_factors[chunk_id][1, 1, cell],
                               collision_data[chunk_id], interaction_data, particles_chunks[chunk_id][1],
                               pia_chunks[chunk_id], cell, 1, Δt, grid.cells[cell].V)

                if pia_chunks[chunk_id].indexer[cell,1].n_local > merge_threshold
                    merge_octree_N2_based!(rng_chunks[chunk_id], oc_chunks[chunk_id], particles_chunks[chunk_id][1], pia_chunks[chunk_id], cell, 1, merge_target, grid)
                    squash_pia!(particles_chunks[chunk_id], pia_chunks[chunk_id])
                end
            end

            convect_particles!(rng_chunks[chunk_id], grid, boundaries,
                                particles_chunks[chunk_id][1], pia_chunks[chunk_id],
                                1, species_data, Δt)
        

            # need to clear the data in the chunk exchanger
            reset!(chunk_exchanger, chunk_id)

            # sort particles
            sort_particles!(gridsorter_chunks[chunk_id], grid, particles_chunks[chunk_id][1], pia_chunks[chunk_id], 1)
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

        # check that total number density is not lost
        @test abs(sum(phys_props.n) - ndens * L) / (ndens * L) < 4*eps() 
    end
end