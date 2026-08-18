function run_chunked_couette(n_chunks, n_timesteps; threaded, factorized_exchange)
    T_wall = 300.0
    v_wall = 500.0
    L = 5e-4
    ndens = 5e22
    nx = 40
    ppc = 200
    Δt = 2.59e-9
    preallocation_margin_multiplier = 1.5

    seed = 1234
    rng_chunks = [StableRNG(seed + i) for i in 0:n_chunks-1]

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")
    interaction_data_path = joinpath(@__DIR__, "..", "data", "vhs.toml")
    interaction_data::Array{Interaction, 2} = load_interaction_data(interaction_data_path, species_data)

    grid = Grid1DUniform(L, nx)
    bc_list = (MaxwellWallBC1D(1, species_data, T_wall, [0.0, -v_wall, 0.0], 1.0),
               MaxwellWallBC1D(1, species_data, T_wall, [0.0, v_wall, 0.0], 1.0))

    cell_indices = Vector(1:nx)
    cell_chunks = chunks(cell_indices; n=n_chunks)

    n_particles_chunks = [floor(Int64, ppc * length(cell_chunk) * preallocation_margin_multiplier) for cell_chunk in cell_chunks]
    particles_chunks = [[ParticleVector{1}(n_particles)] for n_particles in n_particles_chunks]
    pia_chunks = [ParticleIndexerArray(grid.n_cells, 1) for cell_chunk in cell_chunks]
    gridsorter_chunks = [GridSortInPlace(grid, n_particles) for n_particles in n_particles_chunks]

    chunk_exchanger = ChunkExchanger(cell_chunks, nx)

    Fnum = grid.cells[1].V * ndens / ppc

    for (chunk_id, cell_chunk) in enumerate(cell_chunks)
        sample_particles_equal_weight!(rng_chunks[chunk_id], grid, particles_chunks[chunk_id][1],
                                       pia_chunks[chunk_id], 1, species_data, ndens, T_wall, Fnum, cell_chunk)
    end

    collision_data = [CollisionData() for cell_chunk in cell_chunks]
    collision_factors = [create_collision_factors_array(pia, interaction_data, species_data, T_wall, Fnum)
                         for pia in pia_chunks]

    phys_props = PhysProps(pia_chunks[1])
    exchange_list = generate_1_factorization(n_chunks)

    # collide, convect and sort the particles owned by one chunk; identical code runs in the
    # serial and threaded variants so that any difference in the results comes from concurrency
    function advance_chunk!(chunk_id)
        @inbounds for cell in cell_chunks[chunk_id]
            ntc_equal_weight!(rng_chunks[chunk_id], collision_factors[chunk_id][1, 1, cell],
                              collision_data[chunk_id], interaction_data, particles_chunks[chunk_id][1],
                              pia_chunks[chunk_id], cell, 1, Δt, grid.cells[cell].V)
        end

        convect_particles!(rng_chunks[chunk_id], grid, bc_list, particles_chunks[chunk_id][1],
                           pia_chunks[chunk_id], 1, species_data, Δt)

        sort_particles!(gridsorter_chunks[chunk_id], grid, particles_chunks[chunk_id][1], pia_chunks[chunk_id], 1)
    end

    function restore_indexing_and_compute_props!(chunk_id)
        sort_particles_after_exchange!(chunk_exchanger, gridsorter_chunks[chunk_id],
                                       particles_chunks[chunk_id][1], pia_chunks[chunk_id],
                                       cell_chunks[chunk_id], 1)
        compute_props_sorted!(particles_chunks[chunk_id], pia_chunks[chunk_id],
                              species_data, phys_props, cell_chunks[chunk_id])
    end

    for t in 1:n_timesteps
        if threaded
            @threads for chunk_id in 1:n_chunks
                advance_chunk!(chunk_id)
            end
        else
            for chunk_id in 1:n_chunks
                advance_chunk!(chunk_id)
            end
        end

        if factorized_exchange
            for exchange_partial in exchange_list
                if threaded
                    @threads for exchange_pair in exchange_partial
                        exchange_particles!(chunk_exchanger, particles_chunks, pia_chunks, cell_chunks,
                                            1, exchange_pair[1], exchange_pair[2])
                    end
                else
                    for exchange_pair in exchange_partial
                        exchange_particles!(chunk_exchanger, particles_chunks, pia_chunks, cell_chunks,
                                            1, exchange_pair[1], exchange_pair[2])
                    end
                end
            end
        else
            exchange_particles!(chunk_exchanger, particles_chunks, pia_chunks, cell_chunks, 1)
        end

        if threaded
            @threads for chunk_id in 1:n_chunks
                restore_indexing_and_compute_props!(chunk_id)
            end
        else
            for chunk_id in 1:n_chunks
                restore_indexing_and_compute_props!(chunk_id)
            end
        end
    end

    return (phys_props=phys_props,
            n_total=[pia_chunks[chunk_id].n_total[1] for chunk_id in 1:n_chunks],
            pia_correct=[check_pia_is_correct(pia_chunks[chunk_id], 1) for chunk_id in 1:n_chunks],
            index_correct=[check_unique_index(particles_chunks[chunk_id][1], pia_chunks[chunk_id], 1) for chunk_id in 1:n_chunks],
            buffer_correct=[check_unique_buffer(particles_chunks[chunk_id][1]) for chunk_id in 1:n_chunks],
            ndens_total=ndens * L)
end

@testset "chunked couette: threaded results identical to serial" begin
    n_chunks = 4
    n_timesteps = 50

    if Threads.nthreads() == 1
        @info "chunked couette threading test running on 1 thread: no concurrency is exercised"
    end

    # each chunk carries its own RNG and sort_particles_after_exchange! restores a canonical
    # cell-sorted ordering every timestep, so the results depend only on the number of chunks:
    # neither the thread count nor the order in which chunk pairs exchange particles may change
    # a single bit of the answer
    serial = run_chunked_couette(n_chunks, n_timesteps; threaded=false, factorized_exchange=false)

    for (threaded, factorized_exchange) in ((true, false), (false, true), (true, true))
        result = run_chunked_couette(n_chunks, n_timesteps;
                                     threaded=threaded, factorized_exchange=factorized_exchange)

        @test result.phys_props.np == serial.phys_props.np
        @test result.phys_props.n == serial.phys_props.n
        @test result.phys_props.v == serial.phys_props.v
        @test result.phys_props.T == serial.phys_props.T
        @test result.n_total == serial.n_total

        for chunk_id in 1:n_chunks
            @test result.pia_correct[chunk_id] == (1, 0)
            @test result.index_correct[chunk_id] == (1, 0)
            @test result.buffer_correct[chunk_id] == (1, 0)
        end

        @test abs(sum(result.phys_props.n) - result.ndens_total) / result.ndens_total < 4*eps()
    end
end
