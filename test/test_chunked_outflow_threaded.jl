function run_chunked_outflow(n_chunks, n_timesteps; threaded=false)
    T_wall = 300.0
    L = 5e-4
    ndens = 5e22
    nx = 50
    ppc = 400
    Δt = 2.59e-9
    preallocation_margin_multiplier = 1.5

    seed = 1234
    rng_chunks = [StableRNG(seed + i) for i in 0:n_chunks-1]

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")

    grid = Grid1DUniform(L, nx)

    # the right wall lets particles out of the domain
    bc_list = (FullyDiffuseBC1D(1, species_data, T_wall, [0.0, 0.0, 0.0]),
               OutflowBC())

    cell_indices = Vector(1:nx)
    cell_chunks = chunks(cell_indices; n=n_chunks)

    n_particles_chunks = [floor(Int64, ppc * length(cell_chunk) * preallocation_margin_multiplier)
                          for cell_chunk in cell_chunks]
    particles_chunks = [[ParticleVector(n_particles)] for n_particles in n_particles_chunks]
    pia_chunks = [ParticleIndexerArray(grid.n_cells, 1) for cell_chunk in cell_chunks]
    gridsorter_chunks = [GridSortInPlace(grid, n_particles) for n_particles in n_particles_chunks]

    chunk_exchanger = ChunkExchanger(cell_chunks, nx)

    Fnum = grid.cells[1].V * ndens / ppc

    for (chunk_id, cell_chunk) in enumerate(cell_chunks)
        sample_particles_equal_weight!(rng_chunks[chunk_id], grid, particles_chunks[chunk_id][1],
                                       pia_chunks[chunk_id], 1, species_data, ndens, T_wall, Fnum,
                                       cell_chunk)
    end

    phys_props = PhysProps(pia_chunks[1])

    n_particles_initial = sum(pia_chunks[chunk_id].n_total[1] for chunk_id in 1:n_chunks)
    n_deleted_chunks = zeros(Int64, n_chunks)
    n_deleted_total = 0

    function convect_and_sort!(chunk_id)
        @inbounds n_deleted_chunks[chunk_id] =
            convect_particles!(rng_chunks[chunk_id], grid, bc_list,
                               particles_chunks[chunk_id][1], pia_chunks[chunk_id],
                               1, species_data, Δt)

        @inbounds sort_particles!(gridsorter_chunks[chunk_id], grid, particles_chunks[chunk_id][1],
                                  pia_chunks[chunk_id], 1)
        @inbounds update_occupancy_bounds!(chunk_exchanger, gridsorter_chunks[chunk_id],
                                           pia_chunks[chunk_id], chunk_id, 1)
        return nothing
    end

    function restore_indexing_and_compute_props!(chunk_id)
        @inbounds sort_particles_after_exchange!(chunk_exchanger, gridsorter_chunks[chunk_id],
                                                 particles_chunks[chunk_id][1], pia_chunks[chunk_id],
                                                 cell_chunks[chunk_id], 1)
        @inbounds compute_props_sorted!(particles_chunks[chunk_id], pia_chunks[chunk_id],
                                        species_data, phys_props, cell_chunks[chunk_id])
        return nothing
    end

    for t in 1:n_timesteps
        if threaded
            @threads for chunk_id in 1:n_chunks
                convect_and_sort!(chunk_id)
            end
        else
            for chunk_id in 1:n_chunks
                convect_and_sort!(chunk_id)
            end
        end

        n_deleted_total += sum(n_deleted_chunks)

        exchange_particles!(chunk_exchanger, particles_chunks, pia_chunks, cell_chunks, 1)

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
            n_particles_initial=n_particles_initial,
            n_deleted=n_deleted_total,
            n_total=[pia_chunks[chunk_id].n_total[1] for chunk_id in 1:n_chunks],
            pia_correct=[check_pia_is_correct(pia_chunks[chunk_id], 1) for chunk_id in 1:n_chunks],
            index_correct=[check_unique_index(particles_chunks[chunk_id][1], pia_chunks[chunk_id], 1) for chunk_id in 1:n_chunks],
            buffer_correct=[check_unique_buffer(particles_chunks[chunk_id][1]) for chunk_id in 1:n_chunks],
            owns_only_its_cells=[all(pia_chunks[chunk_id].indexer[cell, 1].n_local == 0
                                     for cell in 1:nx if !(cell in cell_chunks[chunk_id]))
                                 for chunk_id in 1:n_chunks])
end

@testset "chunked convection with an outflow boundary" begin
    n_chunks = 4
    n_timesteps = 20

    if Threads.nthreads() == 1
        @info "chunked outflow threading test running on 1 thread: no concurrency is exercised"
    end

    serial = run_chunked_outflow(n_chunks, n_timesteps; threaded=false)

    # particles are deleted only at the outflow boundary: everything else is still in the domain
    @test serial.n_deleted > 0
    @test sum(serial.n_total) + serial.n_deleted == serial.n_particles_initial
    @test sum(serial.n_total) < serial.n_particles_initial

    # deletion leaves every chunk's indexing self-consistent across the whole
    # convect -> sort -> exchange -> sort_after_exchange cycle
    for chunk_id in 1:n_chunks
        @test serial.pia_correct[chunk_id] == (true, 0)
        @test serial.index_correct[chunk_id] == (true, 0)
        @test serial.buffer_correct[chunk_id] == (true, 0)
        @test serial.owns_only_its_cells[chunk_id] == true
    end

    # each chunk deletes only from its own particle vector and indexer, so threading may not
    # change a single bit of the answer
    threaded = run_chunked_outflow(n_chunks, n_timesteps; threaded=true)

    @test threaded.n_deleted == serial.n_deleted
    @test threaded.n_total == serial.n_total
    @test threaded.phys_props.n == serial.phys_props.n
    @test threaded.phys_props.v == serial.phys_props.v
    @test threaded.phys_props.T == serial.phys_props.T

    for chunk_id in 1:n_chunks
        @test threaded.pia_correct[chunk_id] == (true, 0)
        @test threaded.index_correct[chunk_id] == (true, 0)
        @test threaded.buffer_correct[chunk_id] == (true, 0)
        @test threaded.owns_only_its_cells[chunk_id] == true
    end
end
