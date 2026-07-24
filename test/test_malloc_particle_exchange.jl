@testset "malloc: particle exchange between chunks + re-sort" begin

    # 4 cells, 3 chunks
    # [1], [2,3], [4]
    # particle cell-ownership per chunk: [1,1,2,3,4,4], [1,3,4], [1,2,3,3]
    # arrays are sized generously so that no re-allocation happens during exchange / sort
    n_chunks = 3
    cell_chunks = [[1], [2,3], [4]]
    n_cells = 4

    np_actual = [6, 3, 4]
    positions = [[0.0, 0.5, 2.0, 3.0, 4.5, 5.0], [1.5, 3.5, 5.5], [-1.0, 2.5, 4.0, 4.5]]
    cells = [[1,1,2,3,4,4],[1,3,4],[1,2,3,3]]
    np_in_cells = [[2,1,1,2], [1,0,1,1], [1,1,2,0]]  # per-chunk info
    offsets = [[1,3,4,5], [1,1,2,3], [1,2,3,3]]  # per-cell, for easier setup of pia

    # oversize everything so exchange_particles! / sort_particles_after_exchange! never resize
    prealloc = 128

    # build a fresh, fully-indexed state; called between measurements since exchange mutates state
    function build_state()
        gridsorter_chunks = [GridSortInPlace(n_cells, prealloc) for i in 1:n_chunks]
        particles_chunks = [[ParticleVector(prealloc)] for i in 1:n_chunks]
        pia_chunks = [ParticleIndexerArray(n_cells, 1) for i in 1:n_chunks]
        chunk_exchanger = ChunkExchanger(cell_chunks, n_cells)

        for chunk_id in 1:n_chunks
            for np in 1:np_actual[chunk_id]
                Merzbild.add_particle!(particles_chunks[chunk_id][1], np,
                                       cells[chunk_id][np] * 1.0, SVector{3, Float64}(chunk_id, -chunk_id, chunk_id),
                                       SVector{3, Float64}(positions[chunk_id][np], 0.5, 0.0))
            end
            pia_chunks[chunk_id].n_total[1] = np_actual[chunk_id]
            pia_chunks[chunk_id].index_last[1] = np_actual[chunk_id]

            for cell in 1:n_cells
                pia_chunks[chunk_id].indexer[cell,1].n_local = np_in_cells[chunk_id][cell]
                pia_chunks[chunk_id].indexer[cell,1].n_group1 = np_in_cells[chunk_id][cell]

                if np_in_cells[chunk_id][cell] > 0
                    pia_chunks[chunk_id].indexer[cell,1].start1 = offsets[chunk_id][cell]
                    pia_chunks[chunk_id].indexer[cell,1].end1 = offsets[chunk_id][cell] + np_in_cells[chunk_id][cell] - 1
                end
            end
        end

        return chunk_exchanger, particles_chunks, pia_chunks, gridsorter_chunks
    end

    # warmup / compile both functions
    chunk_exchanger, particles_chunks, pia_chunks, gridsorter_chunks = build_state()
    exchange_particles!(chunk_exchanger, particles_chunks, pia_chunks, cell_chunks, 1)
    for chunk_id in 1:n_chunks
        sort_particles_after_exchange!(chunk_exchanger, gridsorter_chunks[chunk_id],
                                       particles_chunks[chunk_id][1], pia_chunks[chunk_id],
                                       cell_chunks[chunk_id], 1)
    end

    # measure exchange_particles! (all-pairs form) on a fresh state
    chunk_exchanger, particles_chunks, pia_chunks, gridsorter_chunks = build_state()
    bytes_exchange = @allocated exchange_particles!(chunk_exchanger, particles_chunks, pia_chunks, cell_chunks, 1)
    @test bytes_exchange == 0

    # sort each chunk after the (already-done) exchange above; measure each
    for chunk_id in 1:n_chunks
        bytes_sort = @allocated sort_particles_after_exchange!(chunk_exchanger, gridsorter_chunks[chunk_id],
                                                               particles_chunks[chunk_id][1], pia_chunks[chunk_id],
                                                               cell_chunks[chunk_id], 1)
        @test bytes_sort == 0
    end

    # measure the pairwise exchange_particles!(..., i, j) form on a fresh state
    chunk_exchanger, particles_chunks, pia_chunks, gridsorter_chunks = build_state()
    bytes_exchange_pair = 0
    for i in 1:n_chunks-1
        for j in i+1:n_chunks
            bytes_exchange_pair += @allocated exchange_particles!(chunk_exchanger, particles_chunks, pia_chunks,
                                                                  cell_chunks, 1, i, j)
        end
    end
    @test bytes_exchange_pair == 0
end
