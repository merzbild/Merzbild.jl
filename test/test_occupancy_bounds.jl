@testset "occupancy bounds for chunk exchange" begin
    # set up a single chunk holding particles in the given cells, sort it, and check that
    # update_occupancy_bounds! finds the same first/last occupied cell as a scan over the pia
    function bounds_for(cells, n_cells)
        np = length(cells)
        particles = ParticleVector(max(np, 1))
        pia = ParticleIndexerArray(n_cells, 1)
        chunk_exchanger = ChunkExchanger([1:n_cells], n_cells)
        gridsorter = GridSortInPlace(n_cells, max(np, 1))

        for i in 1:np
            Merzbild.add_particle!(particles, i, 1.0, SVector{3,Float64}(0.0, 0.0, 0.0),
                                   SVector{3,Float64}(0.0, 0.0, 0.0))
            particles.cell[i] = cells[i]
        end
        pia.n_total[1] = np
        pia.index_last[1] = np
        pia.indexer[1,1].n_local = np
        pia.indexer[1,1].n_group1 = np
        if np > 0
            pia.indexer[1,1].start1 = 1
            pia.indexer[1,1].end1 = np
        end

        sort_particles!(gridsorter, particles, pia, 1)

        update_occupancy_bounds!(chunk_exchanger, gridsorter, pia, 1, 1)  # warmup / compile
        @test (@allocated update_occupancy_bounds!(chunk_exchanger, gridsorter, pia, 1, 1)) == 0

        occupied = [cell for cell in 1:n_cells if pia.indexer[cell,1].n_group1 > 0]
        return chunk_exchanger.occ_lo[1], chunk_exchanger.occ_hi[1], occupied
    end

    n_cells = 8
    # spread over the grid, single cell, first cell, last cell, both ends, all cells
    for cells in ([3,3,4,6], [5,5,5], [1], [8], [1,8], collect(1:8), [2,7,7,2])
        lo, hi, occupied = bounds_for(cells, n_cells)
        @test lo == minimum(occupied)
        @test hi == maximum(occupied)
    end

    # an empty chunk must produce an empty range, so that every pair involving it is rejected
    lo, hi, occupied = bounds_for(Int64[], n_cells)
    @test length(occupied) == 0
    @test lo > hi

    # a fresh ChunkExchanger must claim the whole grid, so that an exchange performed without
    # ever calling update_occupancy_bounds! is still correct
    chunk_exchanger = ChunkExchanger([1:4, 5:8], n_cells)
    for chunk_id in 1:2
        @test chunk_exchanger.occ_lo[chunk_id] == 1
        @test chunk_exchanger.occ_hi[chunk_id] == n_cells
    end

    # the bounds must not change the outcome of an exchange: run the same exchange with
    # correct bounds and with the conservative full-grid bounds, and compare
    n_chunks = 2
    cell_chunks = [[1], [2,3,4]]
    n_cells_ex = 4
    np_actual = [3, 4]
    positions = [[2.0, 3.0, 4.0], [-1.0, -0.5, 0.5, 1.0]]
    cells_ex = [[2,3,4],[1,1,1,1]]
    np_in_cells = [[0,1,1,1], [4,0,0,0]]
    offsets = [[0,1,2,3], [1,0,0,0]]

    function run_exchange(use_bounds)
        particles_chunks = [[ParticleVector(8)] for i in 1:n_chunks]
        pia_chunks = [ParticleIndexerArray(n_cells_ex, 1) for i in 1:n_chunks]
        chunk_exchanger = ChunkExchanger(cell_chunks, n_cells_ex)

        for chunk_id in 1:n_chunks
            for np in 1:np_actual[chunk_id]
                Merzbild.add_particle!(particles_chunks[chunk_id][1], np,
                                       cells_ex[chunk_id][np] * 1.0,
                                       SVector{3,Float64}(chunk_id, -chunk_id, chunk_id),
                                       SVector{3,Float64}(positions[chunk_id][np], 0.5, 0.0))
            end
            pia_chunks[chunk_id].n_total[1] = np_actual[chunk_id]
            pia_chunks[chunk_id].index_last[1] = np_actual[chunk_id]

            for cell in 1:n_cells_ex
                pia_chunks[chunk_id].indexer[cell,1].n_local = np_in_cells[chunk_id][cell]
                pia_chunks[chunk_id].indexer[cell,1].n_group1 = np_in_cells[chunk_id][cell]

                if np_in_cells[chunk_id][cell] > 0
                    pia_chunks[chunk_id].indexer[cell,1].start1 = offsets[chunk_id][cell]
                    pia_chunks[chunk_id].indexer[cell,1].end1 =
                        offsets[chunk_id][cell] + np_in_cells[chunk_id][cell] - 1
                end
            end
        end

        if use_bounds
            # chunk 1 holds particles in cells 2:4, chunk 2 holds particles in cell 1 only
            chunk_exchanger.occ_lo[1] = 2
            chunk_exchanger.occ_hi[1] = 4
            chunk_exchanger.occ_lo[2] = 1
            chunk_exchanger.occ_hi[2] = 1
        end

        exchange_particles!(chunk_exchanger, particles_chunks, pia_chunks, cell_chunks, 1)

        return ([copy(chunk_exchanger.start1), copy(chunk_exchanger.n_group1),
                 copy(chunk_exchanger.start2), copy(chunk_exchanger.n_group2)],
                [[particles_chunks[c][1][np].w for np in 1:pia_chunks[c].n_total[1]]
                 for c in 1:n_chunks])
    end

    ce_bounded, w_bounded = run_exchange(true)
    ce_full, w_full = run_exchange(false)
    @test ce_bounded == ce_full
    @test w_bounded == w_full
end
