@testset "particle re-sorting after exchange between chunks" begin
    
    # case 1
    # 4 cells, 3 chunks
    # [1], [2,3], [4]
    # [1,1,2,3,4,4], [1,3,4], [1,2,3,3]
    n_chunks = 3
    cell_chunks = [[1], [2,3], [4]]
    n_cells = 4

    # 3 chunks, 4 cells, let's allocate for 3 particles in each chunk
    gridsorter_chunks = [GridSortInPlace(n_cells, 3) for i in 1:n_chunks]

    particles_chunks = [[ParticleVector(8)], [ParticleVector(4)], [ParticleVector(5)]]
    pia_chunks = [ParticleIndexerArray(n_cells,1), ParticleIndexerArray(n_cells,1), ParticleIndexerArray(n_cells,1)]
    np_actual = [6, 3, 4]
    positions = [[0.0, 0.5, 2.0, 3.0, 4.5, 5.0], [1.5, 3.5, 5.5], [-1.0, 2.5, 4.0, 4.5]]
    cells = [[1,1,2,3,4,4],[1,3,4],[1,2,3,3]]
    np_in_cells = [[2,1,1,2], [1,0,1,1], [1,1,2,0]]  # per-chunk info
    offsets = [[1,3,4,5], [1,1,2,3], [1,2,3,3]]  # this is for easier setup of pia

    chunk_exchanger = ChunkExchanger(cell_chunks, n_cells)
    for chunk_id in 1:n_chunks
        for np in 1:np_actual[chunk_id]
            Merzbild.add_particle!(particles_chunks[chunk_id][1], np,
                                   cells[chunk_id][np] * 1.0, [chunk_id, -chunk_id, chunk_id],
                                   [positions[chunk_id][np], 0.5, 0.0])
        end
        pia_chunks[chunk_id].n_total[1] = np_actual[chunk_id]

        for cell in 1:n_cells
            pia_chunks[chunk_id].indexer[cell,1].n_local = np_in_cells[chunk_id][cell]
            pia_chunks[chunk_id].indexer[cell,1].n_group1 = np_in_cells[chunk_id][cell]

            if np_in_cells[chunk_id][cell] > 0
                pia_chunks[chunk_id].indexer[cell,1].start1 = offsets[chunk_id][cell]
                pia_chunks[chunk_id].indexer[cell,1].end1 = offsets[chunk_id][cell] + np_in_cells[chunk_id][cell] - 1
            end
        end
    end

    exchange_particles!(chunk_exchanger, particles_chunks, pia_chunks, cell_chunks, 1)

    for chunk_id in 1:n_chunks
        sort_particles_after_exchange!(chunk_exchanger, gridsorter_chunks[chunk_id],
                                       particles_chunks[chunk_id][1], pia_chunks[chunk_id],
                                       cell_chunks[chunk_id], 1)
    end

    np_actual = [4,6,3]  # this does not include particles pushed to another chunk anymore

    # n_total was the old size, including the particles sent away
    length_indices = [6 + Merzbild.DELTA_PARTICLES, 6 + Merzbild.DELTA_PARTICLES, 5 + Merzbild.DELTA_PARTICLES]
    for chunk_id in 1:3
        @test pia_chunks[chunk_id].n_total[1] == np_actual[chunk_id]
        @test length(gridsorter_chunks[chunk_id].sorted_indices) == length_indices[chunk_id]

        # test that pia empty for cells not in chunk
        for cell in 1:n_cells
            if !(cell in cell_chunks[chunk_id])
                @test pia_chunks[chunk_id].indexer[cell,1].n_local == 0
                @test pia_chunks[chunk_id].indexer[cell,1].n_group1 == 0
                @test pia_chunks[chunk_id].indexer[cell,1].start1 == 0
                @test pia_chunks[chunk_id].indexer[cell,1].end1 == -1
            end

            # group 2 is always empty
            @test pia_chunks[chunk_id].indexer[cell,1].n_group2 == 0
            @test pia_chunks[chunk_id].indexer[cell,1].start2 == 0
            @test pia_chunks[chunk_id].indexer[cell,1].end2 == -1
        end 
    end

    for chunk_id in 1:n_chunks
        @test check_unique_index(particles_chunks[chunk_id][1], pia_chunks[chunk_id], 1) == (true, 0)
    end

    # chunk 1 cell 1
    chunk_id = 1
    cell = 1
    @test pia_chunks[chunk_id].indexer[cell,1].n_local == 4
    @test pia_chunks[chunk_id].indexer[cell,1].n_group1 == 4
    @test pia_chunks[chunk_id].indexer[cell,1].start1 == 1
    @test pia_chunks[chunk_id].indexer[cell,1].end1 == 4

    for i in pia_chunks[chunk_id].indexer[cell,1].start1:pia_chunks[chunk_id].indexer[cell,1].end1
        @test particles_chunks[chunk_id][1][i].w == 1.0
    end

    # chunk 2 cells 2,3
    chunk_id = 2
    cell = 2
    @test pia_chunks[chunk_id].indexer[cell,1].n_local == 2
    @test pia_chunks[chunk_id].indexer[cell,1].n_group1 == 2
    @test pia_chunks[chunk_id].indexer[cell,1].start1 == 1
    @test pia_chunks[chunk_id].indexer[cell,1].end1 == 2

    for i in pia_chunks[chunk_id].indexer[cell,1].start1:pia_chunks[chunk_id].indexer[cell,1].end1
        @test particles_chunks[chunk_id][1][i].w == 2.0
    end

    cell = 3
    @test pia_chunks[chunk_id].indexer[cell,1].n_local == 4
    @test pia_chunks[chunk_id].indexer[cell,1].n_group1 == 4
    @test pia_chunks[chunk_id].indexer[cell,1].start1 == 3
    @test pia_chunks[chunk_id].indexer[cell,1].end1 == 6

    for i in pia_chunks[chunk_id].indexer[cell,1].start1:pia_chunks[chunk_id].indexer[cell,1].end1
        @test particles_chunks[chunk_id][1][i].w == 3.0
    end
    
    # chunk 3 cell 4
    chunk_id = 3
    cell = 4
    @test pia_chunks[chunk_id].indexer[cell,1].n_local == 3
    @test pia_chunks[chunk_id].indexer[cell,1].n_group1 == 3
    @test pia_chunks[chunk_id].indexer[cell,1].start1 == 1
    @test pia_chunks[chunk_id].indexer[cell,1].end1 == 3

    for i in pia_chunks[chunk_id].indexer[cell,1].start1:pia_chunks[chunk_id].indexer[cell,1].end1
        @test particles_chunks[chunk_id][1][i].w == 4.0
    end

    # buffer preserved
    # our initial particle array was of length 8, so buffer was of length 2 before
    # we pushed particles
    @test particles_chunks[1][1].nbuffer == 4
    @test particles_chunks[1][1].buffer[1:4] == [8,7,4,6]

    # we increased size of array
    @test particles_chunks[2][1].nbuffer == Merzbild.DELTA_PARTICLES

    # had pv of length 4 and 3 particles, so buffer was of size 1 before we pushed
    @test particles_chunks[3][1].nbuffer == 2
    @test particles_chunks[3][1].buffer[1:2] == [3,4]

    # we also test that re-writing local particle data works OK, no shared Particle instances
    for chunk_id in 1:3
        for cell in 1:4
            s1 = pia_chunks[chunk_id].indexer[cell,1].start1
            e1 = pia_chunks[chunk_id].indexer[cell,1].end1
            for i in s1:e1
                particles_chunks[chunk_id][1][i].w = chunk_id
            end
        end
    end

    for chunk_id in 1:3
        for cell in 1:4
            s1 = pia_chunks[chunk_id].indexer[cell,1].start1
            e1 = pia_chunks[chunk_id].indexer[cell,1].end1
            for i in s1:e1
                @test particles_chunks[chunk_id][1][i].w == chunk_id
            end
        end
    end


    # case 2
    # 4 cells, 2 chunks
    # [2,3,4], [1,1,1,1]
    # becomes [1,1,1,1], [2,3,4]
    n_chunks = 2
    cell_chunks = [[1], [2,3,4]]
    n_cells = 4
    gridsorter_chunks = [GridSortInPlace(n_cells, 3) for i in 1:n_chunks]

    particles_chunks = [[ParticleVector(3)], [ParticleVector(4)]]
    pia_chunks = [ParticleIndexerArray(n_cells,1), ParticleIndexerArray(n_cells,1)]
    np_actual = [3, 4]
    positions = [[2.0, 3.0, 4.0], [-1.0, -0.5, 0.5, 1.0]]
    cells = [[2,3,4],[1,1,1,1]]
    np_in_cells = [[0,1,1,1], [4,0,0,0]]  # per-chunk info
    offsets = [[0,1,2,3], [1,0,0,0]]  # per-cell, this is for easier setup of pia

    chunk_exchanger = ChunkExchanger(cell_chunks, n_cells)
    for chunk_id in 1:n_chunks
        for np in 1:np_actual[chunk_id]
            Merzbild.add_particle!(particles_chunks[chunk_id][1], np,
                                   cells[chunk_id][np] * 1.0, [chunk_id, -chunk_id, chunk_id],
                                   [positions[chunk_id][np], 0.5, 0.0])
        end
        pia_chunks[chunk_id].n_total[1] = np_actual[chunk_id]

        for cell in 1:n_cells
            pia_chunks[chunk_id].indexer[cell,1].n_local = np_in_cells[chunk_id][cell]
            pia_chunks[chunk_id].indexer[cell,1].n_group1 = np_in_cells[chunk_id][cell]

            if np_in_cells[chunk_id][cell] > 0
                pia_chunks[chunk_id].indexer[cell,1].start1 = offsets[chunk_id][cell]
                pia_chunks[chunk_id].indexer[cell,1].end1 = offsets[chunk_id][cell] + np_in_cells[chunk_id][cell] - 1
            end
        end
    end

    exchange_particles!(chunk_exchanger, particles_chunks, pia_chunks, cell_chunks, 1)

    for chunk_id in 1:n_chunks
        sort_particles_after_exchange!(chunk_exchanger, gridsorter_chunks[chunk_id],
                                       particles_chunks[chunk_id][1], pia_chunks[chunk_id],
                                       cell_chunks[chunk_id], 1)
    end

    new_lengths = [4 + Merzbild.DELTA_PARTICLES, 4]
    np_actual = [4,3]  # this does not include particles pushed to another chunk anymore

    length_indices = [4 + Merzbild.DELTA_PARTICLES, 4 + Merzbild.DELTA_PARTICLES]
    for chunk_id in 1:n_chunks
        @test pia_chunks[chunk_id].n_total[1] == np_actual[chunk_id]
        @test length(gridsorter_chunks[chunk_id].sorted_indices) == length_indices[chunk_id]

        # test that pia empty for cells not in chunk
        for cell in 1:n_cells
            if !(cell in cell_chunks[chunk_id])
                @test pia_chunks[chunk_id].indexer[cell,1].n_local == 0
                @test pia_chunks[chunk_id].indexer[cell,1].n_group1 == 0
                @test pia_chunks[chunk_id].indexer[cell,1].start1 == 0
                @test pia_chunks[chunk_id].indexer[cell,1].end1 == -1
            end

            # group 2 is always empty
            @test pia_chunks[chunk_id].indexer[cell,1].n_group2 == 0
            @test pia_chunks[chunk_id].indexer[cell,1].start2 == 0
            @test pia_chunks[chunk_id].indexer[cell,1].end2 == -1
        end 
    end

    for chunk_id in 1:n_chunks
        @test check_unique_index(particles_chunks[chunk_id][1], pia_chunks[chunk_id], 1) == (true, 0)
    end

    # chunk 1 cell 1
    chunk_id = 1
    cell = 1
    @test pia_chunks[chunk_id].indexer[cell,1].n_local == 4
    @test pia_chunks[chunk_id].indexer[cell,1].n_group1 == 4
    @test pia_chunks[chunk_id].indexer[cell,1].start1 == 1
    @test pia_chunks[chunk_id].indexer[cell,1].end1 == 4

    for i in pia_chunks[chunk_id].indexer[cell,1].start1:pia_chunks[chunk_id].indexer[cell,1].end1
        @test particles_chunks[chunk_id][1][i].w == 1.0
    end

    # chunk 2 cells 2,3,4
    chunk_id = 2
    cell = 2
    @test pia_chunks[chunk_id].indexer[cell,1].n_local == 1
    @test pia_chunks[chunk_id].indexer[cell,1].n_group1 == 1
    @test pia_chunks[chunk_id].indexer[cell,1].start1 == 1
    @test pia_chunks[chunk_id].indexer[cell,1].end1 == 1

    for i in pia_chunks[chunk_id].indexer[cell,1].start1:pia_chunks[chunk_id].indexer[cell,1].end1
        @test particles_chunks[chunk_id][1][i].w == 2.0
    end

    cell = 3
    @test pia_chunks[chunk_id].indexer[cell,1].n_local == 1
    @test pia_chunks[chunk_id].indexer[cell,1].n_group1 == 1
    @test pia_chunks[chunk_id].indexer[cell,1].start1 == 2
    @test pia_chunks[chunk_id].indexer[cell,1].end1 == 2

    for i in pia_chunks[chunk_id].indexer[cell,1].start1:pia_chunks[chunk_id].indexer[cell,1].end1
        @test particles_chunks[chunk_id][1][i].w == 3.0
    end

    cell = 4
    @test pia_chunks[chunk_id].indexer[cell,1].n_local == 1
    @test pia_chunks[chunk_id].indexer[cell,1].n_group1 == 1
    @test pia_chunks[chunk_id].indexer[cell,1].start1 == 3
    @test pia_chunks[chunk_id].indexer[cell,1].end1 == 3

    for i in pia_chunks[chunk_id].indexer[cell,1].start1:pia_chunks[chunk_id].indexer[cell,1].end1
        @test particles_chunks[chunk_id][1][i].w == 4.0
    end

    # we also test that re-writing local particle data works OK, no shared Particle instances
    for chunk_id in 1:2
        for cell in 1:4
            s1 = pia_chunks[chunk_id].indexer[cell,1].start1
            e1 = pia_chunks[chunk_id].indexer[cell,1].end1
            for i in s1:e1
                particles_chunks[chunk_id][1][i].w = chunk_id
            end
        end
    end

    for chunk_id in 1:2
        for cell in 1:4
            s1 = pia_chunks[chunk_id].indexer[cell,1].start1
            e1 = pia_chunks[chunk_id].indexer[cell,1].end1
            for i in s1:e1
                @test particles_chunks[chunk_id][1][i].w == chunk_id
            end
        end
    end

    # case 3
    n_chunks = 2
    cell_chunks = [[1], [2]]
    n_cells = 2
    gridsorter_chunks = [GridSortInPlace(n_cells, 5) for i in 1:n_chunks]

    particles_chunks = [[ParticleVector(5)], [ParticleVector(5)]]
    pia_chunks = [ParticleIndexerArray(n_cells,1), ParticleIndexerArray(n_cells,1)]
    np_actual = [3, 4]
    positions = [[2.0, 2.0, 2.0], [-0.5, -0.5, 0.5, -0.5]]
    cells = [[2,2,2],[1,1,1,1]]
    np_in_cells = [[0,3], [4,0]]  # per-chunk info
    offsets = [[0,1], [1,0]]  # per-cell, this is for easier setup of pia

    chunk_exchanger = ChunkExchanger(cell_chunks, n_cells)
    for chunk_id in 1:n_chunks
        for np in 1:np_actual[chunk_id]
            Merzbild.add_particle!(particles_chunks[chunk_id][1], np,
                                   cells[chunk_id][np] * 1.0, [chunk_id, -chunk_id, chunk_id],
                                   [positions[chunk_id][np], 0.5, 0.0])
        end
        pia_chunks[chunk_id].n_total[1] = np_actual[chunk_id]

        for cell in 1:n_cells
            pia_chunks[chunk_id].indexer[cell,1].n_local = np_in_cells[chunk_id][cell]
            pia_chunks[chunk_id].indexer[cell,1].n_group1 = np_in_cells[chunk_id][cell]

            if np_in_cells[chunk_id][cell] > 0
                pia_chunks[chunk_id].indexer[cell,1].start1 = offsets[chunk_id][cell]
                pia_chunks[chunk_id].indexer[cell,1].end1 = offsets[chunk_id][cell] + np_in_cells[chunk_id][cell] - 1
            end
        end
    end

    @test pia_chunks[1].n_total[1] == 3
    @test pia_chunks[2].n_total[1] == 4

    exchange_particles!(chunk_exchanger, particles_chunks, pia_chunks, cell_chunks, 1)

    for chunk_id in 1:n_chunks
        sort_particles_after_exchange!(chunk_exchanger, gridsorter_chunks[chunk_id],
                                       particles_chunks[chunk_id][1], pia_chunks[chunk_id],
                                       cell_chunks[chunk_id], 1)
    end

    @test pia_chunks[1].n_total[1] == 4
    @test pia_chunks[2].n_total[1] == 3

    # add a new particle to chunk 2
    chunk_id = 2
    Merzbild.update_buffer_index_new_particle!(particles_chunks[chunk_id][1], pia_chunks[chunk_id], 2, 1)

    particles_chunks[chunk_id][1][pia_chunks[chunk_id].n_total[1]].w = 2.0
    println(pia_chunks[chunk_id].n_total[1])
    for chunk_id in 1:2
        np = 0
        for cell in 1:2
            s1 = pia_chunks[chunk_id].indexer[cell,1].start1
            e1 = pia_chunks[chunk_id].indexer[cell,1].end1
            for i in s1:e1
                @test particles_chunks[chunk_id][1][i].w == chunk_id
                np += 1
            end

            s2 = pia_chunks[chunk_id].indexer[cell,1].start2
            e2 = pia_chunks[chunk_id].indexer[cell,1].end2
            if pia_chunks[chunk_id].indexer[cell,1].n_group2 > 0
                for i in s2:e2
                    @test particles_chunks[chunk_id][1][i].w == chunk_id
                    np += 1
                end
            end
        end
        @test np == 4
    end
end