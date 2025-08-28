@testset "particle exchange between chunks" begin
    # test particle exchange functionality that is used in multi-threaded code
    # to move particles between chunks

    # first we test swapping particles
    particles_1 = [ParticleVector(4)]
    particles_2 = [ParticleVector(3)]

    particles_1[1][1] = Particle(10.0, [1.0, 2.0, 3.0], [4.0, 85.0, 6.0])
    particles_1[1][2] = Particle(11.0, [11.0, 2.0, 3.0], [24.0, 5.0, 6.0])
    particles_1[1][3] = Particle(12.0, [111.0, 2.0, 3.0], [34.0, 95.0, 6.0])
    particles_1[1][4] = Particle(13.0, [1111.0, 2.0, 3.0], [44.0, 105.0, 6.0])

    particles_2[1][1] = Particle(600.0, [-21.0, -4.0, -8.0], [1.0, 3.0, 16.0])
    particles_2[1][2] = Particle(500.0, [-31.0, -4.0, -7.0], [1.0, 8.0, 9.0])
    particles_2[1][3] = Particle(400.0, [-51.0, -4.0, -6.0], [2.0, 3.0, 44.0])

    Merzbild.swap_particles!(particles_1[1], particles_2[1], 2, 3)
    @test particles_1[1][1].w == 10.0
    @test particles_1[1][1].v[1] == 1.0
    @test particles_1[1][1].x[2] == 85.0

    @test particles_1[1][2].w == 400.0
    @test particles_1[1][2].v[1] == -51.0
    @test particles_1[1][2].x[2] == 3.0

    @test particles_1[1][3].w == 12.0
    @test particles_1[1][3].v[1] == 111.0
    @test particles_1[1][3].x[2] == 95.0

    @test particles_1[1][4].w == 13.0
    @test particles_1[1][4].v[1] == 1111.0
    @test particles_1[1][4].x[2] == 105.0

    @test particles_2[1][1].w == 600.0
    @test particles_2[1][1].v[1] == -21.0
    @test particles_2[1][1].x[2] == 3.0

    @test particles_2[1][2].w == 500.0
    @test particles_2[1][2].v[1] == -31.0
    @test particles_2[1][2].x[2] == 8.0

    @test particles_2[1][3].w == 11.0
    @test particles_2[1][3].v[1] == 11.0
    @test particles_2[1][3].x[2] == 5.0

    # simplest 4-cell tests
    chunks = [[1,2], [3,4]]
    chunk_exchanger = ChunkExchanger(chunks, 4)
    @test chunk_exchanger.n_chunks == 2
    @test chunk_exchanger.n_cells == 4
    @test size(chunk_exchanger.indexer) == (2,4)

    # test that reset! works
    for i in 1:4
        for j in 1:2
            chunk_exchanger.indexer[j,i].n_group1 = 105 - i
            chunk_exchanger.indexer[j,i].start1 = i+j
            chunk_exchanger.indexer[j,i].end1 = 240+j
            chunk_exchanger.indexer[j,i].n_group2 = 777
            chunk_exchanger.indexer[j,i].start2 = 3+i
            chunk_exchanger.indexer[j,i].end2 = 11+j
        end
    end

    for j in 1:2
        reset!(chunk_exchanger, j)
    end

    for i in 1:4
        for j in 1:2
            @test chunk_exchanger.indexer[j,i].n_group1 == 0
            @test chunk_exchanger.indexer[j,i].start1 == 0
            @test chunk_exchanger.indexer[j,i].end1 == -1
            @test chunk_exchanger.indexer[j,i].n_group2 == 0
            @test chunk_exchanger.indexer[j,i].start2 == 0
            @test chunk_exchanger.indexer[j,i].end2 == -1
        end
    end

    # now we test exchange of particles between chunks
    # scenario 1: chunk 1 is "responsible" for cell 1
    # chunk 2 is "responsible" for cell 2
    # all particles in chunk 1 belong to cell 1
    # and all particles in chunk 2 belong to cell 2
    # so nothing happens
    cell_chunks = [[1], [2]]
    n_cells = 2

    chunk_exchanger = ChunkExchanger(cell_chunks, n_cells)

    particles_chunks = [[ParticleVector(4)], [ParticleVector(4)]]
    pia_chunks = [ParticleIndexerArray(2,1), ParticleIndexerArray(2,1)]

    for chunk_id in 1:2
        for np in 1:4
            Merzbild.add_particle!(particles_chunks[chunk_id][1], np,
                                   chunk_id * 1.0, [chunk_id^3 * np, -np + 1.0, np], [chunk_id * 3.0, 0.5, 1.0+chunk_id])
        end
        pia_chunks[chunk_id].n_total[1] = 4
        pia_chunks[chunk_id].indexer[chunk_id,1].n_local = 4  # here cell is same as chunk_id
        pia_chunks[chunk_id].indexer[chunk_id,1].n_group1 = 4  # here cell is same as chunk_id
        pia_chunks[chunk_id].indexer[chunk_id,1].start1 = 1  # here cell is same as chunk_id
        pia_chunks[chunk_id].indexer[chunk_id,1].end1 = 4  # here cell is same as chunk_id
    end
    
    for j in 1:2
        reset!(chunk_exchanger, j)
    end
    exchange_particles!(chunk_exchanger, particles_chunks, pia_chunks, cell_chunks, 1)
    
    # test that array sizes stayed the same, particles remain unchanged
    for chunk_id in 1:2
        @test length(particles_chunks[chunk_id][1]) == 4
        for np in 1:4
            @test particles_chunks[chunk_id][1][np].w == chunk_id
            @test particles_chunks[chunk_id][1][np].v[1] == chunk_id^3 * np
            @test particles_chunks[chunk_id][1][np].x[1] == chunk_id * 3.0
            @test particles_chunks[chunk_id][1][np].x[3] == 1.0 + chunk_id
        end
    end

    # same, 2 chunks, 1,1,2,2; 1,1,2,2 - particles should be swapped
    # imagine cells with bounds [-3, 2],[2,7]
    particles_chunks = [[ParticleVector(4)], [ParticleVector(4)]]
    pia_chunks = [ParticleIndexerArray(2,1), ParticleIndexerArray(2,1)]
    positions = [[0.0, 1.0, 3.0, 4.0], [-2.0, -1.0, 6.0, 5.0]]
    for chunk_id in 1:2
        for np in 1:4
            Merzbild.add_particle!(particles_chunks[chunk_id][1], np,
                                   chunk_id * 1.0, [chunk_id^3 * np, -np + 1.0, np],
                                   [positions[chunk_id][np], 0.5, 1.0+chunk_id])
        end
        pia_chunks[chunk_id].n_total[1] = 4
    end
    for chunk_id in 1:2
        # the indexing is the same for both chunks
        pia_chunks[chunk_id].indexer[1,1].n_local = 2
        pia_chunks[chunk_id].indexer[1,1].n_group1 = 2
        pia_chunks[chunk_id].indexer[1,1].start1 = 1
        pia_chunks[chunk_id].indexer[1,1].end1 = 2
        pia_chunks[chunk_id].indexer[2,1].n_local = 2
        pia_chunks[chunk_id].indexer[2,1].n_group1 = 2
        pia_chunks[chunk_id].indexer[2,1].start1 = 3
        pia_chunks[chunk_id].indexer[2,1].end1 = 4
    end

    for j in 1:2
        reset!(chunk_exchanger, j)
    end
    exchange_particles!(chunk_exchanger, particles_chunks, pia_chunks, cell_chunks, 1)

    new_positions = [[0.0, 1.0, -2.0, -1.0], [3.0, 4.0, 6.0, 5.0]]
    new_weights = [[1.0, 1.0, 2.0, 2.0], [1.0, 1.0, 2.0, 2.0]]
    for chunk_id in 1:2
        @test length(particles_chunks[chunk_id][1]) == 4
        for np in 1:4
            @test particles_chunks[chunk_id][1][np].x[1] == new_positions[chunk_id][np]
            @test particles_chunks[chunk_id][1][np].w == new_weights[chunk_id][np]
        end
    end
    # test that indexing points to the new groups of particles
    # new particles in cell 1 came from chunk 2
    @test chunk_exchanger.indexer[2,1].n_group1 == 2
    @test chunk_exchanger.indexer[2,1].start1 == 3
    @test chunk_exchanger.indexer[2,1].end1 == 4
    # particles don't come from the chunk themselves, obviously
    @test chunk_exchanger.indexer[1,1].start1 == 0
    @test chunk_exchanger.indexer[1,1].end1 == -1

    # new particles in cell 2 came from chunk 1
    @test chunk_exchanger.indexer[1,2].n_group1 == 2
    @test chunk_exchanger.indexer[1,2].start1 == 1
    @test chunk_exchanger.indexer[1,2].end1 == 2
    # particles don't come from the chunk themselves, obviously
    @test chunk_exchanger.indexer[2,2].start1 == 0
    @test chunk_exchanger.indexer[2,2].end1 == -1

    for chunk_id in 1:2
        for cell in 1:2
            @test chunk_exchanger.indexer[chunk_id,cell].n_group2 == 0
            @test chunk_exchanger.indexer[chunk_id,cell].start2 == 0
            @test chunk_exchanger.indexer[chunk_id,cell].end2 == -1
        end
    end

    @test pia_chunks[1].indexer[1,1].n_local == 2
    @test pia_chunks[1].indexer[1,1].n_group1 == 2
    @test pia_chunks[1].indexer[1,1].start1 == 1
    @test pia_chunks[1].indexer[1,1].end1 == 2
    @test pia_chunks[1].indexer[2,1].n_local == 0
    @test pia_chunks[1].indexer[2,1].n_group1 == 0
    @test pia_chunks[1].indexer[2,1].start1 == 0
    @test pia_chunks[1].indexer[2,1].end1 == -1

    @test pia_chunks[2].indexer[1,1].n_local == 0
    @test pia_chunks[2].indexer[1,1].n_group1 == 0
    @test pia_chunks[2].indexer[1,1].start1 == 0
    @test pia_chunks[2].indexer[1,1].end1 == -1
    @test pia_chunks[2].indexer[2,1].n_local == 2
    @test pia_chunks[2].indexer[2,1].n_group1 == 2
    @test pia_chunks[2].indexer[2,1].start1 == 3
    @test pia_chunks[2].indexer[2,1].end1 == 4

    # 1,1,2; 1,1,1,2 - size should be increased and they should be swapped
    # this tests swap + push of remaining particles + push when stopped in the middle of the cell
    # imagine cells with bounds [-3, 2],[2,7]
    particles_chunks = [[ParticleVector(3)], [ParticleVector(4)]]
    pia_chunks = [ParticleIndexerArray(2,1), ParticleIndexerArray(2,1)]
    np_actual = [3,4]
    positions = [[0.0, 1.0, 4.0], [-2.0, -1.0, 1.5, 6.0]]

    for chunk_id in 1:2
        for np in 1:np_actual[chunk_id]
            Merzbild.add_particle!(particles_chunks[chunk_id][1], np,
                                   chunk_id * 1.0, [chunk_id^3 * np, -np + 1.0, np],
                                   [positions[chunk_id][np], 0.5, 1.0+chunk_id])
        end
        pia_chunks[chunk_id].n_total[1] = np_actual[chunk_id]
    end
    pia_chunks[1].indexer[1,1].n_local = 2
    pia_chunks[1].indexer[1,1].n_group1 = 2
    pia_chunks[1].indexer[1,1].start1 = 1
    pia_chunks[1].indexer[1,1].end1 = 2
    pia_chunks[1].indexer[2,1].n_local = 1
    pia_chunks[1].indexer[2,1].n_group1 = 1
    pia_chunks[1].indexer[2,1].start1 = 3
    pia_chunks[1].indexer[2,1].end1 = 3

    pia_chunks[2].indexer[1,1].n_local = 3
    pia_chunks[2].indexer[1,1].n_group1 = 3
    pia_chunks[2].indexer[1,1].start1 = 1
    pia_chunks[2].indexer[1,1].end1 = 3
    pia_chunks[2].indexer[2,1].n_local = 1
    pia_chunks[2].indexer[2,1].n_group1 = 1
    pia_chunks[2].indexer[2,1].start1 = 4
    pia_chunks[2].indexer[2,1].end1 = 4

    @test particles_chunks[1][1].nbuffer == 0
    @test particles_chunks[2][1].nbuffer == 0

    for j in 1:2
        reset!(chunk_exchanger, j)
    end
    exchange_particles!(chunk_exchanger, particles_chunks, pia_chunks, cell_chunks, 1)

    # we increase size by # of particles to add + DELTA_PARTICLES
    new_lengths = [5 + Merzbild.DELTA_PARTICLES, 4]
    np_actual = [5,4]

    # original positions = [[0.0, 1.0, 4.0], [-2.0, -1.0, 1.5, 6.0]]
    # we don't delete the particles that weren't swapped, so they
    # will still show up during direct iteration
    new_positions = [[0.0, 1.0, -2.0, -1.0, 1.5], [4.0, -1.0, 1.5, 6.0]]
    new_weights = [[1.0, 1.0, 2.0, 2.0, 2.0], [1.0, 2.0, 2.0, 2.0]]
    for chunk_id in 1:2
        @test length(particles_chunks[chunk_id][1]) == new_lengths[chunk_id]
        for np in 1:np_actual[chunk_id]
            @test particles_chunks[chunk_id][1][np].x[1] == new_positions[chunk_id][np]
            @test particles_chunks[chunk_id][1][np].w == new_weights[chunk_id][np]
        end
    end
    @test particles_chunks[1][1].nbuffer == Merzbild.DELTA_PARTICLES
    flag = false
    # test that only new elements are in the buffer
    for i in 1:particles_chunks[1][1].nbuffer
        if particles_chunks[1][1].buffer[i] <= 5
            flag = true
            break
        end
    end
    @test flag == false

    # test that the buffer holds the particles that were sent to another chunk
    @test particles_chunks[2][1].nbuffer == 2
    @test particles_chunks[2][1].buffer[1:2] == [2,3]

    # test that in pia we re-set everything related to particles sent to another chunk
    for chunk_id in 1:2
        @test pia_chunks[chunk_id].n_total[1] == np_actual[chunk_id]
    end
    @test pia_chunks[1].indexer[1,1].n_local == 2
    @test pia_chunks[1].indexer[1,1].n_group1 == 2
    @test pia_chunks[1].indexer[1,1].start1 == 1
    @test pia_chunks[1].indexer[1,1].end1 == 2
    @test pia_chunks[1].indexer[2,1].n_local == 0
    @test pia_chunks[1].indexer[2,1].n_group1 == 0
    @test pia_chunks[1].indexer[2,1].start1 == 0
    @test pia_chunks[1].indexer[2,1].end1 == -1

    @test pia_chunks[2].indexer[1,1].n_local == 0
    @test pia_chunks[2].indexer[1,1].n_group1 == 0
    @test pia_chunks[2].indexer[1,1].start1 == 0
    @test pia_chunks[2].indexer[1,1].end1 == -1
    @test pia_chunks[2].indexer[2,1].n_local == 1
    @test pia_chunks[2].indexer[2,1].n_group1 == 1
    @test pia_chunks[2].indexer[2,1].start1 == 4
    @test pia_chunks[2].indexer[2,1].end1 == 4

    # # 1,1,2; 1,1,1,2
    # test that indexing points to the new groups of particles
    # new particles in cell 1 came from chunk 2
    # first the particles due to swapping
    @test chunk_exchanger.indexer[2,1].n_group1 == 1
    @test chunk_exchanger.indexer[2,1].start1 == 3
    @test chunk_exchanger.indexer[2,1].end1 == 3
    # then the extra pushed particles
    @test chunk_exchanger.indexer[2,1].n_group2 == 2
    @test chunk_exchanger.indexer[2,1].start2 == 4
    @test chunk_exchanger.indexer[2,1].end2 == 5

    # new particles in cell 2 came from chunk 1
    @test chunk_exchanger.indexer[1,2].n_group1 == 1
    @test chunk_exchanger.indexer[1,2].start1 == 1
    @test chunk_exchanger.indexer[1,2].end1 == 1
    # but no pushed particles
    @test chunk_exchanger.indexer[1,2].n_group2 == 0
    @test chunk_exchanger.indexer[1,2].start2 == 0
    @test chunk_exchanger.indexer[1,2].end2 == -1

    # particles don't come from the chunk themselves, obviously
    # and in this simple case chunk_id = cell_id
    for i in 1:2
        @test chunk_exchanger.indexer[i,i].n_group1 == 0
        @test chunk_exchanger.indexer[i,i].start1 == 0
        @test chunk_exchanger.indexer[i,i].end1 == -1
        @test chunk_exchanger.indexer[i,i].n_group2 == 0
        @test chunk_exchanger.indexer[i,i].start2 == 0
        @test chunk_exchanger.indexer[i,i].end2 == -1
    end

    # test case where all particles from chunk 2 are moved into chunk 1
    # 1,1; 1,1,1,1
    particles_chunks = [[ParticleVector(2)], [ParticleVector(4)]]
    pia_chunks = [ParticleIndexerArray(2,1), ParticleIndexerArray(2,1)]
    np_actual = [2,4]
    positions = [[0.0, 1.0], [-2.0, -1.0, 1.5, 0.5]]

    for chunk_id in 1:2
        for np in 1:np_actual[chunk_id]
            Merzbild.add_particle!(particles_chunks[chunk_id][1], np,
                                   chunk_id * 1.0, [chunk_id^3 * np, -np + 1.0, np],
                                   [positions[chunk_id][np], 0.5, 1.0+chunk_id])
        end
        pia_chunks[chunk_id].n_total[1] = np_actual[chunk_id]
        pia_chunks[chunk_id].indexer[1,1].n_local = np_actual[chunk_id]
        pia_chunks[chunk_id].indexer[1,1].n_group1 = np_actual[chunk_id]
        pia_chunks[chunk_id].indexer[1,1].start1 = 1
        pia_chunks[chunk_id].indexer[1,1].end1 = np_actual[chunk_id]
    end

    for j in 1:2
        reset!(chunk_exchanger, j)
    end
    exchange_particles!(chunk_exchanger, particles_chunks, pia_chunks, cell_chunks, 1)

    # we increase size by # of particles to add + DELTA_PARTICLES
    new_lengths = [6 + Merzbild.DELTA_PARTICLES, 4]
    np_actual = [6,4]
    new_positions = [[0.0, 1.0, -2.0, -1.0, 1.5, 0.5], [-2.0, -1.0, 1.5, 0.5]]
    new_weights = [[1.0, 1.0, 2.0, 2.0, 2.0, 2.0], [2.0, 2.0, 2.0, 2.0]]

    for chunk_id in 1:2
        @test pia_chunks[chunk_id].n_total[1] == np_actual[chunk_id]
    end

    @test particles_chunks[1][1].nbuffer == Merzbild.DELTA_PARTICLES
    @test particles_chunks[2][1].nbuffer == 4
    @test particles_chunks[2][1].buffer[1:4] == [1,2,3,4]
    for chunk_id in 1:2
        @test length(particles_chunks[chunk_id][1]) == new_lengths[chunk_id]
        for np in 1:np_actual[chunk_id]
            @test particles_chunks[chunk_id][1][np].x[1] == new_positions[chunk_id][np]
            @test particles_chunks[chunk_id][1][np].w == new_weights[chunk_id][np]
        end
    end

    # test part that is valid for both chunks
    np_actual_init = [2,4]

    for chunk_id in 1:2
        @test pia_chunks[chunk_id].indexer[2,1].n_local == 0
        @test pia_chunks[chunk_id].indexer[2,1].n_group1 == 0
        @test pia_chunks[chunk_id].indexer[2,1].start1 == 0
        @test pia_chunks[chunk_id].indexer[2,1].end1 == -1
        for cell in 1:2
            @test pia_chunks[chunk_id].indexer[cell,1].n_group2 == 0
            @test pia_chunks[chunk_id].indexer[cell,1].start2 == 0
            @test pia_chunks[chunk_id].indexer[cell,1].end2 == -1
        end
    end

    @test pia_chunks[1].indexer[1,1].n_local == np_actual_init[1]
    @test pia_chunks[1].indexer[1,1].n_group1 == np_actual_init[1]
    @test pia_chunks[1].indexer[1,1].start1 == 1
    @test pia_chunks[1].indexer[1,1].end1 == np_actual_init[1]

    @test pia_chunks[2].indexer[1,1].n_local == 0
    @test pia_chunks[2].indexer[1,1].n_group1 == 0
    @test pia_chunks[2].indexer[1,1].start1 == 0
    @test pia_chunks[2].indexer[1,1].end1 == -1

    # test chunk_exchanger
    # no new particles in cell 1 from chunk 2 due to swaps
    @test chunk_exchanger.indexer[2,1].n_group1 == 0
    @test chunk_exchanger.indexer[2,1].start1 == 0
    @test chunk_exchanger.indexer[2,1].end1 == -1
    # all particles due to push
    @test chunk_exchanger.indexer[2,1].n_group2 == 4
    @test chunk_exchanger.indexer[2,1].start2 == 3
    @test chunk_exchanger.indexer[2,1].end2 == 6

    # no new particles in cell 2
    @test chunk_exchanger.indexer[1,2].n_group1 == 0
    @test chunk_exchanger.indexer[1,2].start1 == 0
    @test chunk_exchanger.indexer[1,2].end1 == -1
    # but no pushed particles
    @test chunk_exchanger.indexer[1,2].n_group2 == 0
    @test chunk_exchanger.indexer[1,2].start2 == 0
    @test chunk_exchanger.indexer[1,2].end2 == -1

    # particles don't come from the chunk themselves, obviously
    # and in this simple case chunk_id = cell_id
    for i in 1:2
        @test chunk_exchanger.indexer[i,i].n_group1 == 0
        @test chunk_exchanger.indexer[i,i].start1 == 0
        @test chunk_exchanger.indexer[i,i].end1 == -1
        @test chunk_exchanger.indexer[i,i].n_group2 == 0
        @test chunk_exchanger.indexer[i,i].start2 == 0
        @test chunk_exchanger.indexer[i,i].end2 == -1
    end

    # 3 cells, 2 chunks
    # [1,1,2,3,3] [empty] (but pre-allocated) - all particles moved, buffer
    cell_chunks = [[1,2], [3]]
    n_cells = 3

    particles_chunks = [[ParticleVector(5)], [ParticleVector(4)]]
    pia_chunks = [ParticleIndexerArray(3,1), ParticleIndexerArray(3,1)]
    np_actual = [5, 0]
    positions = [[0.0, 1.0, 2.0, 6.0, 6.5], []]
    cells = [[1,1,2,3,3],[]]
    np_in_cells = [[2,1,2], [0,0,0]]  # per-chunk info
    offsets = [[1,3,4], [0,0,0]]  # this is for easier setup of pia

    chunk_exchanger = ChunkExchanger(cell_chunks, n_cells)
    for chunk_id in 1:2
        for np in 1:np_actual[chunk_id]
            Merzbild.add_particle!(particles_chunks[chunk_id][1], np,
                                   cells[chunk_id][np] * 1.0, [chunk_id, -chunk_id, chunk_id],
                                   [positions[chunk_id][np], 0.5, 0.0])
        end
        pia_chunks[chunk_id].n_total[1] = np_actual[chunk_id]

        for cell in 1:n_cells
            pia_chunks[chunk_id].indexer[cell,1].n_local = np_in_cells[chunk_id][cell]
            pia_chunks[chunk_id].indexer[cell,1].n_group1 = np_in_cells[chunk_id][cell]
            pia_chunks[chunk_id].indexer[cell,1].start1 = offsets[chunk_id][cell]
            pia_chunks[chunk_id].indexer[cell,1].end1 = offsets[chunk_id][cell] + np_in_cells[chunk_id][cell] - 1
        end
    end

    exchange_particles!(chunk_exchanger, particles_chunks, pia_chunks, cell_chunks, 1)

    # [1,1,2] [3,3]
    new_lengths = [5, 4]
    np_actual = [5,2]  # this includes particles pushed to another chunk!
    new_positions = [[0.0, 1.0, 2.0, 6.0, 6.5], [6.0, 6.5]]
    new_weights = [[1.0, 1.0, 2.0, 3.0, 3.0], [3.0, 3.0]]
    new_velocities = [[1.0, 1.0, 1.0, 1.0, 1.0], [1.0, 1.0]]

    for chunk_id in 1:2
        @test pia_chunks[chunk_id].n_total[1] == np_actual[chunk_id]
    end

    @test particles_chunks[1][1].nbuffer == 2
    @test particles_chunks[1][1].buffer[1:2] == [4,5]
    @test particles_chunks[2][1].nbuffer == 2

    for chunk_id in 1:2
        @test length(particles_chunks[chunk_id][1]) == new_lengths[chunk_id]
        for np in 1:np_actual[chunk_id]
            @test particles_chunks[chunk_id][1][np].x[1] == new_positions[chunk_id][np]
            @test particles_chunks[chunk_id][1][np].w == new_weights[chunk_id][np]
            @test particles_chunks[chunk_id][1][np].v[1] == new_velocities[chunk_id][np]
        end
    end

    # test that pia in chunk 1 does not point to particles sent to chunk 2
    np_in_cells_new = [[2,1,0], [0,0,0]]
    for chunk_id in 1:2
        for cell in 1:n_cells
            @test pia_chunks[chunk_id].indexer[cell,1].n_local == np_in_cells_new[chunk_id][cell]
            @test pia_chunks[chunk_id].indexer[cell,1].n_group1 == np_in_cells_new[chunk_id][cell]
            if np_in_cells_new[chunk_id][cell] > 0
                @test pia_chunks[chunk_id].indexer[cell,1].start1 == offsets[chunk_id][cell]
                @test pia_chunks[chunk_id].indexer[cell,1].end1 == offsets[chunk_id][cell] + np_in_cells_new[chunk_id][cell] - 1
            else
                @test pia_chunks[chunk_id].indexer[cell,1].start1 == 0
                @test pia_chunks[chunk_id].indexer[cell,1].end1 == -1
            end
            @test pia_chunks[chunk_id].indexer[cell,1].n_group2 == 0
            @test pia_chunks[chunk_id].indexer[cell,1].start2 == 0
            @test pia_chunks[chunk_id].indexer[cell,1].end2 == -1
        end
    end

    # no new particles in cells 1, 2 in both chunks
    for chunk_id in 1:2
        for cell in 1:2
            @test chunk_exchanger.indexer[chunk_id,cell].n_group1 == 0
            @test chunk_exchanger.indexer[chunk_id,cell].start1 == 0
            @test chunk_exchanger.indexer[chunk_id,cell].end1 == -1

            @test chunk_exchanger.indexer[chunk_id,cell].n_group2 == 0
            @test chunk_exchanger.indexer[chunk_id,cell].start2 == 0
            @test chunk_exchanger.indexer[chunk_id,cell].end2 == -1
        end
    end

    # cell 3 belongs to chunk 2, therefore nothing 
    @test chunk_exchanger.indexer[2,3].n_group1 == 0
    @test chunk_exchanger.indexer[2,3].start1 == 0
    @test chunk_exchanger.indexer[2,3].end1 == -1

    @test chunk_exchanger.indexer[2,3].n_group2 == 0
    @test chunk_exchanger.indexer[2,3].start2 == 0
    @test chunk_exchanger.indexer[2,3].end2 == -1

    # no new particles in cell 3 due to swap
    @test chunk_exchanger.indexer[1,3].n_group1 == 0
    @test chunk_exchanger.indexer[1,3].start1 == 0
    @test chunk_exchanger.indexer[1,3].end1 == -1
    # and everything due to push
    @test chunk_exchanger.indexer[1,3].n_group2 == 2
    @test chunk_exchanger.indexer[1,3].start2 == 1
    @test chunk_exchanger.indexer[1,3].end2 == 2


    # 4 cells, 3 chunks
    # [1], [2,3], [4]
    # [1,1,2,3,4,4], [1,3,4], [1,2,3,3]
    cell_chunks = [[1], [2,3], [4]]
    n_cells = 4

    particles_chunks = [[ParticleVector(8)], [ParticleVector(4)], [ParticleVector(5)]]
    pia_chunks = [ParticleIndexerArray(n_cells,1), ParticleIndexerArray(n_cells,1), ParticleIndexerArray(n_cells,1)]
    np_actual = [6, 3, 4]
    positions = [[0.0, 0.5, 2.0, 3.0, 4.5, 5.0], [1.5, 3.5, 5.5], [-1.0, 2.5, 4.0, 4.5]]
    cells = [[1,1,2,3,4,4],[1,3,4],[1,2,3,3]]
    np_in_cells = [[2,1,1,2], [1,0,1,1], [1,1,2,0]]  # per-chunk info
    offsets = [[1,3,4,5], [1,1,2,3], [1,2,3,3]]  # this is for easier setup of pia

    chunk_exchanger = ChunkExchanger(cell_chunks, n_cells)
    for chunk_id in 1:3
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

    # denote specific particles
    # [1_1,1_2,2_1,3_1,4_1,4_2], [1_3,3_2,4_3], [1_4,2_2,3_3,3_4]
    # becomes [1_1,1_2,1_3,3_1,4_1,4_2], [2_1,3_2,4_3,3_1], [1_4,2_2,3_3,3_4] after one pass (chunks 1 & 2)
    # then chunks 1 & 3
    # [1_1,1_2,1_3,3_1,1_4,4_2], [2_1,3_2,4_3,3_1], [4_1,2_2,3_3,3_4,4_2]
    # then chunks 2 & 3
    # [1_1,1_2,1_3,3_1,1_4,4_2], [2_1,3_2,2_2,3_1,3_3,3_4], [4_1,4_3,3_3,3_4,4_2]
    new_lengths = [8, 6 + Merzbild.DELTA_PARTICLES, 5]
    np_actual = [6,6,5]  # this includes particles pushed to another chunk!

    for chunk_id in 1:3
        @test pia_chunks[chunk_id].n_total[1] == np_actual[chunk_id]
    end

    # our initial particle array was of length 8, so buffer was of length 2 before
    # we pushed particles
    @test particles_chunks[1][1].nbuffer == 4
    @test particles_chunks[1][1].buffer[1:4] == [8,7,4,6]

    # we increased size of array
    @test particles_chunks[2][1].nbuffer == Merzbild.DELTA_PARTICLES
    # @test particles_chunks[2][1].buffer[1] == 4

    # had pv of length 4 and 3 particles, so buffer was of size 1 before we pushed
    @test particles_chunks[3][1].nbuffer == 2
    @test particles_chunks[3][1].buffer[1:2] == [3,4]

    # initial positions:
    # [[0.0, 0.5, 2.0, 3.0, 4.5, 5.0], [1.5, 3.5, 5.5], [-1.0, 2.5, 4.0, 4.5]]
    # [[1_1, 1_2, 2_1, 3_1, 4_1, 4_2], [1_3, 3_2, 4_3], [1_4,  2_2, 3_3, 3_4]]
    # new particle ordering:
    # [1_1,1_2,1_3,3_1,1_4,4_2], [2_1,3_2,2_2,3_1,3_3,3_4], [4_1,4_3,3_3,3_4,4_2]
    new_positions = [[0.0, 0.5, 1.5, 3.0, -1.0, 5.0],
                     [2.0, 3.5, 2.5, 3.0, 4.0, 4.5],
                     [4.5, 5.5, 4.0, 4.5, 5.0]]
    new_weights = [[1.0, 1.0, 1.0, 3.0, 1.0, 4.0], [2.0, 3.0, 2.0, 3.0, 3.0, 3.0], [4.0, 4.0, 3.0, 3.0, 4.0]]
    new_velocities = [[1.0, 1.0, 2.0, 1.0, 3.0, 1.0],
                      [1.0, 2.0, 3.0, 1.0, 3.0, 3.0],
                      [1.0, 2.0, 3.0, 3.0, 1.0]]

    for chunk_id in 1:3
        @test length(particles_chunks[chunk_id][1]) == new_lengths[chunk_id]
        for np in 1:np_actual[chunk_id]
            @test particles_chunks[chunk_id][1][np].x[1] == new_positions[chunk_id][np]
            @test particles_chunks[chunk_id][1][np].w == new_weights[chunk_id][np]
            @test particles_chunks[chunk_id][1][np].v[1] == new_velocities[chunk_id][np]
        end
    end

    # these are the actual particle counts if we disregard
    # the particles that have been pushed to another chunk
    # so basically exchange_particles! can only lead to a decrease
    # in the indexing values in pia, except for n_total
    np_in_cells_new = [[2,0,0,0], [0,0,1,0], [0,0,0,0]]
    for chunk_id in 1:3
        for cell in 1:n_cells
            @test pia_chunks[chunk_id].indexer[cell,1].n_local == np_in_cells_new[chunk_id][cell]
            @test pia_chunks[chunk_id].indexer[cell,1].n_group1 == np_in_cells_new[chunk_id][cell]
            @test pia_chunks[chunk_id].indexer[cell,1].n_group2 == 0
            @test pia_chunks[chunk_id].indexer[cell,1].start2 == 0
            @test pia_chunks[chunk_id].indexer[cell,1].end2 == -1

            if np_in_cells_new[chunk_id][cell] > 0
                @test pia_chunks[chunk_id].indexer[cell,1].start1 == offsets[chunk_id][cell]
                @test pia_chunks[chunk_id].indexer[cell,1].end1 == offsets[chunk_id][cell] + np_in_cells_new[chunk_id][cell] - 1
            end
        end
    end

    # cell 1
    @test chunk_exchanger.indexer[1,1].n_group1 == 0
    @test chunk_exchanger.indexer[1,1].start1 == 0
    @test chunk_exchanger.indexer[1,1].end1 == -1

    @test chunk_exchanger.indexer[2,1].n_group1 == 1
    @test chunk_exchanger.indexer[2,1].start1 == 3
    @test chunk_exchanger.indexer[2,1].end1 == 3

    @test chunk_exchanger.indexer[3,1].n_group1 == 1
    @test chunk_exchanger.indexer[3,1].start1 == 5
    @test chunk_exchanger.indexer[3,1].end1 == 5
    
    # cell 1, no pushes, only swaps, therefore group2 is empty
    for chunk_id in 1:3
        @test chunk_exchanger.indexer[chunk_id,1].n_group2 == 0
        @test chunk_exchanger.indexer[chunk_id,1].start2 == 0
        @test chunk_exchanger.indexer[chunk_id,1].end2 == -1
    end
    
    # cell 2
    @test chunk_exchanger.indexer[1,2].n_group1 == 1
    @test chunk_exchanger.indexer[1,2].start1 == 1
    @test chunk_exchanger.indexer[1,2].end1 == 1

    @test chunk_exchanger.indexer[1,2].n_group2 == 0
    @test chunk_exchanger.indexer[1,2].start2 == 0
    @test chunk_exchanger.indexer[1,2].end2 == -1

    # no self-pushes/swaps
    for cell in [2,3]
        @test chunk_exchanger.indexer[2,cell].n_group1 == 0
        @test chunk_exchanger.indexer[2,cell].start1 == 0
        @test chunk_exchanger.indexer[2,cell].end1 == -1

        @test chunk_exchanger.indexer[2,cell].n_group2 == 0
        @test chunk_exchanger.indexer[2,cell].start2 == 0
        @test chunk_exchanger.indexer[2,cell].end2 == -1
    end

    @test chunk_exchanger.indexer[3,2].n_group1 == 1
    @test chunk_exchanger.indexer[3,2].start1 == 3
    @test chunk_exchanger.indexer[3,2].end1 == 3

    @test chunk_exchanger.indexer[3,2].n_group2 == 0
    @test chunk_exchanger.indexer[3,2].start2 == 0
    @test chunk_exchanger.indexer[3,2].end2 == -1

    # cell 3
    @test chunk_exchanger.indexer[1,3].n_group1 == 0
    @test chunk_exchanger.indexer[1,3].start1 == 0
    @test chunk_exchanger.indexer[1,3].end1 == -1

    @test chunk_exchanger.indexer[1,3].n_group2 == 1
    @test chunk_exchanger.indexer[1,3].start2 == 4
    @test chunk_exchanger.indexer[1,3].end2 == 4

    # from chunk 3
    @test chunk_exchanger.indexer[3,3].n_group1 == 0
    @test chunk_exchanger.indexer[3,3].start1 == 0
    @test chunk_exchanger.indexer[3,3].end1 == -1

    @test chunk_exchanger.indexer[3,3].n_group2 == 2
    @test chunk_exchanger.indexer[3,3].start2 == 5
    @test chunk_exchanger.indexer[3,3].end2 == 6

    # cell 4
    @test chunk_exchanger.indexer[1,4].n_group1 == 1
    @test chunk_exchanger.indexer[1,4].start1 == 1
    @test chunk_exchanger.indexer[1,4].end1 == 1

    @test chunk_exchanger.indexer[1,4].n_group2 == 1
    @test chunk_exchanger.indexer[1,4].start2 == 5
    @test chunk_exchanger.indexer[1,4].end2 == 5

    @test chunk_exchanger.indexer[2,4].n_group1 == 1
    @test chunk_exchanger.indexer[2,4].start1 == 2
    @test chunk_exchanger.indexer[2,4].end1 == 2

    @test chunk_exchanger.indexer[2,4].n_group2 == 0
    @test chunk_exchanger.indexer[2,4].start2 == 0
    @test chunk_exchanger.indexer[2,4].end2 == -1

    @test chunk_exchanger.indexer[3,4].n_group1 == 0
    @test chunk_exchanger.indexer[3,4].start1 == 0
    @test chunk_exchanger.indexer[3,4].end1 == -1

    @test chunk_exchanger.indexer[3,4].n_group2 == 0
    @test chunk_exchanger.indexer[3,4].start2 == 0
    @test chunk_exchanger.indexer[3,4].end2 == -1

    # number of particles in a cell
    # [1,1,2,3,4,4], [1,3,4], [1,2,3,3] - what we had
    np_correct = [4, 2, 4, 3]
    # weight is the same as the cell the particle belongs to
    w_correct = [4.0, 2.0*2, 3.0*4, 4.0*3]

    # which chunk a cell belongs to
    chunk_id_map_cell = [1, 2, 2, 3]
    for cell in 1:4
        np = 0
        w = 0.0

        chunk_id = chunk_id_map_cell[cell]
        pia = pia_chunks[chunk_id]
        for i in pia.indexer[cell,1].start1:pia.indexer[cell,1].end1
            np += 1
            w += particles_chunks[chunk_id][1][i].w
        end

        for chunk_id_2 in 1:3
            # we iterate over all chunks
            s1 = chunk_exchanger.indexer[chunk_id_2,cell].start1
            e1 = chunk_exchanger.indexer[chunk_id_2,cell].end1
            s2 = chunk_exchanger.indexer[chunk_id_2,cell].start2
            e2 = chunk_exchanger.indexer[chunk_id_2,cell].end2
            for i in s1:e1
                np += 1
                w += particles_chunks[chunk_id][1][i].w
            end
            for i in s2:e2
                np += 1
                w += particles_chunks[chunk_id][1][i].w
            end
        end

        @test np == np_correct[cell]
        @test w == w_correct[cell]
    end

    # 4 cells, 2 chunks
    # [2,3,4], [1,1,1,1]
    # becomes [1,1,1,1], [2,3,4]
    cell_chunks = [[1], [2,3,4]]
    n_cells = 4

    particles_chunks = [[ParticleVector(3)], [ParticleVector(4)]]
    pia_chunks = [ParticleIndexerArray(n_cells,1), ParticleIndexerArray(n_cells,1)]
    np_actual = [3, 4]
    positions = [[2.0, 3.0, 4.0], [-1.0, -0.5, 0.5, 1.0]]
    cells = [[2,3,4],[1,1,1,1]]
    np_in_cells = [[0,1,1,1], [4,0,0,0]]  # per-chunk info
    offsets = [[0,1,2,3], [1,0,0,0]]  # per-cell, this is for easier setup of pia

    chunk_exchanger = ChunkExchanger(cell_chunks, n_cells)
    for chunk_id in 1:2
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

    new_lengths = [4 + Merzbild.DELTA_PARTICLES, 4]
    np_actual = [4,4]  # this includes particles pushed to another chunk!

    for chunk_id in 1:2
        @test pia_chunks[chunk_id].n_total[1] == np_actual[chunk_id]
    end

    new_positions = [[-1.0, -0.5, 0.5, 1.0],
                     [2.0, 3.0, 4.0, 1.0]]
    new_weights = [[1.0, 1.0, 1.0, 1.0], [2.0, 3.0, 4.0, 1.0]]
    new_velocities = [[2.0, 2.0, 2.0, 2.0],
                      [1.0, 1.0, 1.0, 2.0]]

    for chunk_id in 1:2
        @test length(particles_chunks[chunk_id][1]) == new_lengths[chunk_id]
        for np in 1:np_actual[chunk_id]
            @test particles_chunks[chunk_id][1][np].x[1] == new_positions[chunk_id][np]
            @test particles_chunks[chunk_id][1][np].w == new_weights[chunk_id][np]
            @test particles_chunks[chunk_id][1][np].v[1] == new_velocities[chunk_id][np]
        end
    end

    # cell 1
    # nothing from self
    @test chunk_exchanger.indexer[1,1].n_group1 == 0
    @test chunk_exchanger.indexer[1,1].start1 == 0
    @test chunk_exchanger.indexer[1,1].end1 == -1

    @test chunk_exchanger.indexer[1,1].n_group2 == 0
    @test chunk_exchanger.indexer[1,1].start2 == 0
    @test chunk_exchanger.indexer[1,1].end2 == -1

    # 3 via swaps
    @test chunk_exchanger.indexer[2,1].n_group1 == 3
    @test chunk_exchanger.indexer[2,1].start1 == 1
    @test chunk_exchanger.indexer[2,1].end1 == 3

    # 1 via push
    @test chunk_exchanger.indexer[2,1].n_group2 == 1
    @test chunk_exchanger.indexer[2,1].start2 == 4
    @test chunk_exchanger.indexer[2,1].end2 == 4


    # cells 2,3,4
    for cell in [2,3,4]
        # nothing from self
        @test chunk_exchanger.indexer[2,cell].n_group1 == 0
        @test chunk_exchanger.indexer[2,cell].start1 == 0
        @test chunk_exchanger.indexer[2,cell].end1 == -1

        @test chunk_exchanger.indexer[2,cell].n_group2 == 0
        @test chunk_exchanger.indexer[2,cell].start2 == 0
        @test chunk_exchanger.indexer[2,cell].end2 == -1

        # and nothing pushed from chunk 1
        @test chunk_exchanger.indexer[1,cell].n_group2 == 0
        @test chunk_exchanger.indexer[1,cell].start2 == 0
        @test chunk_exchanger.indexer[1,cell].end2 == -1
    end

    # test swaps from chunk 1
    @test chunk_exchanger.indexer[1,2].n_group1 == 1
    @test chunk_exchanger.indexer[1,2].start1 == 1
    @test chunk_exchanger.indexer[1,2].end1 == 1

    @test chunk_exchanger.indexer[1,3].n_group1 == 1
    @test chunk_exchanger.indexer[1,3].start1 == 2
    @test chunk_exchanger.indexer[1,3].end1 == 2

    @test chunk_exchanger.indexer[1,4].n_group1 == 1
    @test chunk_exchanger.indexer[1,4].start1 == 3
    @test chunk_exchanger.indexer[1,4].end1 == 3
end