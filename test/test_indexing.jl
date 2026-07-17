@testset "test_indexing" begin
    pia = ParticleIndexerArray(20)

    @test pia.n_cells == 1
    @test pia.n_species == 1
    @test pia.indexer[1,1].n_local == 20
    @test pia.indexer[1,1].start1 == 1
    @test pia.indexer[1,1].end1 == 20
    @test pia.indexer[1,1].n_group1 == 20
    @test pia.indexer[1,1].n_group2 == 0

    Merzbild.update_particle_indexer_new_particle!(pia, 1, 1)
    Merzbild.update_particle_indexer_new_particle!(pia, 1, 1)

    @test pia.n_cells == 1
    @test pia.n_species == 1
    @test pia.n_total[1] == 22
    @test pia.indexer[1,1].n_local == 22
    @test pia.indexer[1,1].start1 == 1
    @test pia.indexer[1,1].end1 == 20
    @test pia.indexer[1,1].n_group1 == 20
    @test pia.indexer[1,1].start2 == 21
    @test pia.indexer[1,1].end2 == 22
    @test pia.indexer[1,1].n_group2 == 2

    pia = ParticleIndexerArray(10)

    pia.n_total[1] += 5
    pia.indexer[1,1].n_local += 5
    pia.indexer[1,1].start2 = 31
    pia.indexer[1,1].end2 = 35
    pia.indexer[1,1].n_group2 = 5

    i0_arr = [1, 4, 10]
    for i0 in i0_arr
        index = Merzbild.map_cont_index(pia, 1, 1, i0-1)  # maps from 0, n_local - 1 to actual particle
        @test index == 1 + i0 - 1
    end

    i0_arr = [11, 12, 15]
    for i0 in i0_arr
        index = Merzbild.map_cont_index(pia, 1, 1, i0-1)
        @test index == 31 + i0 - 10 - 1
    end

    pia.n_total[1] = 10
    pia.indexer[1,1].n_local = 10
    pia.indexer[1,1].start1 = 24
    pia.indexer[1,1].end1 = 28
    pia.indexer[1,1].n_group1 = 5
    pia.indexer[1,1].start2 = 31
    pia.indexer[1,1].end2 = 35
    pia.indexer[1,1].n_group2 = 5

    i0_arr = [1, 4, 5]
    for i0 in i0_arr
        index = Merzbild.map_cont_index(pia, 1, 1, i0-1)  # maps from 0, n_local - 1 to actual particle
        @test index == 24 + i0 - 1
    end

    i0_arr = [6, 8, 10]
    for i0 in i0_arr
        index = Merzbild.map_cont_index(pia, 1, 1, i0-1)
        @test index == 31 + i0 - 5 - 1
    end

    pia = ParticleIndexerArray(10, 2)
    @test pia.n_cells == 10
    @test pia.n_species == 2
    @test length(pia.n_total) == 2
    @test size(pia.indexer) == (10, 2)
end