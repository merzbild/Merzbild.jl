@testset "test_indexing" begin
    particle_indexer = create_particle_indexer_array(20)

    @test particle_indexer.n_total[1] == 20
    @test particle_indexer.indexer[1,1].n_local == 20
    @test particle_indexer.indexer[1,1].start1 == 1
    @test particle_indexer.indexer[1,1].end1 == 20
    @test particle_indexer.indexer[1,1].n_group1 == 20
    @test particle_indexer.indexer[1,1].n_group2 == 0

    Merzbild.update_particle_indexer_new_particle(1, 1, particle_indexer)
    Merzbild.update_particle_indexer_new_particle(1, 1, particle_indexer)

    @test particle_indexer.n_total[1] == 22
    @test particle_indexer.indexer[1,1].n_local == 22
    @test particle_indexer.indexer[1,1].start1 == 1
    @test particle_indexer.indexer[1,1].end1 == 20
    @test particle_indexer.indexer[1,1].n_group1 == 20
    @test particle_indexer.indexer[1,1].start2 == 21
    @test particle_indexer.indexer[1,1].end2 == 22
    @test particle_indexer.indexer[1,1].n_group2 == 2

    particle_indexer = create_particle_indexer_array(10)

    particle_indexer.n_total[1] += 5
    particle_indexer.indexer[1,1].n_local += 5
    particle_indexer.indexer[1,1].start2 = 31
    particle_indexer.indexer[1,1].end2 = 35
    particle_indexer.indexer[1,1].n_group2 = 5

    i0_arr = [1, 4, 10]
    for i0 in i0_arr
        index = Merzbild.map_cont_index(1, 1, particle_indexer, i0-1)  # maps from 0, n_local - 1 to actual particle
        @test index == 1 + i0 - 1
    end

    i0_arr = [11, 12, 15]
    for i0 in i0_arr
        index = Merzbild.map_cont_index(1, 1, particle_indexer, i0-1)
        @test index == 31 + i0 - 10 - 1
    end


    particle_indexer.n_total[1] = 10
    particle_indexer.indexer[1,1].n_local = 10
    particle_indexer.indexer[1,1].start1 = 24
    particle_indexer.indexer[1,1].end1 = 28
    particle_indexer.indexer[1,1].n_group1 = 5
    particle_indexer.indexer[1,1].start2 = 31
    particle_indexer.indexer[1,1].end2 = 35
    particle_indexer.indexer[1,1].n_group2 = 5

    i0_arr = [1, 4, 5]
    for i0 in i0_arr
        index = Merzbild.map_cont_index(1, 1, particle_indexer, i0-1)  # maps from 0, n_local - 1 to actual particle
        @test index == 24 + i0 - 1
    end

    i0_arr = [6, 8, 10]
    for i0 in i0_arr
        index = Merzbild.map_cont_index(1, 1, particle_indexer, i0-1)
        @test index == 31 + i0 - 5 - 1
    end

    Merzbild.update_particle_indexer_new_lower_count(1, 1, particle_indexer, 8)
    @test particle_indexer.n_total[1] == 8
    @test particle_indexer.indexer[1,1].n_local == 8
    @test particle_indexer.indexer[1,1].start1 == 24
    @test particle_indexer.indexer[1,1].end1 == 28
    @test particle_indexer.indexer[1,1].n_group1 == 5
    @test particle_indexer.indexer[1,1].start2 == 31
    @test particle_indexer.indexer[1,1].end2 == 33
    @test particle_indexer.indexer[1,1].n_group2 == 3


    Merzbild.update_particle_indexer_new_lower_count(1, 1, particle_indexer, 5)
    @test particle_indexer.n_total[1] == 5
    @test particle_indexer.indexer[1,1].n_local == 5
    @test particle_indexer.indexer[1,1].start1 == 24
    @test particle_indexer.indexer[1,1].end1 == 28
    @test particle_indexer.indexer[1,1].n_group1 == 5
    @test particle_indexer.indexer[1,1].start2 == 0
    @test particle_indexer.indexer[1,1].end2 == 0
    @test particle_indexer.indexer[1,1].n_group2 == 0


    particle_indexer.n_total[1] = 10
    particle_indexer.indexer[1,1].n_local = 10
    particle_indexer.indexer[1,1].start1 = 24
    particle_indexer.indexer[1,1].end1 = 28
    particle_indexer.indexer[1,1].n_group1 = 5
    particle_indexer.indexer[1,1].start2 = 31
    particle_indexer.indexer[1,1].end2 = 35
    particle_indexer.indexer[1,1].n_group2 = 5

    Merzbild.update_particle_indexer_new_lower_count(1, 1, particle_indexer, 2)
    @test particle_indexer.n_total[1] == 2
    @test particle_indexer.indexer[1,1].n_local == 2
    @test particle_indexer.indexer[1,1].start1 == 24
    @test particle_indexer.indexer[1,1].end1 == 25
    @test particle_indexer.indexer[1,1].n_group1 == 2
    @test particle_indexer.indexer[1,1].start2 == 0
    @test particle_indexer.indexer[1,1].end2 == 0
    @test particle_indexer.indexer[1,1].n_group2 == 0
end