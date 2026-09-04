@testset "last index and particle deletion/addition" begin

    function generate_pia_and_particles(n_particles, n_group1, n_group2)
        particles = ParticleVector{0}(n_particles)
        pia = ParticleIndexerArray(length(n_group1), 1)

        start = 1
        offset = 0
        for (cell, ng1) in enumerate(n_group1)
            if ng1 > 0
                pia.indexer[cell, 1].start1 = start
                pia.indexer[cell, 1].end1 = start + ng1 - 1
                pia.indexer[cell, 1].n_group1 = ng1

                start += ng1

                for _ in 1:ng1
                    offset += 1
                    Merzbild.add_particle!(particles, offset, cell * 1.0,
                                           SVector{3, Float64}(1.0, 1.0, 1.0),
                                           SVector{0, Float64}())
                end
            end
        end

        for (cell, ng2) in enumerate(n_group2)
            if ng2 > 0
                pia.indexer[cell, 1].start2 = start
                pia.indexer[cell, 1].end2 = start + ng2 - 1
                pia.indexer[cell, 1].n_group2 = ng2

                start += ng2

                for _ in 1:ng2
                    offset += 1
                    Merzbild.add_particle!(particles, offset, cell * 1.0,
                                           SVector{3, Float64}(2.0, 2.0, 2.0),
                                           SVector{0, Float64}())
                end
            end
        end

        n_total = 0
        for (cell, ng1) in enumerate(n_group1)
            pia.indexer[cell,1].n_local = pia.indexer[cell,1].n_group1 + pia.indexer[cell,1].n_group2
            n_total += pia.indexer[cell,1].n_local
        end

        pia.n_total[1] = n_total
        pia.index_last[1] = n_total
        pia.contiguous[1] = true

        return pia, particles
    end

    pia, particles = generate_pia_and_particles(44, [6, 5, 5, 0, 5, 0, 3], [5, 0, 3, 3, 2, 0, 4])

    n_total = pia.n_total[1]
    last_index = pia.index_last[1]

    # start by checking that everything has been initialized correctly
    @test check_pia_is_correct(pia, 1) == (true, 0)
    @test check_unique_index(particles, pia, 1) == (true, 0)
    @test check_unique_buffer(particles) == (true, 0)
    @test last_index == n_total

    # delete only from group 1
    Merzbild.delete_particle!(particles, pia, 7, 1, 23)
    Merzbild.delete_particle!(particles, pia, 5, 1, 21)
    Merzbild.delete_particle!(particles, pia, 4, 1, 33)
    Merzbild.delete_particle!(particles, pia, 3, 1, 13)
    Merzbild.delete_particle!(particles, pia, 2, 1, 8)
    Merzbild.delete_particle!(particles, pia, 1, 1, 6)

    # everything still correct
    @test check_pia_is_correct(pia, 1) == (true, 0)
    @test check_unique_index(particles, pia, 1) == (true, 0)
    @test check_unique_buffer(particles) == (true, 0)
    @test last_index == pia.index_last[1]
    @test pia.index_last[1] > pia.n_total[1]
    @test pia.n_total[1] == 35
    @test pia.index_last[1] == 41

    # create particle in cell 7 in group 2
    # group2 was [38, 41] before
    Merzbild.update_buffer_index_new_particle!(particles, pia, 7, 1)

    @test pia.indexer[7,1].n_group2 == 5
    @test pia.indexer[7,1].start2 == 38
    @test pia.indexer[7,1].end2 == 42
    @test pia.index_last[1] == 42
    @test pia.n_total[1] == 36

    @test check_pia_is_correct(pia, 1) == (true, 0)
    @test check_unique_index(particles, pia, 1) == (true, 0)
    @test check_unique_buffer(particles) == (true, 0)

    # case
    pia, particles = generate_pia_and_particles(30, [4, 0, 0], [0, 2, 1])

    @test pia.index_last[1] == 7
    @test pia.n_total[1] == 7

    Merzbild.delete_particle_end_group2!(particles, pia, 3, 1)
    @test pia.index_last[1] == 6
    @test pia.n_total[1] == 6

    @test check_pia_is_correct(pia, 1) == (true, 0)
    @test check_unique_index(particles, pia, 1) == (true, 0)
    @test check_unique_buffer(particles) == (true, 0)

    # case
    pia, particles = generate_pia_and_particles(30, [4, 0, 0], [0, 2, 1])

    @test pia.index_last[1] == 7
    @test pia.n_total[1] == 7

    Merzbild.delete_particle_end_group2!(particles, pia, 2, 1)
    Merzbild.delete_particle_end_group2!(particles, pia, 2, 1)
    @test pia.index_last[1] == 7
    @test pia.n_total[1] == 5
    Merzbild.delete_particle_end_group2!(particles, pia, 3, 1)
    @test pia.index_last[1] == 4
    @test pia.n_total[1] == 4

    @test check_pia_is_correct(pia, 1) == (true, 0)
    @test check_unique_index(particles, pia, 1) == (true, 0)
    @test check_unique_buffer(particles) == (true, 0)

    # case
    pia, particles = generate_pia_and_particles(30, [4, 0, 0], [0, 2, 1])

    @test pia.index_last[1] == 7
    @test pia.n_total[1] == 7

    Merzbild.delete_particle_end_group1!(particles, pia, 1, 1)
    Merzbild.delete_particle_end_group1!(particles, pia, 1, 1)
    @test pia.index_last[1] == 7
    @test pia.n_total[1] == 5
    Merzbild.delete_particle_end_group2!(particles, pia, 2, 1)
    Merzbild.delete_particle_end_group2!(particles, pia, 2, 1)
    @test pia.index_last[1] == 7
    @test pia.n_total[1] == 3

    @test check_pia_is_correct(pia, 1) == (true, 0)
    @test check_unique_index(particles, pia, 1) == (true, 0)
    @test check_unique_buffer(particles) == (true, 0)

    # case
    pia, particles = generate_pia_and_particles(30, [3, 2, 1], [1, 2, 3])
    @test pia.index_last[1] == 12
    @test pia.n_total[1] == 12
    
    Merzbild.delete_particle_end!(particles, pia, 1, 1)

    # [3, 2, 1], [0, 2, 3]
    @test pia.index_last[1] == 12
    @test pia.n_total[1] == 11

    Merzbild.delete_particle_end!(particles, pia, 1, 1)
    Merzbild.delete_particle_end!(particles, pia, 1, 1)
    Merzbild.delete_particle_end!(particles, pia, 1, 1)

    # [0, 2, 1], [0, 2, 3]
    @test pia.index_last[1] == 12
    @test pia.n_total[1] == 8

    @test check_pia_is_correct(pia, 1) == (true, 0)
    @test check_unique_index(particles, pia, 1) == (true, 0)
    @test check_unique_buffer(particles) == (true, 0)

    Merzbild.delete_particle_end!(particles, pia, 2, 1)
    Merzbild.delete_particle_end!(particles, pia, 2, 1)
    Merzbild.delete_particle_end!(particles, pia, 2, 1)

    # [0, 1, 1], [0, 0, 3]
    @test pia.index_last[1] == 12
    @test pia.n_total[1] == 5

    @test check_pia_is_correct(pia, 1) == (true, 0)
    @test check_unique_index(particles, pia, 1) == (true, 0)
    @test check_unique_buffer(particles) == (true, 0)

    Merzbild.delete_particle_end!(particles, pia, 3, 1)
    Merzbild.delete_particle_end!(particles, pia, 3, 1)

    # [0, 1, 1], [0, 0, 1]
    @test pia.index_last[1] == 10
    @test pia.n_total[1] == 3

    @test check_pia_is_correct(pia, 1) == (true, 0)
    @test check_unique_index(particles, pia, 1) == (true, 0)
    @test check_unique_buffer(particles) == (true, 0)

    # [0, 1, 1], [0, 0, 0]
    Merzbild.delete_particle_end!(particles, pia, 3, 1)
    @test pia.index_last[1] == 6
    @test pia.n_total[1] == 2

    @test check_pia_is_correct(pia, 1) == (true, 0)
    @test check_unique_index(particles, pia, 1) == (true, 0)
    @test check_unique_buffer(particles) == (true, 0)

    # check setting of index last once we have to iterate over preceding cells
    pia, particles = generate_pia_and_particles(10, [7, 2, 1], [0, 0, 0])
    @test pia.index_last[1] == 10
    @test pia.n_total[1] == 10

    Merzbild.delete_particle_end!(particles, pia, 2, 1)
    @test pia.index_last[1] == 10
    @test pia.n_total[1] == 9

    Merzbild.delete_particle_end!(particles, pia, 2, 1)
    @test pia.index_last[1] == 10
    @test pia.n_total[1] == 8
    @test pia.indexer[2,1].n_local == 0

    Merzbild.delete_particle_end!(particles, pia, 3, 1)
    @test pia.index_last[1] == 7

    # now we test batch deletion
    pia, particles = generate_pia_and_particles(30, [4, 0, 3], [0, 2, 2])

    @test pia.index_last[1] == 11
    @test pia.n_total[1] == 11

    Merzbild.delete_batch_end!(particles, pia, 3, 1, 1)

    @test pia.index_last[1] == 10
    @test pia.n_total[1] == 10

    @test pia.indexer[1,1].n_local == 4
    @test pia.indexer[1,1].n_group1 == 4
    @test pia.indexer[1,1].start1 == 1
    @test pia.indexer[1,1].end1 == 4
    @test pia.indexer[1,1].n_group2 == 0
    @test pia.indexer[1,1].start2 == 0
    @test pia.indexer[1,1].end2 == -1

    @test pia.indexer[2,1].n_local == 2
    @test pia.indexer[2,1].n_group1 == 0
    @test pia.indexer[2,1].start1 == 0
    @test pia.indexer[2,1].end1 == -1
    @test pia.indexer[2,1].n_group2 == 2
    @test pia.indexer[2,1].start2 == 8
    @test pia.indexer[2,1].end2 == 9

    @test pia.indexer[3,1].n_local == 4
    @test pia.indexer[3,1].n_group1 == 3
    @test pia.indexer[3,1].n_group2 == 1
    @test pia.indexer[3,1].start1 == 5
    @test pia.indexer[3,1].end1 == 7
    @test pia.indexer[3,1].start2 == 10
    @test pia.indexer[3,1].end2 == 10

    Merzbild.delete_batch_end!(particles, pia, 2, 1, 3)

    @test pia.index_last[1] == 10
    @test pia.n_total[1] == 8

    @test pia.indexer[1,1].n_local == 4
    @test pia.indexer[1,1].n_group1 == 4
    @test pia.indexer[1,1].start1 == 1
    @test pia.indexer[1,1].end1 == 4
    @test pia.indexer[1,1].n_group2 == 0
    @test pia.indexer[1,1].start2 == 0
    @test pia.indexer[1,1].end2 == -1

    @test pia.indexer[2,1].n_local == 0
    @test pia.indexer[2,1].n_group1 == 0
    @test pia.indexer[2,1].start1 == 0
    @test pia.indexer[2,1].end1 == -1
    @test pia.indexer[2,1].n_group2 == 0
    @test pia.indexer[2,1].start2 == 0
    @test pia.indexer[2,1].end2 == -1

    @test pia.indexer[3,1].n_local == 4
    @test pia.indexer[3,1].n_group1 == 3
    @test pia.indexer[3,1].n_group2 == 1
    @test pia.indexer[3,1].start1 == 5
    @test pia.indexer[3,1].end1 == 7
    @test pia.indexer[3,1].start2 == 10
    @test pia.indexer[3,1].end2 == 10

    # current state is [4, 0, 3], [0, 0, 1]
    Merzbild.delete_batch_end!(particles, pia, 1, 1, 2)

    @test pia.index_last[1] == 10
    @test pia.n_total[1] == 6

    @test pia.indexer[1,1].n_local == 2
    @test pia.indexer[1,1].n_group1 == 2
    @test pia.indexer[1,1].start1 == 1
    @test pia.indexer[1,1].end1 == 2
    @test pia.indexer[1,1].n_group2 == 0
    @test pia.indexer[1,1].start2 == 0
    @test pia.indexer[1,1].end2 == -1

    @test pia.indexer[2,1].n_local == 0
    @test pia.indexer[2,1].n_group1 == 0
    @test pia.indexer[2,1].start1 == 0
    @test pia.indexer[2,1].end1 == -1
    @test pia.indexer[2,1].n_group2 == 0
    @test pia.indexer[2,1].start2 == 0
    @test pia.indexer[2,1].end2 == -1

    @test pia.indexer[3,1].n_local == 4
    @test pia.indexer[3,1].n_group1 == 3
    @test pia.indexer[3,1].n_group2 == 1
    @test pia.indexer[3,1].start1 == 5
    @test pia.indexer[3,1].end1 == 7
    @test pia.indexer[3,1].start2 == 10
    @test pia.indexer[3,1].end2 == 10

    # current state is [2, 0, 3], [0, 0, 1]
    # will become [2, 0, 1], [0, 0, 0]
    Merzbild.delete_batch_end!(particles, pia, 3, 1, 3)

    @test pia.index_last[1] == 5
    @test pia.n_total[1] == 3

    @test pia.indexer[1,1].n_local == 2
    @test pia.indexer[1,1].n_group1 == 2
    @test pia.indexer[1,1].start1 == 1
    @test pia.indexer[1,1].end1 == 2
    @test pia.indexer[1,1].n_group2 == 0
    @test pia.indexer[1,1].start2 == 0
    @test pia.indexer[1,1].end2 == -1

    @test pia.indexer[2,1].n_local == 0
    @test pia.indexer[2,1].n_group1 == 0
    @test pia.indexer[2,1].start1 == 0
    @test pia.indexer[2,1].end1 == -1
    @test pia.indexer[2,1].n_group2 == 0
    @test pia.indexer[2,1].start2 == 0
    @test pia.indexer[2,1].end2 == -1

    @test pia.indexer[3,1].n_local == 1
    @test pia.indexer[3,1].n_group1 == 1
    @test pia.indexer[3,1].n_group2 == 0
    @test pia.indexer[3,1].start1 == 5
    @test pia.indexer[3,1].end1 == 5
    @test pia.indexer[3,1].start2 == 0
    @test pia.indexer[3,1].end2 == -1

    Merzbild.delete_batch_end!(particles, pia, 3, 1, 3)

    @test pia.index_last[1] == 2
    @test pia.n_total[1] == 2

    @test pia.indexer[1,1].n_local == 2
    @test pia.indexer[1,1].n_group1 == 2
    @test pia.indexer[1,1].start1 == 1
    @test pia.indexer[1,1].end1 == 2
    @test pia.indexer[1,1].n_group2 == 0
    @test pia.indexer[1,1].start2 == 0
    @test pia.indexer[1,1].end2 == -1

    @test pia.indexer[2,1].n_local == 0
    @test pia.indexer[2,1].n_group1 == 0
    @test pia.indexer[2,1].start1 == 0
    @test pia.indexer[2,1].end1 == -1
    @test pia.indexer[2,1].n_group2 == 0
    @test pia.indexer[2,1].start2 == 0
    @test pia.indexer[2,1].end2 == -1

    @test pia.indexer[3,1].n_local == 0
    @test pia.indexer[3,1].n_group1 == 0
    @test pia.indexer[3,1].n_group2 == 0
    @test pia.indexer[3,1].start1 == 0
    @test pia.indexer[3,1].end1 == -1
    @test pia.indexer[3,1].start2 == 0
    @test pia.indexer[3,1].end2 == -1

    # test that index_last correctly falls back to 0 in find_index_last_after_group2_delete!
    # when the deleted group2 held the last index and no group1 particles exist anywhere
    pia, particles = generate_pia_and_particles(5, [0], [2])
    @test pia.index_last[1] == 2
    @test pia.n_total[1] == 2

    Merzbild.delete_particle_end_group2!(particles, pia, 1, 1)
    @test pia.index_last[1] == 1
    @test pia.n_total[1] == 1

    Merzbild.delete_particle_end_group2!(particles, pia, 1, 1)
    @test pia.index_last[1] == 0
    @test pia.n_total[1] == 0

    @test check_pia_is_correct(pia, 1) == (true, 0)

    # test that delete_batch_end_group1! and delete_batch_end_group2! are no-ops when n == 0
    # (delete_batch_end! itself never calls them with n == 0, so this exercises the guard directly)
    particles_zero = ParticleVector{0}(5)
    indexer_zero = ParticleIndexer(5)
    nbuffer_before = particles_zero.nbuffer

    Merzbild.delete_batch_end_group1!(particles_zero, indexer_zero, 0)
    @test indexer_zero.n_group1 == 5
    @test indexer_zero.n_local == 5
    @test particles_zero.nbuffer == nbuffer_before

    Merzbild.delete_batch_end_group2!(particles_zero, indexer_zero, 0)
    @test indexer_zero.n_group2 == 0
    @test indexer_zero.n_local == 5
    @test particles_zero.nbuffer == nbuffer_before
end