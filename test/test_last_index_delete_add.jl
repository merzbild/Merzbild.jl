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
    # pretty_print_pia(pia, 1)

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
    pretty_print_pia(pia, 1)
    println(particles.buffer[1:particles.nbuffer])

    # create particle in cell 7 in group 2
    # group2 was [38, 41] before
    Merzbild.update_buffer_index_new_particle!(particles, pia, 7, 1)

    @test pia.indexer[7,1].n_group2 == 5
    @test pia.indexer[7,1].start2 == 38
    @test pia.indexer[7,1].end2 == 42
    @test pia.index_last[1] == 42
    @test pia.n_total[1] == 36

    println(particles.index[42])

    println(particles.buffer[1:particles.nbuffer])

    @test check_pia_is_correct(pia, 1) == (true, 0)
    @test check_unique_index(particles, pia, 1) == (true, 0)
    @test check_unique_buffer(particles) == (true, 0)
end