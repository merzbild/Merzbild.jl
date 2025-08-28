@testset "particle vector and multidimensional indexing" begin
    # the tests of indexing ParticleVectors and the buffer part of the ParticleVectors
    # are not that easy to split nicely, so some duplication of tests is possible
    
    pia = ParticleIndexerArray(8, 3)  # 8 cells 3 species
    @test size(pia.indexer) == (8, 3)

    particles = [ParticleVector(10)]

    @test particles[1].index == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    @test length(particles[1].particles) == 10

    for i in 1:10
        particles[1].particles[i] = Particle(i, [0.0, 0.0, 0.0], [10.0 - i, 0.0, 1.0])
    end

    for i in 1:10
        @test particles[1][i].w == i
        @test particles[1][i].x[1] == 10.0 - i
    end

    # now we reverse index
    particles[1].index = [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]

    # check that we now get the particles in a different order
    # but the underlying particle vector hasn't changed
    for i in 1:10
        @test particles[1][i].w == 11.0 - i
        @test particles[1][i].x[1] == i - 1.0

        @test particles[1].particles[i].w == i
        @test particles[1].particles[i].x[1] == 10.0 - i

        @test particles[1].cell[i] == 0
    end

    old_buffer_length = particles[1].nbuffer
    resize!(particles[1], length(particles[1]) + 3)
    @test length(particles[1]) == 13
    @test length(particles[1].particles) == 13
    @test length(particles[1].index) == 13
    @test length(particles[1].cell) == 13
    @test length(particles[1].buffer) == 13

    @test particles[1].nbuffer == old_buffer_length + 3
    @test particles[1].buffer == [13,12,11,10,9,8,7,6,5,4,3,2,1]

    particles = [ParticleVector(3)]
    @test_throws UndefRefError particles[1][1]
    @test_throws UndefRefError particles[1][2]
    @test_throws UndefRefError particles[1][3]

    Merzbild.add_particle!(particles[1], 1, 20.0, [2.0, 2.0, 3.0], [11.0, 12.0, 14.0])

    @test_throws UndefRefError particles[1][2]
    @test_throws UndefRefError particles[1][3]

    @test particles[1][1].w == 20.0
    @test particles[1][1].v[1] == 2.0
    @test particles[1][1].v[2] == 2.0
    @test particles[1][1].v[3] == 3.0
    @test particles[1][1].x[1] == 11.0
    @test particles[1][1].x[2] == 12.0
    @test particles[1][1].x[3] == 14.0

    Merzbild.add_particle!(particles[1], 2, 40.0, [2.0, 2.0, 3.0], [11.0, 12.0, 14.0])
    @test_throws UndefRefError particles[1][3]

    @test particles[1][1].w == 20.0
    @test particles[1][2].w == 40.0
    @test particles[1].nbuffer == 1

    particles = [ParticleVector(8)]
    for i in 1:4
        Merzbild.add_particle!(particles[1], i, i*1.0, [2.0, 2.0, 3.0], [11.0, 12.0, 14.0])
    end

    pia = ParticleIndexerArray(1, 1)
    pia.n_total[1] = 4
    pia.indexer[1,1].n_group1 = 1
    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = 4

    @test count_disordered_particles(particles[1], pia, 1; use_offset=false) == 0

    particles[1].index[2] = 3
    particles[1].index[3] = 2
    @test particles[1][1].w == 1
    @test particles[1][2].w == 3
    @test particles[1][3].w == 2
    @test particles[1][4].w == 4

    @test count_disordered_particles(particles[1], pia, 1; use_offset=false) == 2

    particles[1].index[1] = 4
    particles[1].index[4] = 1

    # reduce n_total
    pia.n_total[1] = 3
    @test count_disordered_particles(particles[1], pia, 1; use_offset=false) == 3


    particles = [ParticleVector(10)]
    for i in 1:10
        Merzbild.add_particle!(particles[1], i, i*1.0, [2.0, 2.0, 3.0], [11.0, 12.0, 14.0])
    end


    pia = ParticleIndexerArray(3, 1)

    pia.n_total[1] = 10

    pia.indexer[1,1].n_group1 = 3
    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = 3

    pia.indexer[2,1].n_group1 = 4
    pia.indexer[2,1].start1 = 4
    pia.indexer[2,1].end1 = 7

    pia.indexer[3,1].n_group1 = 2
    pia.indexer[3,1].start1 = 8
    pia.indexer[3,1].end1 = 10

    # set indexing so that group 1 and 2 still work correctly
    particles[1].index = [1, 2, 3, 5, 6, 7, 8, 4, 10, 9]

    # without offset only 3 placed correctly: 10-3=6
    @test count_disordered_particles(particles[1], pia, 1; use_offset=false) == 7

    # with offset 7 placed correctly
    # but we skip the first particle of the cell in the offset-based approach
    # so it's only 2 particles out-of-place
    @test count_disordered_particles(particles[1], pia, 1; use_offset=true) == 2
end