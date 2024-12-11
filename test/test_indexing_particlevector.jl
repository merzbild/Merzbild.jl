@testset "particle vector and multidimensional indexing" begin
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
    @test particles[1].buffer == [11,12,13,10,9,8,7,6,5,4,3,2,1]
end