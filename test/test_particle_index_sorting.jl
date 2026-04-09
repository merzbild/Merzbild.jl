@testset "particle index re-sorting" begin
    # tests for restore_particle_ordering!

    
    particles = ParticleVector(10)
    
    # test that for non-initialized particles the buffer is correct
    @test particles.index == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    @test particles.buffer == [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
    @test particles.nbuffer == 10

    inv_map = zeros(Int64, 4)

    @testset "test that used_indices are resized, in case of ordered input nothing changes" begin
        for i in 1:10
            particles[i] = Particle(i, [i, 2 * i, 3 * i], [-i, -4 * i, -5 * i])
        end

        restore_particle_ordering!(particles, inv_map)

        # check that nothing changed
        @test particles.index == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        @test particles.buffer == [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
        @test particles.nbuffer == 10
        @test length(inv_map) == 10

        for i in 1:10
            @test particles[i].w == i
            @test particles[i].v[1] == i
            @test particles[i].v[2] == 2 * i
            @test particles[i].v[3] == 3 * i
            @test particles[i].x[1] == -i
            @test particles[i].x[2] == -4 * i
            @test particles[i].x[3] == -5 * i
        end
    end

    @testset "test that indexing and buffer are correctly re-ordered 1" begin
        particles.index = [5, 1, 7, 3, 2, 4, 10, 8, 6, 9]

        # 2 particles in use, 8 in buffer
        for i in 1:2
            particles[i] = Particle(i, [i, 2 * i, 3 * i], [-i, -4 * i, -5 * i])
        end

        # these particles are not in use
        for i in 3:10
            particles[i] = Particle(0, [0, 0, 0], [0, 0, 0])
        end

        particles.buffer = [7, 3, 2, 4, 10, 8, 6, 9, 5, 1]  # only first 8 elements matter
        particles.nbuffer = 8

        restore_particle_ordering!(particles, inv_map)

        # check that everything is correct
        @test particles.index == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

        # test only relevant part of buffer
        @test particles.buffer[1:8] == [10, 9, 8, 7, 6, 5, 4, 3]
        @test particles.nbuffer == 8

        for i in 1:2
            @test particles[i].w == i
            @test particles[i].v[1] == i
            @test particles[i].v[2] == 2 * i
            @test particles[i].v[3] == 3 * i
            @test particles[i].x[1] == -i
            @test particles[i].x[2] == -4 * i
            @test particles[i].x[3] == -5 * i
        end
    end

    @testset "test that indexing and buffer are correctly re-ordered 2" begin
        particles.index = [10, 9, 5, 6, 1, 8, 3, 7, 2, 4]

        # 7 particles in use, 3 in buffer
        for i in 1:7
            particles[i] = Particle(i, [i, 2 * i, 3 * i], [-i, -4 * i, -5 * i])
        end

        # these particles are not in use
        for i in 8:10
            particles[i] = Particle(0, [0, 0, 0], [0, 0, 0])
        end

        particles.buffer = [4, 7, 2, 10, 9, 5, 6, 1, 8, 3]  # only first 3 elements matter
        particles.nbuffer = 3

        restore_particle_ordering!(particles, inv_map)

        # check that everything is correct
        @test particles.index == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

        # test only relevant part of buffer
        @test particles.buffer[1:3] == [10, 9, 8]
        @test particles.nbuffer == 3

        for i in 1:7
            @test particles[i].w == i
            @test particles[i].v[1] == i
            @test particles[i].v[2] == 2 * i
            @test particles[i].v[3] == 3 * i
            @test particles[i].x[1] == -i
            @test particles[i].x[2] == -4 * i
            @test particles[i].x[3] == -5 * i
        end
    end

    @testset "test that indexing and buffer are correctly re-ordered 3" begin
        particles.index = [4, 10, 7, 3, 1, 2, 9, 8, 5, 6]

        # all 10 particles in use, 0 in buffer
        for i in 1:10
            particles[i] = Particle(i, [i, 2 * i, 3 * i], [-i, -4 * i, -5 * i])
        end

        particles.buffer = [4, 7, 2, 10, 9, 5, 6, 1, 8, 3]  # doesn't matter
        particles.nbuffer = 0

        restore_particle_ordering!(particles, inv_map)

        # check that everything is correct
        @test particles.index == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

        # test only relevant part of buffer
        @test particles.nbuffer == 0

        for i in 1:10
            @test particles[i].w == i
            @test particles[i].v[1] == i
            @test particles[i].v[2] == 2 * i
            @test particles[i].v[3] == 3 * i
            @test particles[i].x[1] == -i
            @test particles[i].x[2] == -4 * i
            @test particles[i].x[3] == -5 * i
        end
    end

    @testset "test that indexing and buffer are correctly re-ordered 4" begin
        particles.index = [10, 5, 7, 6, 2, 1, 3, 5, 4, 6]

        # 5 particles in use, 5 in buffer, 5 last elements of index are duplicate rubbish
        for i in 1:5
            particles[i] = Particle(i, [i, 2 * i, 3 * i], [-i, -4 * i, -5 * i])
        end

        particles.buffer = [3, 7, 4, 1, 9, 10, 8, 2, 6, 5]  # doesn't matter
        particles.nbuffer = 5

        restore_particle_ordering!(particles, inv_map)

        # check that everything is correct
        @test particles.index == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

        # test only relevant part of buffer
        @test particles.nbuffer == 5
        @test particles.buffer[1:5] == [10, 9, 8, 7, 6]

        for i in 1:5
            @test particles[i].w == i
            @test particles[i].v[1] == i
            @test particles[i].v[2] == 2 * i
            @test particles[i].v[3] == 3 * i
            @test particles[i].x[1] == -i
            @test particles[i].x[2] == -4 * i
            @test particles[i].x[3] == -5 * i
        end
    end
end