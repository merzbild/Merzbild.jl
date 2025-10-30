@testset "debugging functions" begin

    function create_pia()
        pia_ = ParticleIndexerArray(4, 1)

        pia_.n_total[1] = 9

        pia_.indexer[1,1].n_local = 5
        pia_.indexer[1,1].n_group1 = 3
        pia_.indexer[1,1].start1 = 1
        pia_.indexer[1,1].end1 = 3
        pia_.indexer[1,1].n_group2 = 2
        pia_.indexer[1,1].start2 = 4
        pia_.indexer[1,1].end2 = 5

        pia_.indexer[3,1].n_local = 3
        pia_.indexer[3,1].n_group2 = 3
        pia_.indexer[3,1].start2 = 6
        pia_.indexer[3,1].end2 = 8

        pia_.indexer[4,1].n_local = 1
        pia_.indexer[4,1].n_group1 = 1
        pia_.indexer[4,1].start1 = 9
        pia_.indexer[4,1].end1 = 9
        return pia_
    end

    particles = [ParticleVector(10)]
    for i in 1:8
        Merzbild.add_particle!(particles[1], i, i*1.0, [2.0, 2.0, 3.0], [11.0, 12.0, 14.0])
    end

    pia = ParticleIndexerArray(1, 1)

    pia.n_total[1] = 8

    pia.indexer[1,1].n_group1 = 8
    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = 8

    # should be OK
    @test check_unique_index(particles[1], pia, 1) == (true, 0)

    # test that buffer particles are checked correctly
    particles[1].nbuffer = 3
    particles[1].buffer[3] = particles[1].index[1]
    @test check_unique_index(particles[1], pia, 1) == (false, -1)

    # add 2 other indices pointing to particle 3
    particles[1].index[4] = 3
    particles[1].index[5] = 3
    @test check_unique_index(particles[1], pia, 1) == (false, -3)

    # fix buffer and one index
    particles[1].nbuffer = 2
    particles[1].index[5] = 5
    @test check_unique_index(particles[1], pia, 1) == (false, 2)

    # fix last index
    particles[1].index[4] = 4
    @test check_unique_index(particles[1], pia, 1) == (true, 0)

    
    particles = [ParticleVector(10)]
    for i in 1:8
        Merzbild.add_particle!(particles[1], i, i*1.0, [2.0, 2.0, 3.0], [11.0, 12.0, 14.0])
    end
    pia = ParticleIndexerArray(1, 1)

    pia.n_total[1] = 8

    pia.indexer[1,1].n_group1 = 4
    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = 4
    pia.indexer[1,1].n_group2 = 4
    pia.indexer[1,1].start2 = 5
    pia.indexer[1,1].end2 = 8

    # should be OK
    @test check_unique_index(particles[1], pia, 1) == (true, 0)

    # test pia debugging
    pia = create_pia()

    @test check_pia_is_correct(pia, 1) == (true, 0)

    # n_total is not correct
    pia.n_total[1] = 8

    @test check_pia_is_correct(pia, 1) == (false, 0)

    # incorrect start of empty group
    pia.n_total[1] = 9
    pia.indexer[4,1].start2 = 10
    @test check_pia_is_correct(pia, 1) == (false, 4)

    pia.indexer[4,1].start2 = 0

    # incorrect group count
    pia.indexer[3,1].n_group2 = 2
    @test check_pia_is_correct(pia, 1) == (false, 3)

    pia.indexer[3,1].n_group2 = 3
    # incorrect n_local
    pia.indexer[1,1].n_local = 6
    @test check_pia_is_correct(pia, 1) == (false, 1)

    pia = create_pia()

    # n_group1 = 0 but start1 != 0, end1 != -1
    pia.indexer[3,1].start1 = 2
    pia.indexer[3,1].end1 = 3
    @test check_pia_is_correct(pia, 1) == (false, 3)

    # n_group1 is not equal to e1 - s1 + 1
    pia = create_pia()
    pia.indexer[4,1].n_group2 = 2
    @test check_pia_is_correct(pia, 1) == (false, 4)

    # n_group2 is not equal to e2 - s2 + 1
    pia = create_pia()
    pia.indexer[3,1].n_group2 = 2
    @test check_pia_is_correct(pia, 1) == (false, 3)

    # test pretty_print_pia
    pia = create_pia()

    out = stdout
    stream = redirect_stdout()
    
    pretty_print_pia(pia, 1)

    redirect_stdout(out)
    @test String(readavailable(stream)) == """
    Total: 9
    Cell 1: group1: [1, 3] group2: [4, 5]
    Cell 3: group2: [6, 8]
    Cell 4: group1: [9, 9]
    """
end