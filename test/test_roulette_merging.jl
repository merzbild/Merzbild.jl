@testset "roulette merge test" begin
    seed = 1234
    rng = StableRNG(seed)

    vp = ParticleVector(12)

    x_cell1 = [0.05, 0.05, 0.05, 0.45]
    for i in 1:4
        Merzbild.update_particle_buffer_new_particle!(vp, i)
        vp[i] = Particle(i, [0.5, 3.0, 4.0], [x_cell1[i], 0.0, 1.0])
    end

    x_cell2 = [0.55, 0.95, 0.8, 0.85, 0.99]
    for i in 5:9
        Merzbild.update_particle_buffer_new_particle!(vp, i)
        vp[i] = Particle(3.0, [-0.5, -3.0, -4.0], [x_cell2[i-4], 0.0, 2.0])
    end

    x_cell3 = [10.0, 11.0, 12.0]
    for i in 10:12
        Merzbild.update_particle_buffer_new_particle!(vp, i)
        vp[i] = Particle(3.0, [20.0, 20.0, 20.0], [x_cell3[i-9], 0.0, 3.0])
    end

    pia = ParticleIndexerArray(3, 1)

    pia.n_total[1] = 12

    pia.indexer[1,1].n_local = 4
    pia.indexer[1,1].n_group1 = 4
    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = 4

    pia.indexer[1,1].n_group2 = 0
    pia.indexer[1,1].start2 = 0
    pia.indexer[1,1].end2 = -1

    pia.indexer[2,1].n_local = 5
    pia.indexer[2,1].n_group1 = 5
    pia.indexer[2,1].start1 = 5
    pia.indexer[2,1].end1 = 9

    pia.indexer[2,1].n_group2 = 0
    pia.indexer[2,1].start2 = 0
    pia.indexer[2,1].end2 = -1

    pia.indexer[3,1].n_local = 3
    pia.indexer[3,1].n_group1 = 3
    pia.indexer[3,1].start1 = 10
    pia.indexer[3,1].end1 = 12

    pia.indexer[3,1].n_group2 = 0
    pia.indexer[3,1].start2 = 0
    pia.indexer[3,1].end2 = -1

    pia.contiguous[1] = true

    merge_roulette!(rng, vp, pia, 2, 1, 3)

    @test pia.contiguous[1] == false

    @test pia.indexer[1,1].end1 == 4

    @test pia.indexer[2,1].n_local == 3
    @test pia.indexer[2,1].n_group1 == 3
    @test pia.indexer[2,1].start1 == 5
    @test pia.indexer[2,1].end1 == 7

    @test pia.indexer[2,1].n_group2 == 0
    @test pia.indexer[2,1].start2 == 0
    @test pia.indexer[2,1].end2 == -1

    @test pia.indexer[3,1].start1 == 10

    for i in pia.indexer[1,1].start1:pia.indexer[1,1].end1
        @test vp[i].w == i
        @test vp[i].v[1] == 0.5
        @test vp[i].x[3] == 1.0
    end
    for i in pia.indexer[3,1].start1:pia.indexer[3,1].end1
        @test vp[i].w == 3.0
        @test vp[i].v[1] == 20.0
        @test vp[i].x[3] == 3.0
    end
    for i in pia.indexer[2,1].start1:pia.indexer[2,1].end1
        @test abs(vp[i].w - 5.0) < 1e-15
        @test vp[i].v[1] == -0.5  # velocity unchanged
        @test vp[i].x[1] in x_cell2
        @test vp[i].x[3] == 2.0
    end


    vp = ParticleVector(12)

    x_cell1 = [0.05, 0.05, 0.05, 0.45]
    for i in 1:4
        Merzbild.update_particle_buffer_new_particle!(vp, i)
        vp[i] = Particle(i, [0.5, 3.0, 4.0], [x_cell1[i], 0.0, 1.0])
    end

    x_cell2 = [0.55, 0.95, 0.8, 0.85, 0.99]
    for i in 5:5
        Merzbild.update_particle_buffer_new_particle!(vp, i)
        vp[i] = Particle(2.0, [-0.5, -3.0, -4.0], [x_cell2[i-4], 0.0, 2.0])
    end

    x_cell3 = [10.0, 11.0, 12.0]
    for i in 6:8
        Merzbild.update_particle_buffer_new_particle!(vp, i)
        vp[i] = Particle(3.0, [20.0, 20.0, 20.0], [x_cell3[i-5], 0.0, 3.0])
    end

    for i in 9:12
        Merzbild.update_particle_buffer_new_particle!(vp, i)
        vp[i] = Particle(2.0, [-0.5, -3.0, -4.0], [x_cell2[i-7], 0.0, 2.0])
    end

    pia = ParticleIndexerArray(3, 1)

    pia.n_total[1] = 12

    pia.indexer[1,1].n_local = 4
    pia.indexer[1,1].n_group1 = 4
    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = 4

    pia.indexer[1,1].n_group2 = 0
    pia.indexer[1,1].start2 = 0
    pia.indexer[1,1].end2 = -1

    pia.indexer[2,1].n_local = 5
    pia.indexer[2,1].n_group1 = 1
    pia.indexer[2,1].start1 = 5
    pia.indexer[2,1].end1 = 5

    pia.indexer[2,1].n_group2 = 4
    pia.indexer[2,1].start2 = 9
    pia.indexer[2,1].end2 = 12

    pia.indexer[3,1].n_local = 3
    pia.indexer[3,1].n_group1 = 3
    pia.indexer[3,1].start1 = 6
    pia.indexer[3,1].end1 = 8

    pia.indexer[3,1].n_group2 = 0
    pia.indexer[3,1].start2 = 0
    pia.indexer[3,1].end2 = -1

    pia.contiguous[1] = true

    merge_roulette!(rng, vp, pia, 2, 1, 1)

    @test pia.contiguous[1] == false

    @test pia.indexer[1,1].end1 == 4

    @test pia.indexer[2,1].n_local == 1
    @test pia.indexer[2,1].n_group1 == 1
    @test pia.indexer[2,1].start1 == 5
    @test pia.indexer[2,1].end1 == 5

    @test pia.indexer[2,1].n_group2 == 0
    @test pia.indexer[2,1].start2 == 0
    @test pia.indexer[2,1].end2 == -1

    @test pia.indexer[3,1].start1 == 6

    for i in pia.indexer[1,1].start1:pia.indexer[1,1].end1
        @test vp[i].w == i
        @test vp[i].v[1] == 0.5
        @test vp[i].x[3] == 1.0
    end
    for i in pia.indexer[3,1].start1:pia.indexer[3,1].end1
        @test vp[i].w == 3.0
        @test vp[i].v[1] == 20.0
        @test vp[i].x[3] == 3.0
    end
    for i in pia.indexer[2,1].start1:pia.indexer[2,1].end1
        @test abs(vp[i].w - 10.0) < 1e-15
        @test vp[i].v[1] == -0.5  # velocity unchanged
        @test vp[i].x[1] in x_cell2
        @test vp[i].x[3] == 2.0
    end

    # check that if we have post-merge particles in group2, the weight is still conserved
    vp = ParticleVector(12)

    x_cell1 = [0.05, 0.05, 0.05, 0.45]
    for i in 1:4
        Merzbild.update_particle_buffer_new_particle!(vp, i)
        vp[i] = Particle(i, [0.5, 3.0, 4.0], [x_cell1[i], 0.0, 1.0])
    end

    x_cell2 = [0.55, 0.95, 0.8, 0.85, 0.99]
    for i in 5:5
        Merzbild.update_particle_buffer_new_particle!(vp, i)
        vp[i] = Particle(2.0, [-0.5, -3.0, -4.0], [x_cell2[i-4], 0.0, 2.0])
    end

    x_cell3 = [10.0, 11.0, 12.0]
    for i in 6:8
        Merzbild.update_particle_buffer_new_particle!(vp, i)
        vp[i] = Particle(3.0, [20.0, 20.0, 20.0], [x_cell3[i-5], 0.0, 3.0])
    end

    for i in 9:12
        Merzbild.update_particle_buffer_new_particle!(vp, i)
        vp[i] = Particle(2.0, [-0.5, -3.0, -4.0], [x_cell2[i-7], 0.0, 2.0])
    end

    pia = ParticleIndexerArray(3, 1)

    pia.n_total[1] = 12

    pia.indexer[1,1].n_local = 4
    pia.indexer[1,1].n_group1 = 4
    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = 4

    pia.indexer[1,1].n_group2 = 0
    pia.indexer[1,1].start2 = 0
    pia.indexer[1,1].end2 = -1

    pia.indexer[2,1].n_local = 5
    pia.indexer[2,1].n_group1 = 1
    pia.indexer[2,1].start1 = 5
    pia.indexer[2,1].end1 = 5

    pia.indexer[2,1].n_group2 = 4
    pia.indexer[2,1].start2 = 9
    pia.indexer[2,1].end2 = 12

    pia.indexer[3,1].n_local = 3
    pia.indexer[3,1].n_group1 = 3
    pia.indexer[3,1].start1 = 6
    pia.indexer[3,1].end1 = 8

    pia.indexer[3,1].n_group2 = 0
    pia.indexer[3,1].start2 = 0
    pia.indexer[3,1].end2 = -1

    pia.contiguous[1] = true

    merge_roulette!(rng, vp, pia, 2, 1, 3)

    @test pia.indexer[2,1].n_local == 3
    @test pia.indexer[2,1].n_group2 > 0

    new_w = 0.0
    for i in pia.indexer[2,1].start1:pia.indexer[2,1].end1
        new_w += vp[i].w
    end
    for i in pia.indexer[2,1].start2:pia.indexer[2,1].end2
        new_w += vp[i].w
    end
    @test abs(new_w - 10.0) < 2*eps()
end