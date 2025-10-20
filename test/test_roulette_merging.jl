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
    # and we also don't conserve momentum and energy
    vp = ParticleVector(12)

    x_cell1 = [0.05, 0.05, 0.05, 0.45]
    for i in 1:4
        Merzbild.update_particle_buffer_new_particle!(vp, i)
        vp[i] = Particle(i, [0.5 - 0.1*i, 3.0 + 2.0*i^2, 4.0-0.3*i], [x_cell1[i], 0.0, 1.0])
    end

    x_cell2 = [0.55, 0.95, 0.8, 0.85, 0.99]
    for i in 5:5
        Merzbild.update_particle_buffer_new_particle!(vp, i)
        vp[i] = Particle(2.0, [-0.5+0.2*i^2, -3.0+(i*0.2)^3, -4.0+0.3*i], [x_cell2[i-4], 0.0, 2.0])
    end

    x_cell3 = [10.0, 11.0, 12.0]
    for i in 6:8
        Merzbild.update_particle_buffer_new_particle!(vp, i)
        vp[i] = Particle(3.0, [20.0-i, 20.0+0.5*i, 20.0-0.2*i], [x_cell3[i-5], 0.0, 3.0])
    end

    # add more particles to cell2
    for i in 9:12
        Merzbild.update_particle_buffer_new_particle!(vp, i)
        vp[i] = Particle(2.0, [-0.5-i, -3.0-0.05*i^2, -4.0+0.5*i], [x_cell2[i-7], 0.0, 2.0])
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

    weight_per_cell = [0.0 for i in 1:3]
    v0_per_cell = [zeros(3) for i in 1:3]
    E_per_cell = [0.0 for i in 1:3]
    for cell in 1:3
        wt = 0.0
        for i in pia.indexer[cell,1].start1:pia.indexer[cell,1].end1
            wt += vp[i].w
            v0_per_cell[cell] += vp[i].w * vp[i].v
        end
        for i in pia.indexer[cell,1].start2:pia.indexer[cell,1].end2
            wt += vp[i].w
            v0_per_cell[cell] += vp[i].w * vp[i].v
        end
        weight_per_cell[cell] = wt
        v0_per_cell[cell] /= wt

        for i in pia.indexer[cell,1].start1:pia.indexer[cell,1].end1
            E_per_cell[cell] += vp[i].w * sum((vp[i].v - v0_per_cell[cell]).^2)
        end
        for i in pia.indexer[cell,1].start2:pia.indexer[cell,1].end2
            E_per_cell[cell] += vp[i].w * sum((vp[i].v - v0_per_cell[cell]).^2)
        end
        E_per_cell[cell] /= wt
    end

    merge_roulette!(rng, vp, pia, 2, 1, 3)

    @test pia.indexer[2,1].n_local == 3
    @test pia.indexer[2,1].n_group2 > 0

    for cell in 1:3
        w_new = 0.0
        v0_new = zeros(3)
        for i in pia.indexer[cell,1].start1:pia.indexer[cell,1].end1
            w_new += vp[i].w
            v0_new += vp[i].w * vp[i].v
        end
        for i in pia.indexer[cell,1].start2:pia.indexer[cell,1].end2
            w_new += vp[i].w
            v0_new += vp[i].w * vp[i].v
        end
        v0_new /= w_new

        E_new = 0.0
        for i in pia.indexer[cell,1].start1:pia.indexer[cell,1].end1
            E_new += vp[i].w * sum((vp[i].v - v0_per_cell[cell]).^2)
        end
        for i in pia.indexer[cell,1].start2:pia.indexer[cell,1].end2
            E_new += vp[i].w * sum((vp[i].v - v0_per_cell[cell]).^2)
        end
        E_new /= w_new

        @test abs(w_new - weight_per_cell[cell]) < 2*eps()
        if cell == 2
            @test maximum(abs.(v0_new - v0_per_cell[cell])) > 1e-10
            @test abs(E_new - E_per_cell[cell]) > 1e-10
        else
            @test maximum(abs.(v0_new - v0_per_cell[cell])) < 2*eps()
            @test abs(E_new - E_per_cell[cell]) < 2*eps()
        end
    end

    # similar test but now we have the conservative version
    # and we also don't conserve momentum and energy
    vp = ParticleVector(12)

    x_cell1 = [0.05, 0.05, 0.05, 0.45]
    for i in 1:4
        Merzbild.update_particle_buffer_new_particle!(vp, i)
        vp[i] = Particle(i, [0.5 - 0.1*i, 3.0 + 2.0*i^2, 4.0-0.3*i], [x_cell1[i], 0.0, 1.0])
    end

    x_cell2 = [0.55, 0.95, 0.8, 0.85, 0.99]
    for i in 5:5
        Merzbild.update_particle_buffer_new_particle!(vp, i)
        vp[i] = Particle(2.0, [-0.5+0.2*i^2, -3.0+(i*0.2)^3, -4.0+0.3*i], [x_cell2[i-4], 0.0, 2.0])
    end

    x_cell3 = [10.0, 11.0, 12.0]
    for i in 6:8
        Merzbild.update_particle_buffer_new_particle!(vp, i)
        vp[i] = Particle(3.0+0.5*i, [20.0-i, 20.0+0.5*i, 20.0-0.2*i], [x_cell3[i-5], 0.0, 3.0])
    end

    # add more particles to cell2
    for i in 9:12
        Merzbild.update_particle_buffer_new_particle!(vp, i)
        vp[i] = Particle(2.0+i, [-0.5-i, -3.0-0.05*i^2, -4.0+0.5*i], [x_cell2[i-7], 0.0, 2.0])
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

    weight_per_cell = [0.0 for i in 1:3]
    v0_per_cell = [zeros(3) for i in 1:3]
    E_per_cell = [0.0 for i in 1:3]
    for cell in 1:3
        wt = 0.0
        for i in pia.indexer[cell,1].start1:pia.indexer[cell,1].end1
            wt += vp[i].w
            v0_per_cell[cell] += vp[i].w * vp[i].v
        end
        for i in pia.indexer[cell,1].start2:pia.indexer[cell,1].end2
            wt += vp[i].w
            v0_per_cell[cell] += vp[i].w * vp[i].v
        end
        weight_per_cell[cell] = wt
        v0_per_cell[cell] /= wt

        for i in pia.indexer[cell,1].start1:pia.indexer[cell,1].end1
            E_per_cell[cell] += vp[i].w * sum((vp[i].v - v0_per_cell[cell]).^2)
        end
        for i in pia.indexer[cell,1].start2:pia.indexer[cell,1].end2
            E_per_cell[cell] += vp[i].w * sum((vp[i].v - v0_per_cell[cell]).^2)
        end
        E_per_cell[cell] /= wt
    end

    merge_roulette!(rng, vp, pia, 2, 1, 3; conservative=true)

    @test pia.indexer[2,1].n_local == 3
    @test pia.indexer[2,1].n_group2 > 0

    for cell in 1:3
        w_new = 0.0
        v0_new = zeros(3)
        for i in pia.indexer[cell,1].start1:pia.indexer[cell,1].end1
            w_new += vp[i].w
            v0_new += vp[i].w * vp[i].v
        end
        for i in pia.indexer[cell,1].start2:pia.indexer[cell,1].end2
            w_new += vp[i].w
            v0_new += vp[i].w * vp[i].v
        end
        v0_new /= w_new

        E_new = 0.0
        for i in pia.indexer[cell,1].start1:pia.indexer[cell,1].end1
            E_new += vp[i].w * sum((vp[i].v - v0_per_cell[cell]).^2)
        end
        for i in pia.indexer[cell,1].start2:pia.indexer[cell,1].end2
            E_new += vp[i].w * sum((vp[i].v - v0_per_cell[cell]).^2)
        end
        E_new /= w_new

        @test abs(w_new - weight_per_cell[cell]) < 2*eps()
        @test maximum(abs.(v0_new - v0_per_cell[cell])) < 2e-15
        @test abs(E_new - E_per_cell[cell]) < 2*eps()
    end
end