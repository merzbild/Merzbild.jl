@testset "octree_bounds_and_splitting" begin

    function create_particle_in_octant(octant, v_val, v0)
        # assume octants symmetric around v0
        # v_val > 0
        if octant >= 5
            v_z = v_val
        else
            v_z = -v_val
        end
        if octant % 2 == 1
            v_x = -v_val
        else
            v_x = v_val
        end
        if (octant == 3) || (octant == 4) || (octant == 7) || (octant == 8)
            v_y = v_val
        else
            v_y = -v_val
        end

        return Particle(1.0, [v_x, v_y, v_z] + v0, [0.0, 0.0, 0.0])
    end

    function create_15particles_nested()
        # suboctant symmetric around -2.0, 2.0, -2.0: that is octant # 3
        # so there particles have velocities as some combination of (-1.0, 1.0, -1.0), (-3.0, 3.0, -3.0)

        # mixed index array, created using Random.shuffle(1:15)
        # octants are
        #      3,  3,  3, 3,  3,  3, 3, 3, 8, 7,  6, 5, 4, 2, 1
        # or sorted:
        #      11, 3, 12, 15, 2, 14, 13, 7, 8, 10, 4, 9, 1, 6, 5
        mia = [12, 15, 2, 14, 13, 7, 8, 10, 5, 6, 1, 9, 4, 3, 11]

        vp = Vector{Particle}(undef, 15)
        i = 1
        for v_x in [-1.0, -3.0]
            for v_y in [1.0, 3.0]
                for v_z in [-3.0, -1.0]
                    vp[mia[i]] = Particle(1.0, [v_x, v_y, v_z], [0.0, 0.0, 0.0])
                    i += 1
                end
            end
        end

        rr = 8:-1:1

        for octant in rr
            if octant != 3
                vp[mia[i]] = create_particle_in_octant(octant, 1.0, [0.0, 0.0, 0.0])
                i += 1
            end
        end

        return vp
    end

    function create_8particles(v0)
        # create 1 particle per octant

        vp = Vector{Particle}(undef, 8)


        for i in 1:8
            vp[i] = create_particle_in_octant(i, 1.0, v0)
        end
        
        return vp
    end

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, "Ar")

    seed = 1234
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    particles1::Vector{Vector{Particle}} = [create_8particles([1.0, 1.0, 1.0])]
    pia1 = ParticleIndexerArray(8)
    octree1 = OctreeN2Merge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinC)

    # test bounds that are from -speed of light to +speed of light
    Merzbild.init_octree!(octree1, particles1[1], pia1, 1, 1)
    @test maximum(abs.(octree1.bins[1].v_min - [-299_792_458.0, -299_792_458.0, -299_792_458.0])) < 1e-12
    @test maximum(abs.(octree1.bins[1].v_max - [299_792_458.0, 299_792_458.0, 299_792_458.0])) < 1e-12
    # various tests for bin splitting, computing mean velocities, etc.

    Merzbild.split_bin!(octree1, 1, particles1[1])
    @test maximum(abs.(octree1.vel_middle)) < 1e-12  # should be almost 0 since everything is symmetric


    octree2 = OctreeN2Merge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVelSym)
    particles2::Vector{Vector{Particle}} = [create_15particles_nested()]
    pia2 = ParticleIndexerArray(15)
    Merzbild.init_octree!(octree2, particles2[1], pia2, 1, 1)
    @test maximum(abs.(octree2.bins[1].v_min + octree2.bins[1].v_max)) < 1e-12  # check that init bin is symmetric
    @test maximum(abs.(octree2.bins[1].v_max - [3.0, 3.0, 3.0])) < 1e-12  # check bin bounds

    # test bounds, non-symmetric octree bin
    octree3 = OctreeN2Merge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel)
    Merzbild.init_octree!(octree3, particles2[1], pia2, 1, 1)
    @test maximum(abs.(octree3.bins[1].v_min - [-3.0, -1.0, -3.0])) < 1e-11  # check bin bounds, non-symmetrized octree bin
    @test maximum(abs.(octree3.bins[1].v_max - [1.0, 3.0, 1.0])) < 1e-11  # check bin bounds, non-symmetrized octree bin


    # Test: refine bin in 3rd octant
    # test that the new suboctant bin bounds make sense


    octree4 = OctreeN2Merge(OctreeBinMeanSplit; init_bin_bounds=OctreeInitBinMinMaxVel)
    v_mean = SVector{3, Float64}(0.0, 0.0, 0.0)
    for i in 1:15
        v_mean = v_mean + particles2[1][i].v
    end
    v_mean = v_mean / 15.0

    Merzbild.init_octree!(octree4, particles2[1], pia2, 1, 1)
    Merzbild.compute_v_mean!(octree4, 1, 15, particles2[1])
    @test maximum(abs.(octree4.vel_middle - v_mean)) < 1e-11  # check computation of v_mean

    particles2[1][15].v = SVector{3, Float64}(120_000.0, -440_000.0, 920_000.0)

    vx_previous = octree4.bins[1].v_min[1]
    vy_previous = octree4.bins[1].v_max[2]
    vz_previous = octree4.bins[1].v_min[3]
    Merzbild.bin_bounds_recompute!(octree4, 1, 1, 15, particles2[1])

    @test octree4.bins[1].v_min[1] == vx_previous
    @test octree4.bins[1].v_max[1] == particles2[1][15].v[1]

    @test octree4.bins[1].v_min[2] == particles2[1][15].v[2]
    @test octree4.bins[1].v_max[2] == vy_previous

    @test octree4.bins[1].v_min[3] == vz_previous
    @test octree4.bins[1].v_max[3] == particles2[1][15].v[3]


    # Test: Merzbild.bin_bounds_inherit!
end