@testset "octree_sorting" begin

    function create_particle_in_octant(octant, v_val; w=1.0)
        # assume octants symmetric around (0, 0, 0)
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

        return Particle(w, [v_x, v_y, v_z], [0.0, 0.0, 0.0])
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
                    vp[mia[i]] = Particle(0.5, [v_x, v_y, v_z], [0.0, 0.0, 0.0])
                    i += 1
                end
            end
        end

        rr = 8:-1:1

        for octant in rr
            if octant != 3
                vp[mia[i]] = create_particle_in_octant(octant, 1.0)
                i += 1
            end
        end

        return vp
    end

    function create_9particles(extra_octant, reverse)
        # create 1 particle per octant, plus one extra particle in specified octant

        vp = Vector{Particle}(undef, 9)

        i = 1

        if reverse
            rr = 8:-1:1 # iterate in reverse so that we actually have to sort the particles
        else
            rr = 1:8
        end

        for octant in rr
            vp[i] = create_particle_in_octant(octant, 1.0, w=octant)
            i += 1
            if octant == extra_octant
                vp[i] = create_particle_in_octant(octant, 2.0, w=octant)
                i += 1
            end
        end
        
        return vp
    end

    function create_7particles(missing_octant, reverse)
        # create 1 particle per octant, except one empty octant
        
        vp = Vector{Particle}(undef, 7)

        i = 1

        if reverse
            rr = 8:-1:1 # iterate in reverse so that we actually have to sort the particles
        else
            rr = 1:8
        end

        for octant in rr
            if octant != missing_octant
                vp[i] = create_particle_in_octant(octant, 1.0, w=octant)
                i += 1
            end
        end
        
        return vp
    end

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, "Ar")

    seed = 1234
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    # first we test function that computes octants
    vmid = SVector{3, Float64}(-1.0, -1.0, -1.0)
    
    v0 = SVector{3, Float64}(-2.0, -2.0, -2.0)  # - - -
    @test Merzbild.compute_octant(v0, vmid) == 1

    v0 = SVector{3, Float64}(2.0, -2.0, -2.0)  # + - -
    @test Merzbild.compute_octant(v0, vmid) == 2
    
    v0 = SVector{3, Float64}(-2.0, 2.0, -2.0)  # - + -
    @test Merzbild.compute_octant(v0, vmid) == 3
    
    v0 = SVector{3, Float64}(2.0, 2.0, -2.0)  # + + -
    @test Merzbild.compute_octant(v0, vmid) == 4

    v0 = SVector{3, Float64}(-2.0, -2.0, 2.0)  # - - +
    @test Merzbild.compute_octant(v0, vmid) == 5

    v0 = SVector{3, Float64}(2.0, -2.0, 2.0)  # + - +
    @test Merzbild.compute_octant(v0, vmid) == 6
    
    v0 = SVector{3, Float64}(-2.0, 2.0, 2.0)  # - + +
    @test Merzbild.compute_octant(v0, vmid) == 7

    v0 = SVector{3, Float64}(2.0, 2.0, 2.0)  # + + +
    @test Merzbild.compute_octant(v0, vmid) == 8

    # then we test particle sorting with 2 particles in a specific octant
    # octants are 8, 7, 6, 5, 4, 3, 3, 2, 1
    particles9::Vector{Vector{Particle}} = [create_9particles(3, true)]
    pia = ParticleIndexerArray(9)
    octree = OctreeN2Merge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinC)

    Merzbild.init_octree!(octree, particles9[1], pia, 1, 1)
    for i in 1:9
        @test i in octree.particle_indexes_sorted[1:9]
    end
    
    Merzbild.split_bin!(octree, 1, particles9[1])
    @test octree.Nbins == 8
    for i in 1:9
        @test i in octree.particle_indexes_sorted[1:9]
    end

    for i in 1:9
        @test octree.particle_indexes_sorted[i] == 9 - i + 1
    end

    # more convoluted test case for the sorting
    # particles:  1, 2, 3, 4, 5, 6, 7, 8, 9
    # octants are 4, 5, 6, 7, 8, 3, 3, 2, 1
    # end results should be 9, 8, 7, 6, 1, 2, 3, 4, 5
    particles9[1][1:5] = create_9particles(3, false)[5:9]
    particles9[1][6:9] = create_9particles(3, true)[6:9]
    Merzbild.init_octree!(octree, particles9[1], pia, 1, 1)
    Merzbild.split_bin!(octree, 1, particles9[1])
    @test octree.Nbins == 8
    expected = [9, 8, 7, 6, 1, 2, 3, 4, 5]
    @test octree.particle_indexes_sorted[1:9] == expected

    @test octree.bin_start[1:8] == [1, 2, 3, 5, 6, 7, 8, 9]
    @test octree.bin_end[1:8] == [1, 2, 4, 5, 6, 7, 8, 9]

    # the same but with pia constructed by hand with particle indices split into 2 groups
    pia2 = ParticleIndexerArray(3)
    pia2.n_total[1] = 9
    pia2.indexer[1,1].n_local = 9

    pia2.indexer[1,1].start1 = 1
    pia2.indexer[1,1].end1 = 3
    pia2.indexer[1,1].n_group1 = 3
    
    pia2.indexer[1,1].start2 = 4
    pia2.indexer[1,1].end2 = 9
    pia2.indexer[1,1].n_group2 = 6

    Merzbild.init_octree!(octree, particles9[1], pia2, 1, 1)
    Merzbild.split_bin!(octree, 1, particles9[1])

    @test octree.Nbins == 8
    expected = [9, 8, 7, 6, 1, 2, 3, 4, 5]
    @test octree.particle_indexes_sorted[1:9] == expected
    @test octree.bin_start[1:8] == [1, 2, 3, 5, 6, 7, 8, 9]
    @test octree.bin_end[1:8] == [1, 2, 4, 5, 6, 7, 8, 9]

    num_per_bin = [1, 1, 2, 1, 1, 1, 1, 1]
    w_per_bin = [1, 2, 6, 4, 5, 6, 7, 8]
    for i in 1:8
        @test octree.bins[i].np == num_per_bin[i]
        @test octree.bins[i].w == w_per_bin[i]
    end

    # now we test for 7 particles, i.e. with one empty octree bin
    # particles:  1, 2, 3, 4, 5, 6, 7
    # octants are 8, 7, 6, 4, 3, 2, 1
    particles7::Vector{Vector{Particle}} = [create_7particles(5, true)]
    pia3 = ParticleIndexerArray(7)

    Merzbild.init_octree!(octree, particles7[1], pia3, 1, 1)
    Merzbild.split_bin!(octree, 1, particles7[1])
    @test octree.Nbins == 7
    expected = [7, 6, 5, 4, 3, 2, 1]
    @test octree.particle_indexes_sorted[1:7] == expected
    @test octree.bin_start[1:7] == [1, 2, 3, 4, 5, 6, 7]
    @test octree.bin_end[1:7] == [1, 2, 3, 4, 5, 6, 7]

    num_per_bin = [1, 1, 1, 1, 1, 1, 1]
    w_per_bin = [1, 2, 3, 4, 6, 7, 8]
    for i in 1:7
        @test octree.bins[i].np == num_per_bin[i]
        @test octree.bins[i].w == w_per_bin[i]
    end


    # test double-splitting, i.e. we have 1 p/bin in 7 bins, the 8th bin has 4 particles
    # we do splitting of the main 0 bin, and then do a second split
    
    particles15::Vector{Vector{Particle}} = [create_15particles_nested()]
    pia4 = ParticleIndexerArray(15)

    octree2 = OctreeN2Merge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVelSym)

    # the initial bin will have bounds [-3.0, -3.0, -3.0], [3.0, 3.0, 3.0] 
    Merzbild.init_octree!(octree2, particles15[1], pia4, 1, 1)

    # test bounds, symmetric octree bin
    @test octree2.bin_start[1] == 1
    @test octree2.bin_end[1] == 15

    Merzbild.split_bin!(octree2, 1, particles15[1])
    @test octree2.Nbins == 8


    # we take the following octants and assign them to corresponding particle indices from mia
    #        3,  3,  3, 3,  3,  3, 3, 3,  1, 2, 4, 5, 6, 7, 8
    # mia = [12, 15, 2, 14, 13, 7, 8, 10, 5, 6, 1, 9, 4, 3, 11]
    # 
    # so sorted array would be:
    #      11, 3, 12, 15, 2, 14, 13, 7, 8, 10, 4, 9, 1, 6, 5
    expected12 = [11, 3] # in octants 1/2
    expected1115 = [4, 9, 1, 6, 5] # in octants 1/2

    @test octree2.particle_indexes_sorted[1:2] == expected12
    @test octree2.bin_start[3] == 3 
    @test octree2.bin_end[3] == 10

    octind3 = sort(copy(octree2.particle_indexes_sorted[3:10]))  # sort for easier comparison
    @test octind3 == [2, 7, 8, 10, 12, 13, 14, 15]

    @test octree2.particle_indexes_sorted[11:15] == expected1115

    num_per_bin = [1, 1, 8, 1, 1, 1, 1, 1]
    w_per_bin = [1, 1, 4.0, 1, 1, 1, 1, 1]
    for i in 1:8
        @test octree2.bins[i].np == num_per_bin[i]
        @test octree2.bins[i].w == w_per_bin[i]
    end
    # the sub-bin will have bounds [-3.0, 0.0], [0.0, 3.0], [-3.0, 0.0]
    # and will be split along velocity of [-1.5, 1.5, 1.5]
    # each sub-bin will have 1 particle
    Merzbild.split_bin!(octree2, 3, particles15[1])
    @test octree2.Nbins == 15

    # particle indices stayed in the same order
    @test octree2.particle_indexes_sorted[1:2] == expected12
    octind3 = sort(copy(octree2.particle_indexes_sorted[3:10]))  # sort for easier comparison
    @test octind3 == [2, 7, 8, 10, 12, 13, 14, 15]

    @test octree2.particle_indexes_sorted[11:15] == expected1115
    
    # 1 particle per bin
    # first 3 octants/suboctants remain in place
    for i in 1:3
        @test octree2.bin_start[i] == i
    end
    # the old octant pointers skip the particles in suboctant 3 (1 particle in 3, 7 in next suboctants)
    for i in 4:8
        @test octree2.bin_start[i] == i + 7
    end
    # the new suboctants are tacked on to the end
    # so the indices point backwards to where the particles in octant 3 were/are
    for i in 9:15
        @test octree2.bin_start[i] == i - 5
    end

    num_per_bin = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    w_per_bin = [1, 1, 0.5, 1, 1, 1, 1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
    for i in 1:15
        @test octree2.bins[i].np == num_per_bin[i]
        @test octree2.bins[i].w == w_per_bin[i]
    end
end