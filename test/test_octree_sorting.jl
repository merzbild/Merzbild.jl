@testset "octree_sorting" begin

    function create_particle_in_octant(octant, v_val)
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

        return Particle(1.0, [v_x, v_y, v_z], [0.0, 0.0, 0.0])
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
            vp[i] = create_particle_in_octant(octant, 1.0)
            i += 1
            if octant == extra_octant
                vp[i] = create_particle_in_octant(octant, 2.0)
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
                vp[i] = create_particle_in_octant(octant, 1.0)
                i += 1
            end
        end
        
        return vp
    end

    species_list::Vector{Species} = load_species_list("data/particles.toml", "Ar")

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
    pia = create_particle_indexer_array(9)
    octree = create_merging_octree(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinC)

    Merzbild.init_octree!(1, 1, octree, pia)
    for i in 1:9
        @test i in octree.particle_indexes_sorted[1:9]
    end
    
    Merzbild.split_bin!(octree, 1, particles9[1])
    @test maximum(abs.(octree.vel_middle)) < 1e-12  # should be almost 0 since everything is symmetric
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
    Merzbild.init_octree!(1, 1, octree, pia)
    Merzbild.split_bin!(octree, 1, particles9[1])
    @test octree.Nbins == 8
    expected = [9, 8, 7, 6, 1, 2, 3, 4, 5]
    @test octree.particle_indexes_sorted[1:9] == expected

    @test octree.bin_start[1:8] == [1, 2, 3, 5, 6, 7, 8, 9]
    @test octree.bin_end[1:8] == [1, 2, 4, 5, 6, 7, 8, 9]

    # the same but with pia constructed by hand with particle indices split into 2 groups
    pia2 = create_particle_indexer_array(3)
    pia2.n_total[1] = 9
    pia2.indexer[1,1].n_local = 9

    pia2.indexer[1,1].start1 = 1
    pia2.indexer[1,1].end1 = 3
    pia2.indexer[1,1].n_group1 = 3
    
    pia2.indexer[1,1].start2 = 4
    pia2.indexer[1,1].end2 = 9
    pia2.indexer[1,1].n_group2 = 6

    Merzbild.init_octree!(1, 1, octree, pia2)
    Merzbild.split_bin!(octree, 1, particles9[1])

    @test octree.Nbins == 8
    expected = [9, 8, 7, 6, 1, 2, 3, 4, 5]
    @test octree.particle_indexes_sorted[1:9] == expected
    @test octree.bin_start[1:8] == [1, 2, 3, 5, 6, 7, 8, 9]
    @test octree.bin_end[1:8] == [1, 2, 4, 5, 6, 7, 8, 9]

    # now we test for 7 particles, i.e. with one empty octree bin
    # particles:  1, 2, 3, 4, 5, 6, 7
    # octants are 8, 7, 6, 4, 3, 2, 1
    particles7::Vector{Vector{Particle}} = [create_7particles(5, true)]
    pia3 = create_particle_indexer_array(7)

    Merzbild.init_octree!(1, 1, octree, pia3)
    Merzbild.split_bin!(octree, 1, particles7[1])
    @test octree.Nbins == 7
    expected = [7, 6, 5, 4, 3, 2, 1]
    @test octree.particle_indexes_sorted[1:7] == expected
    @test octree.bin_start[1:7] == [1, 2, 3, 4, 5, 6, 7]
    @test octree.bin_end[1:7] == [1, 2, 3, 4, 5, 6, 7]

    # TODO: test double-splitting, i.e. we have 1 p/bin in 7 bins, the 8th bin has 4 particles
    # we do splitting of the main 0 bin, and then do a second split
end