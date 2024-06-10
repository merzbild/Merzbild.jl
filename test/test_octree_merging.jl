@testset "octree_merging" begin

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

        return Particle(w, [v_x, v_y, v_z], [1.0, -10.0, 3.0])
    end
    
    function create_24_3particles_in_octant()
    # create 3 particles in octant, each with weight == octant
    # and velocity = 9.0 - octant - 0.5 / 9.0 - octant + 0.5 / 9.0 - octant
        vp = Vector{Particle}(undef, 24)

        i = 0
        for octant in 1:8
            i += 1
            vp[i] = create_particle_in_octant(octant, 9.0 - octant - 0.5, w=octant * 1.0)
            i += 1
            vp[i] = create_particle_in_octant(octant, 9.0 - octant + 0.5, w=octant * 1.0)
            i += 1
            vp[i] = create_particle_in_octant(octant, 9.0 - octant, w=octant * 1.0)
        end

        return vp
    end


    species_list::Vector{Species} = load_species_list("data/particles.toml", "Ar")

    seed = 1234
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)
    
    particles24::Vector{Vector{Particle}} = [create_24_3particles_in_octant()]
    pia = create_particle_indexer_array(24)

    # symmetric octree with split at v0 = (0.0, 0.0, 0.0)
    octree = create_merging_octree(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinC)

    Merzbild.init_octree!(1, 1, octree, particles24[1], pia)
    
    Merzbild.split_bin!(octree, 1, particles24[1])
    @test octree.Nbins == 8
    
    for i in 1:8
        @test octree.bins[i].np == 3
        @test octree.bins[i].w == 3 * i
    end

    for i in 1:8
        Merzbild.compute_bin_props!(octree, i, particles24[1])
    end

    for i in 1:8

        # compute signs
        if i >= 5
            v_zs = 1
        else
            v_zs = -1
        end
        if i % 2 == 1
            v_xs = -1
        else
            v_xs = 1
        end
        if (i == 3) || (i == 4) || (i == 7) || (i == 8)
            v_ys = 1
        else
            v_ys = -1
        end

        @test abs(octree.full_bins[i].v_mean[1] - (9.0 - i) * v_xs) < 1e-14 
        @test abs(octree.full_bins[i].v_mean[2] - (9.0 - i) * v_ys) < 1e-14 
        @test abs(octree.full_bins[i].v_mean[3] - (9.0 - i) * v_zs) < 1e-14 

        @test abs(octree.full_bins[i].x_mean[1] - (1.0)) < 1e-14
        @test abs(octree.full_bins[i].x_mean[2] - (-10.0)) < 1e-14
        @test abs(octree.full_bins[i].x_mean[3] - (3.0)) < 1e-14
    end
end