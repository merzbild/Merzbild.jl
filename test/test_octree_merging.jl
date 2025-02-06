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
        vp = ParticleVector(24)

        i = 0
        for octant in 1:8
            i += 1
            Merzbild.update_particle_buffer_new_particle(vp, i)
            vp[i] = create_particle_in_octant(octant, 9.0 - octant - 0.5, w=octant * 1.0)
            i += 1
            Merzbild.update_particle_buffer_new_particle(vp, i)
            vp[i] = create_particle_in_octant(octant, 9.0 - octant + 0.5, w=octant * 1.0)
            i += 1
            Merzbild.update_particle_buffer_new_particle(vp, i)
            vp[i] = create_particle_in_octant(octant, 9.0 - octant, w=octant * 1.0)
        end

        return vp
    end

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, "Ar")

    seed = 1234
    rng = StableRNG(seed)

    phys_props::PhysProps = PhysProps(1, 1, [], Tref=1)
    
    # particles24::Vector{Vector{Particle}} = [create_24_3particles_in_octant()]
    particles24 = [create_24_3particles_in_octant()]
    pia = ParticleIndexerArray(24)

    # symmetric octree with split at v0 = (0.0, 0.0, 0.0)
    octree = OctreeN2Merge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinC)

    Merzbild.init_octree!(octree, particles24[1], pia, 1, 1)
    
    Merzbild.split_bin!(octree, 1, particles24[1])
    @test octree.Nbins == 8
    
    total_w = sum([3 * i for i in 1:8])
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

        @test Merzbild.get_bin_post_merge_np(octree, i) == 2
    end

    compute_props!(particles24, pia, species_data, phys_props)

    n0_computed = phys_props.n[1,1]
    np0_computed = phys_props.np[1,1]
    T0_computed = phys_props.T[1,1]
    v0_computed = phys_props.v[:,1,1]
    @test n0_computed == total_w

    octree2 = OctreeN2Merge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinC)
    merge_octree_N2_based!(rng, octree2, particles24[1], pia, 1, 1, 16)

    @test octree2.Nbins == 8
    @test pia.n_total[1] == 16
    @test pia.n_total[1] == pia.indexer[1,1].n_local
    for i in 1:8
        @test octree2.bins[i].depth == 1
    end

    compute_props!(particles24, pia, species_data, phys_props)
    @test pia.n_total[1] == phys_props.np[1,1]
    @test abs(phys_props.n[1,1] - n0_computed) < eps()
    @test abs(phys_props.T[1,1] - T0_computed) < 1e-14
    @test abs(v0_computed[1] - phys_props.v[1,1,1]) < 1e-14
    @test abs(v0_computed[2] - phys_props.v[2,1,1]) < 1e-14
    @test abs(v0_computed[3] - phys_props.v[3,1,1]) < 1e-14

    merge_octree_N2_based!(rng,octree2, particles24[1], pia, 1, 1, 2)

    @test octree2.Nbins == 1
    @test pia.n_total[1] == 2
    for i in 1:1
        @test octree2.bins[i].depth == 0
    end

    compute_props!(particles24, pia, species_data, phys_props)
    @test pia.n_total[1] == phys_props.np[1,1]
    @test particles24[1][1].w == 0.5 * total_w
    @test particles24[1][1].w == particles24[1][2].w
    @test abs(phys_props.n[1,1] - n0_computed) < eps()
    @test abs(phys_props.T[1,1] - T0_computed) < 1e-14
    @test abs(v0_computed[1] - phys_props.v[1,1,1]) < 1e-14
    @test abs(v0_computed[2] - phys_props.v[2,1,1]) < 1e-14
    @test abs(v0_computed[3] - phys_props.v[3,1,1]) < 1e-14
end