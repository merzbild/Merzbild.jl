@testset "conservative N:1 octree_merging, 3D particles" begin

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

        return Particle(Float64(w), [v_x, v_y, v_z], [1.0, -10.0, 3.0])
    end
    
    function create_24_3particles_in_octant(; weights=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    # create 3 particles in octant, each with weight == octant * weights[octant]
    # and velocity = 9.0 - octant - 0.5 / 9.0 - octant + 0.5 / 9.0 - octant + 0.01 * i
    # so that we avoid having same particles (otherwise refinement acts weird because we have identical particles and it tries to refine to bins with 1 particle)
        vp = ParticleVector(24)

        i = 0
        for octant in 1:8
            i += 1
            Merzbild.update_particle_buffer_new_particle!(vp, i)
            vp[i] = create_particle_in_octant(octant, 9.0 - octant - 0.5 + 0.01 * i, w=octant*weights[octant])
            i += 1
            Merzbild.update_particle_buffer_new_particle!(vp, i)
            vp[i] = create_particle_in_octant(octant, 9.0 - octant + 0.5 + 0.01 * i, w=octant*weights[octant])
            i += 1
            Merzbild.update_particle_buffer_new_particle!(vp, i)
            vp[i] = create_particle_in_octant(octant, 9.0 - octant + 0.01 * i, w=octant*weights[octant])
        end

        return vp
    end

    function create_2particles_total()
    # create just 2 particles
        vp = ParticleVector(2)

        i = 0
        i += 1
        Merzbild.update_particle_buffer_new_particle!(vp, i)
        vp[i] = create_particle_in_octant(1, 9.0 - 1 - 0.5, w=1)
        i += 1
        Merzbild.update_particle_buffer_new_particle!(vp, i)
        vp[i] = create_particle_in_octant(1, 9.0 - 1 + 0.5, w=2)

        return vp
    end

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, "Ar")

    seed = 1234
    rng = StableRNG(seed)

    phys_props::PhysProps = PhysProps(1, 1)
    
    # particles24::Vector{Vector{Particle}} = [create_24_3particles_in_octant()]
    particles24 = [create_24_3particles_in_octant()]
    pia = ParticleIndexerArray(24)

    # symmetric octree with split at v0 = (0.0, 0.0, 0.0)
    octree = OctreeMerge{3,1}(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel)
    
    total_w = sum([3 * i for i in 1:8])
    
    for i in 1:8
        Merzbild.compute_bin_props!(octree, i, particles24[1])
    end

    compute_props!(particles24, pia, species_data, phys_props)

    n0_computed = phys_props.n[1,1]
    np0_computed = phys_props.np[1,1]
    T0_computed = phys_props.T[1,1]
    v0_computed = phys_props.v[:,1,1]
    @test n0_computed == total_w

    octree2 = OctreeMerge{3,1}(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel)
    merge_octree!(rng, octree2, particles24[1], pia, 1, 1, 16)

    @test octree2.Nbins == 10
    @test pia.n_total[1] == 10
    @test pia.n_total[1] == pia.indexer[1,1].n_local

    compute_props!(particles24, pia, species_data, phys_props)
    @test pia.n_total[1] == phys_props.np[1,1]
    @test abs(phys_props.n[1,1] - n0_computed) < eps()
    @test abs(phys_props.T[1,1] - T0_computed) < 1e-14
    @test abs(v0_computed[1] - phys_props.v[1,1,1]) < 1e-14
    @test abs(v0_computed[2] - phys_props.v[2,1,1]) < 1e-14
    @test abs(v0_computed[3] - phys_props.v[3,1,1]) < 1e-14

    merge_octree!(rng,octree2, particles24[1], pia, 1, 1, 2)
    # no conservation possible here
    @test octree2.Nbins == 1
    @test pia.n_total[1] == 1
    for i in 1:1
        @test octree2.bins[i].depth == 0
    end

    compute_props!(particles24, pia, species_data, phys_props)
    @test pia.n_total[1] == phys_props.np[1,1]
    @test particles24[1][1].w == total_w
    @test abs(phys_props.n[1,1] - n0_computed) < eps()
    @test abs(v0_computed[1] - phys_props.v[1,1,1]) < 1e-14
    @test abs(v0_computed[2] - phys_props.v[2,1,1]) < 1e-14
    @test abs(v0_computed[3] - phys_props.v[3,1,1]) < 1e-14

    weights_arr = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0]
    particles24 = [create_24_3particles_in_octant(weights=weights_arr)]
    pia = ParticleIndexerArray(24)
    compute_props!(particles24, pia, species_data, phys_props)

    n0_computed = phys_props.n[1,1]
    np0_computed = phys_props.np[1,1]
    T0_computed = phys_props.T[1,1]
    v0_computed = phys_props.v[:,1,1]

    merge_octree!(rng, octree2, particles24[1], pia, 1, 1, 8)
    compute_props!(particles24, pia, species_data, phys_props)
    @test pia.n_total[1] == sum(weights_arr .> 0.0) # we should skip bins with weight 0.0
    @test pia.n_total[1] == phys_props.np[1,1]
    @test abs(phys_props.n[1,1] - n0_computed) < eps()
    @test abs(phys_props.T[1,1] - T0_computed) < 1e-14
    @test abs(v0_computed[1] - phys_props.v[1,1,1]) < 1e-14
    @test abs(v0_computed[2] - phys_props.v[2,1,1]) < 1e-14
    @test abs(v0_computed[3] - phys_props.v[3,1,1]) < 1e-14

    s1 = pia.indexer[1,1].start1
    e1 = pia.indexer[1,1].end1

    @test e1 == pia.n_total[1]

    for i in s1:e1
        @test isnan(particles24[1][i].x[1]) == false
        @test particles24[1][i].w > 0
    end
    s2 = pia.indexer[1,1].start2
    e2 = pia.indexer[1,1].end2

    @test e2 == -1
    @test s2 == 0

    for i in s2:e2
        @test isnan(particles24[1][i].x[1]) == false
        @test particles24[1][i].w > 0
    end

    # now we split into 2 indexing groups
    weights_arr = [1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0]
    particles24 = [create_24_3particles_in_octant(weights=weights_arr)]
    pia = ParticleIndexerArray(24)

    pia.indexer[1,1].n_group1 = 7
    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = 7

    pia.indexer[1,1].n_group2 = 17
    pia.indexer[1,1].start2 = 8
    pia.indexer[1,1].end2 = 24

    n_zero = 0

    s1 = pia.indexer[1,1].start1
    e1 = pia.indexer[1,1].end1
    s2 = pia.indexer[1,1].start2
    e2 = pia.indexer[1,1].end2

    for i in s1:e1
        if particles24[1][i].w == 0
            n_zero += 1
        end
    end
    for i in s2:e2
        if particles24[1][i].w == 0
            n_zero += 1
        end
    end
    @test n_zero == sum(weights_arr .== 0.0) * 3

    compute_props!(particles24, pia, species_data, phys_props)

    n0_computed = phys_props.n[1,1]
    np0_computed = phys_props.np[1,1]
    T0_computed = phys_props.T[1,1]
    v0_computed = phys_props.v[:,1,1]

    merge_octree!(rng, octree2, particles24[1], pia, 1, 1, 8)
    compute_props!(particles24, pia, species_data, phys_props)
    @test pia.n_total[1] == sum(weights_arr .> 0.0) # we should skip bins with weight 0.0
    @test pia.n_total[1] == phys_props.np[1,1]
    @test abs(phys_props.n[1,1] - n0_computed) < eps()
    @test abs(phys_props.T[1,1] - T0_computed) < 1e-14
    @test abs(v0_computed[1] - phys_props.v[1,1,1]) < 1e-14
    @test abs(v0_computed[2] - phys_props.v[2,1,1]) < 1e-14
    @test abs(v0_computed[3] - phys_props.v[3,1,1]) < 1e-14

    s1 = pia.indexer[1,1].start1
    e1 = pia.indexer[1,1].end1
    s2 = pia.indexer[1,1].start2
    e2 = pia.indexer[1,1].end2

    @test e1 - s1 + 1 + e2 - s2 + 1 == pia.n_total[1]

    for i in s1:e1
        @test isnan(particles24[1][i].x[1]) == false
        @test particles24[1][i].w > 0
    end

    @test e2 == -1
    @test s2 == 0

    particles2 = [create_2particles_total()]
    pia = ParticleIndexerArray(2)

    pia.indexer[1,1].n_local = 2
    pia.indexer[1,1].n_group1 = 2
    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = 2

    pia.indexer[1,1].n_group2 = 0
    pia.indexer[1,1].start2 = 0
    pia.indexer[1,1].end2 = -1

    # even the top-level bin cannot be refined
    merge_octree!(rng, octree2, particles2[1], pia, 1, 1, 16)
    @test octree2.bins[1].np == 1
    @test octree2.n_particles == 2
    @test octree2.bins[1].can_be_refined == false

    # non-inheriting bin bounds version
    particles24 = [create_24_3particles_in_octant()]
    pia = ParticleIndexerArray(24)

    # symmetric octree with split at v0 = (0.0, 0.0, 0.0)
    octree = OctreeMerge{3,1}(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel, bin_bounds_compute=OctreeBinBoundsInherit)
    merge_octree!(rng, octree, particles24[1], pia, 1, 1, 16)
    compute_props!(particles24, pia, species_data, phys_props)

    n0_computed = phys_props.n[1,1]
    np0_computed = phys_props.np[1,1]
    T0_computed = phys_props.T[1,1]
    v0_computed = phys_props.v[:,1,1]

    compute_props!(particles24, pia, species_data, phys_props)

    @test pia.n_total[1] == phys_props.np[1,1]
    @test abs(phys_props.n[1,1] - n0_computed) < eps()
    @test abs(phys_props.T[1,1] - T0_computed) < 1e-14
    @test abs(v0_computed[1] - phys_props.v[1,1,1]) < 1e-14
    @test abs(v0_computed[2] - phys_props.v[2,1,1]) < 1e-14
    @test abs(v0_computed[3] - phys_props.v[3,1,1]) < 1e-14
end