@testset "octree_merging_buffer_sorting" begin
    # test that buffers and indices work correctly in octree merging
    # and that sorting restores proper indexing in pia
    function create_particles_and_pia(n_gr1)
        # n_gr1 particles in cell 1
        # 10 particles in cell 2, particles in cell 1 are split into 2 groups: [1,n_gr1] and [n_gr1+11, 100]
        pv = ParticleVector(100)
        for i in 1:100
            Merzbild.update_particle_buffer_new_particle(pv, i)

            if i >= n_gr1+1 && i <= n_gr1+10
                w = 1000.0
                x_pos = 5.0
            else
                w = 1.0
                x_pos = 1.0
            end
            pv[i] = Particle(w, [i - 1.0, (i^2) + 3.0, sqrt(i + 5.0)^(2 * (i%2) - 1.0)], [x_pos, 0.0, 0.0])
        end
        # 2 cells, 1 species
        pia = ParticleIndexerArray(2,1)

        pia.n_total[1] = 100
        pia.indexer[1,1].n_local = 90
        pia.indexer[1,1].start1 = 1
        pia.indexer[1,1].end1 = n_gr1
        pia.indexer[1,1].n_group1 = n_gr1

        pia.indexer[1,1].start2 = n_gr1+11
        pia.indexer[1,1].end2 = 100
        pia.indexer[1,1].n_group2 = 90-n_gr1

        pia.indexer[2,1].n_local = 10
        pia.indexer[2,1].start1 = n_gr1+1
        pia.indexer[2,1].end1 = n_gr1+10
        pia.indexer[2,1].n_group1 = 10

        pia.indexer[2,1].start2 = -1
        pia.indexer[2,1].end2 = -1
        pia.indexer[2,1].n_group2 = 0

        return [pv], pia
    end

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")

    seed = 1234
    rng = StableRNG(seed)

    phys_props = PhysProps(2, 1, [], Tref=1)

    particles, pia = create_particles_and_pia(50)
    
    compute_props!(particles, pia, species_data, phys_props)

    @test phys_props.np[1, 1] == 90
    @test phys_props.np[2, 1] == 10

    @test phys_props.n[1, 1] == 90 * 1.0
    @test phys_props.n[2, 1] == 10 * 1000.0

    @test particles[1].nbuffer == 0

    @test pia.contiguous[1] == true

    octree = OctreeN2Merge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel)
    merge_octree_N2_based!(rng, octree, particles[1], pia, 1, 1, 6)

    compute_props!(particles, pia, species_data, phys_props)

    @test pia.contiguous[1] == false

    @test phys_props.np[1, 1] == 2
    @test phys_props.np[2, 1] == 10

    @test phys_props.n[1, 1] == 90 * 1.0
    @test phys_props.n[2, 1] == 10 * 1000.0

    @test particles[1].nbuffer == 88

    # now we test what happend to the buffer
    buffer_part1 = true
    for i in 1:40
        if particles[1].buffer[i] != 100-i+1
            buffer_part1 = false
        end
    end
    @test buffer_part1 == true

    buffer_part2 = true
    for i in 41:88
        if particles[1].buffer[i] != 50-i+41
            buffer_part1 = false
        end
    end
    @test buffer_part2 == true

    @test pia.indexer[1,1].n_local == 2
    @test pia.indexer[1,1].n_group1 == 2
    @test pia.indexer[1,1].start1 == 1
    @test pia.indexer[1,1].end1 == 2
    @test pia.indexer[1,1].n_group2 <= 0
    @test pia.indexer[1,1].start2 <= 0

    @test pia.indexer[2,1].n_local == 10
    @test pia.indexer[2,1].n_group1 == 10
    @test pia.indexer[2,1].start1 == 51
    @test pia.indexer[2,1].end1 == 60
    @test pia.indexer[2,1].n_group2 <= 0
    @test pia.indexer[2,1].start2 <= 0

    # finally, we check that if we sort particles, the pia indexing becomes contiguous again
    # currently, the end1 + 1 in cell1 != start1 in cell2, so there's a "hole" in the indexing array

    # 2 cells, [0,4], [4,8]
    grid_coarse = Grid1DUniform(8.0, 2)
    
    gridsorter = GridSortInPlace(grid_coarse, 100)
    sort_particles!(gridsorter, grid_coarse, particles[1], pia, 1)

    compute_props!(particles, pia, species_data, phys_props)

    @test pia.contiguous[1] == true

    @test phys_props.np[1, 1] == 2
    @test phys_props.np[2, 1] == 10

    @test phys_props.n[1, 1] == 90 * 1.0
    @test phys_props.n[2, 1] == 10 * 1000.0

    @test pia.indexer[1,1].n_local == 2
    @test pia.indexer[1,1].n_group1 == 2
    @test pia.indexer[1,1].start1 == 1
    @test pia.indexer[1,1].end1 == 2
    @test pia.indexer[1,1].n_group2 == 0
    @test pia.indexer[1,1].start2 == 0

    @test pia.indexer[2,1].n_local == 10
    @test pia.indexer[2,1].n_group1 == 10
    @test pia.indexer[2,1].start1 == 3
    @test pia.indexer[2,1].end1 == 12
    @test pia.indexer[2,1].n_group2 == 0
    @test pia.indexer[2,1].start2 == 0


    particles, pia = create_particles_and_pia(5)
    
    compute_props!(particles, pia, species_data, phys_props)

    @test phys_props.np[1, 1] == 90
    @test phys_props.np[2, 1] == 10

    @test phys_props.n[1, 1] == 90 * 1.0
    @test phys_props.n[2, 1] == 10 * 1000.0


    @test particles[1].nbuffer == 0

    @test pia.contiguous[1] == true

    octree = OctreeN2Merge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel)

    # 16 is the target particle number but we will get 12 post-merge particles
    merge_octree_N2_based!(rng, octree, particles[1], pia, 1, 1, 16)

    compute_props!(particles, pia, species_data, phys_props)

    @test pia.contiguous[1] == false

    @test phys_props.np[1, 1] == 12
    @test phys_props.np[2, 1] == 10

    @test phys_props.n[1, 1] == 90 * 1.0
    @test phys_props.n[2, 1] == 10 * 1000.0

    @test particles[1].nbuffer == 78

    # now we test what happend to the buffer
    buffer_part1 = true
    for i in 1:78
        if particles[1].buffer[i] != 100-i+1
            buffer_part1 = false
        end
    end
    @test buffer_part1 == true

    @test pia.indexer[1,1].n_local == 12
    @test pia.indexer[1,1].n_group1 == 5
    @test pia.indexer[1,1].start1 == 1
    @test pia.indexer[1,1].end1 == 5
    @test pia.indexer[1,1].n_group2 == 7
    @test pia.indexer[1,1].start2 == 16
    @test pia.indexer[1,1].end2 == 22

    @test pia.indexer[2,1].n_local == 10
    @test pia.indexer[2,1].n_group1 == 10
    @test pia.indexer[2,1].start1 == 6
    @test pia.indexer[2,1].end1 == 15
    @test pia.indexer[2,1].n_group2 <= 0
    @test pia.indexer[2,1].start2 <= 0
end