@testset "octree_merging 1D" begin
    
    function create_particles_in_2cells()
    # create 4 particles in 2 cells
    # such that mean(x)+-std(x) is outside of the domain
        vp = ParticleVector(8)

        x_cell1 = [0.05, 0.05, 0.05, 0.45]
        for i in 1:4
            Merzbild.update_particle_buffer_new_particle(vp, i)
            vp[i] = Particle(2.0, [0.5, -3.0, 4.0], [x_cell1[i], 0.0, 0.0])
        end

        x_cell2 = [0.55, 0.95, 0.85, 0.99]
        for i in 5:8
            Merzbild.update_particle_buffer_new_particle(vp, i)
            vp[i] = Particle(3.0, [0.5, -3.0, 4.0], [x_cell2[i-4], 0.0, 0.0])
        end

        pia_ = ParticleIndexerArray(grid.n_cells, 1)

        pia_.n_total[1] = 8

        pia_.indexer[1,1].n_local = 4
        pia_.indexer[1,1].n_group1 = 4
        pia_.indexer[1,1].start1 = 1
        pia_.indexer[1,1].end1 = 4

        pia_.indexer[1,1].n_group2 = 0
        pia_.indexer[1,1].start2 = 0
        pia_.indexer[1,1].end2 = 0

        pia_.indexer[2,1].n_local = 4
        pia_.indexer[2,1].n_group1 = 4
        pia_.indexer[2,1].start1 = 5
        pia_.indexer[2,1].end1 = 8

        pia_.indexer[2,1].n_group2 = 0
        pia_.indexer[2,1].start2 = 0
        pia_.indexer[2,1].end2 = 0

        return [vp], pia_
    end

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, "Ar")

    grid = Grid1DUniform(1.0, 2)

    seed = 1234
    rng = StableRNG(seed)

    phys_props::PhysProps = PhysProps(grid.n_cells, 1, [], Tref=1)
    
    particles, pia = create_particles_in_2cells()
    
    octree = OctreeN2Merge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinC)
    
    # test that merging without accounting for domain bounds leads to out-of-domain particles
    merge_octree_N2_based!(rng, octree, particles[1], pia, 1, 1, 2)
    merge_octree_N2_based!(rng, octree, particles[1], pia, 2, 1, 2)

    @test pia.indexer[1,1].n_local == 2
    @test pia.indexer[1,1].n_group1 == 2
    @test pia.indexer[1,1].start1 == 1
    @test pia.indexer[1,1].end1 == 2

    @test pia.indexer[2,1].n_local == 2
    @test pia.indexer[2,1].n_group1 == 2
    @test pia.indexer[2,1].start1 == 5
    @test pia.indexer[2,1].end1 == 6

    out_of_bounds = false
    for i in pia.indexer[1,1].start1:pia.indexer[1,1].end1
        if (particles[1][i].x[1] < grid.min_x) || (particles[1][i].x[1] > grid.max_x)
            out_of_bounds = true
        end
    end
    @test out_of_bounds == true


    out_of_bounds = false
    for i in pia.indexer[2,1].start1:pia.indexer[2,1].end1
        if (particles[1][i].x[1] < grid.min_x) || (particles[1][i].x[1] > grid.max_x)
            out_of_bounds = true
        end
    end
    @test out_of_bounds == true

    # reset
    particles, pia = create_particles_in_2cells()
    
    merge_octree_N2_based!(rng, octree, particles[1], pia, 1, 1, 2)
    merge_octree_N2_based!(rng, octree, particles[1], pia, 2, 1, 2)

    @test pia.indexer[1,1].n_local == 2
    @test pia.indexer[1,1].n_group1 == 2
    @test pia.indexer[1,1].start1 == 1
    @test pia.indexer[1,1].end1 == 2

    @test pia.indexer[2,1].n_local == 2
    @test pia.indexer[2,1].n_group1 == 2
    @test pia.indexer[2,1].start1 == 5
    @test pia.indexer[2,1].end1 == 6

    # compute_props!(particles, pia, species_data, phys_props)

    # test that merging without accounting for domain bounds leads to out-of-domain particles
    merge_octree_N2_based!(rng, octree, particles[1], pia, 1, 1, 2, grid)
    merge_octree_N2_based!(rng, octree, particles[1], pia, 2, 1, 2, grid)

    @test pia.indexer[1,1].n_local == 2
    @test pia.indexer[1,1].n_group1 == 2
    @test pia.indexer[1,1].start1 == 1
    @test pia.indexer[1,1].end1 == 2

    @test pia.indexer[2,1].n_local == 2
    @test pia.indexer[2,1].n_group1 == 2
    @test pia.indexer[2,1].start1 == 5
    @test pia.indexer[2,1].end1 == 6

    out_of_bounds = false
    for i in pia.indexer[1,1].start1:pia.indexer[1,1].end1
        if (particles[1][i].x[1] < grid.min_x) || (particles[1][i].x[1] > grid.max_x)
            out_of_bounds = true
        end
    end
    @test out_of_bounds == false

    out_of_bounds = false
    for i in pia.indexer[2,1].start1:pia.indexer[2,1].end1
        if (particles[1][i].x[1] < grid.min_x) || (particles[1][i].x[1] > grid.max_x)
            out_of_bounds = true
        end
    end
    @test out_of_bounds == false

    @test pia.contiguous[1] == false

    gridsorter = GridSortInPlace(grid, pia.n_total[1])
    sort_particles!(gridsorter, grid, particles[1], pia, 1)

    compute_props!(particles, pia, species_data, phys_props)
    @test phys_props.np[1,1] == 2
    @test phys_props.np[2,1] == 2

    @test phys_props.n[1,1] == 8.0
    @test phys_props.n[2,1] == 12.0
end