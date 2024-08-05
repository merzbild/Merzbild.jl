@testset "grid sorting" begin
    Δlarge = 24.0

    seed = 1234
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    # domain length of 8, 2 cells
    grid_coarse = create_grid1D_uniform(8, 2)

    n_per_cell = 1e10
    ppc = 4
    Fnum::Float64 = n_per_cell / ppc
    T = 500.0

    particles = [create_particle_vector(ppc * grid_coarse.n_cells)]

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_list(particles_data_path, "Ar")

    pia_coarse = create_particle_indexer_array(grid_coarse.n_cells, 1)

    sample_particles_equal_weight!(rng, grid_coarse, particles[1], pia_coarse, 1,
                                   species_data, ppc, T, Fnum)

    # reverse particles
    # x = [7.92 6.93 5.94 4.95] [3.96 2.97 1.98 0.99]
    for i in 1:pia_coarse.n_total[1]
        particles[1][i].x = SVector{3,Float64}(0.99 * 8.0 * (9.0 - i) / 8, 0.0, 0.0)
    end

    gridsorter = create_grid_sort_inplace(grid_coarse, pia_coarse.n_total[1])
    sort_particles!(gridsorter, grid_coarse, particles[1], pia_coarse, 1)
    
    # x (sorted) = [3.96 2.97 1.98 0.99] [7.92 6.93 5.94 4.95] 

    # test that we correctly set stuff in the particle indexer
    @test pia_coarse.indexer[1,1].n_group1 == 4
    @test pia_coarse.indexer[1,1].start1 == 1
    @test pia_coarse.indexer[1,1].end1 == 4
    @test pia_coarse.indexer[1,1].n_group2 == 0
    @test pia_coarse.indexer[1,1].start2 == 0
    @test pia_coarse.indexer[1,1].end2 == -1
    @test pia_coarse.indexer[1,1].n_local == 4

    @test pia_coarse.indexer[2,1].n_group1 == 4
    @test pia_coarse.indexer[2,1].start1 == 5
    @test pia_coarse.indexer[2,1].end1 == 8
    @test pia_coarse.indexer[2,1].n_group2 == 0
    @test pia_coarse.indexer[2,1].start2 == 0
    @test pia_coarse.indexer[2,1].end2 == -1
    @test pia_coarse.indexer[2,1].n_local == 4

    @test particles[1].index == [5, 6, 7, 8, 1, 2, 3, 4]

    # underlying particles did not change
    for i in 1:8
        particles[1].particles[i].x[1] == 0.99 * 8.0 * (9.0 - i) / 8
    end

    # but when we access them through the ParticleVector we get the expected result
    for i in 1:4
        @test particles[1][i].x[1] == 0.99 * 8.0 * (9.0 - i - 4) / 8
    end
    for i in 5:8
        @test particles[1][i].x[1] == 0.99 * 8.0 * (9.0 - i + 4) / 8
    end

    # now we try a finer grid
    grid_fine = create_grid1D_uniform(8, 4)

    pia_fine = create_particle_indexer_array(grid_fine.n_cells, 1)

    # we don't need to fill in per-cell info since it will be overwritten anyway
    pia_fine.n_total[1] = 8

    gridsorter_fine = create_grid_sort_inplace(grid_fine, pia_fine.n_total[1])
    sort_particles!(gridsorter_fine, grid_fine, particles[1], pia_fine, 1)

    # x (sorted) = [1.98 0.99] [3.96 2.97] [5.94 4.95] [7.92 6.93]
    @test particles[1].index == [7, 8, 5, 6, 3, 4, 1, 2]

    for cell in 1:4
        @test pia_fine.indexer[cell,1].n_group1 == 2
        @test pia_fine.indexer[cell,1].start1 == 1 + 2 * (cell - 1)
        @test pia_fine.indexer[cell,1].end1 == 2 + 2 * (cell - 1)
        @test pia_fine.indexer[cell,1].n_group2 == 0
        @test pia_fine.indexer[cell,1].start2 == 0
        @test pia_fine.indexer[cell,1].end2 == -1
        @test pia_fine.indexer[cell,1].n_local == 2
    end

    # underlying particles did not change
    for i in 1:8
        particles[1].particles[i].x[1] == 0.99 * 8.0 * (9.0 - i) / 8
    end

    for i in 1:2
        @test particles[1][i].x[1] == 0.99 * 8.0 * (9.0 - i - 6) / 8
    end
    for i in 3:4
        @test particles[1][i].x[1] == 0.99 * 8.0 * (9.0 - i - 2) / 8
    end
    for i in 5:6
        @test particles[1][i].x[1] == 0.99 * 8.0 * (9.0 - i + 2) / 8
    end
    for i in 7:8
        @test particles[1][i].x[1] == 0.99 * 8.0 * (9.0 - i + 6) / 8
    end

    # now we split the particles unevenly across 3 cells
    # [0.0,2.0] - 3 particles
    # [2.0, 4.0] - 1 particle
    # [4.0, 6.0] - 0 particles
    # [6.0, 8.0] - 4 particles

    # we reset the index
    particles[1].index = [1, 2, 3, 4, 5, 6, 7, 8]

    for i in 1:4
        particles[1][i].x = SVector{3,Float64}(6.75, 0.0, 0.0)
    end
    for i in 5:5
        particles[1][i].x = SVector{3,Float64}(2.5, 0.0, 0.0)
    end
    for i in 6:8
        particles[1][i].x = SVector{3,Float64}(0.5, 0.0, 0.0)
    end

    sort_particles!(gridsorter_fine, grid_fine, particles[1], pia_fine, 1)
    @test particles[1].index == [6, 7, 8, 5, 1, 2, 3, 4]
    counts = [3, 1, 0, 4]
    starts = [1, 4, 0, 5]
    ends = [3, 4, -1, 8]

    for cell in 1:4
        @test pia_fine.indexer[cell,1].n_group1 == counts[cell]
        @test pia_fine.indexer[cell,1].start1 == starts[cell]
        @test pia_fine.indexer[cell,1].end1 == ends[cell]
        @test pia_fine.indexer[cell,1].n_group2 == 0
        @test pia_fine.indexer[cell,1].start2 == 0
        @test pia_fine.indexer[cell,1].end2 == -1
        @test pia_fine.indexer[cell,1].n_local == counts[cell]
    end

    phys_props::PhysProps = create_props(grid_fine.n_cells, 1, [], Tref=1)
    compute_props!(particles, pia_fine, species_data, phys_props)
    n_per_cell = Fnum * counts

    for i in 1:grid_fine.n_cells
        @test abs(phys_props.n[i, 1] - n_per_cell[i]) < 2*eps()
        @test phys_props.np[i, 1] == counts[i]
    end
end 