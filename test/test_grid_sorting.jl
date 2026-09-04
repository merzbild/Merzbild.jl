@testset "grid sorting" begin
    Δlarge = 24.0

    seed = 1234
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    # domain length of 8, 2 cells
    grid_coarse = Grid1DUniform(8, 2)

    n_per_cell = 1e10
    ppc = 4
    Fnum::Float64 = n_per_cell / ppc
    T = 500.0

    particles = [ParticleVector(ppc * grid_coarse.n_cells)]

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, "Ar")

    pia_coarse = ParticleIndexerArray(grid_coarse.n_cells, 1)

    sample_particles_equal_weight!(rng, grid_coarse, particles[1], pia_coarse, 1,
                                   species_data, ppc, T, Fnum)

    # reverse particles
    # x = [7.92 6.93 5.94 4.95] [3.96 2.97 1.98 0.99]
    for i in 1:pia_coarse.n_total[1]
        particles[1][i].x = SVector{3,Float64}(0.99 * 8.0 * (9.0 - i) / 8, 0.0, 0.0)
    end

    gridsorter = GridSortInPlace(grid_coarse, pia_coarse.n_total[1])
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

    # test resizing
    particles2 = [ParticleVector(ppc * grid_coarse.n_cells)]
    pia_coarse2 = ParticleIndexerArray(grid_coarse.n_cells, 1)
    sample_particles_equal_weight!(rng, grid_coarse, particles2[1], pia_coarse2, 1,
                                   species_data, ppc, T, Fnum)
    for i in 1:pia_coarse.n_total[1]
        particles2[1][i].x = SVector{3,Float64}(0.99 * 8.0 * (9.0 - i) / 8, 0.0, 0.0)
    end

    gridsorter_small = GridSortInPlace(grid_coarse, 1)
    @test length(gridsorter_small.sorted_indices) < pia_coarse2.n_total[1]
    sort_particles!(gridsorter_small, grid_coarse, particles2[1], pia_coarse2, 1)
    @test length(gridsorter_small.sorted_indices) > pia_coarse2.n_total[1]

    @test particles2[1].index == [5, 6, 7, 8, 1, 2, 3, 4]

    # underlying particles did not change
    for i in 1:8
        particles2[1].particles[i].x[1] == 0.99 * 8.0 * (9.0 - i) / 8
    end

    # but when we access them through the ParticleVector we get the expected result
    for i in 1:4
        @test particles2[1][i].x[1] == 0.99 * 8.0 * (9.0 - i - 4) / 8
    end
    for i in 5:8
        @test particles2[1][i].x[1] == 0.99 * 8.0 * (9.0 - i + 4) / 8
    end

    # now we try a finer grid
    grid_fine = Grid1DUniform(8, 4)

    pia_fine = ParticleIndexerArray(grid_fine.n_cells, 1)

    # we just need to set indexing so that we can iterate through the particles,
    # exact cell counts don't matter: we're sorting anyway
    pia_fine.n_total[1] = 8
    pia_fine.indexer[1,1].n_local = 8
    pia_fine.indexer[1,1].start1 = 1
    pia_fine.indexer[1,1].end1 = 8
    pia_fine.indexer[1,1].n_group1 = 8

    gridsorter_fine = GridSortInPlace(grid_fine, pia_fine.n_total[1])
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
    @test gridsorter_fine.occ_lo == 1
    @test gridsorter_fine.occ_hi == 4
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

    phys_props = PhysProps(grid_fine.n_cells, 1)
    compute_props!(particles, pia_fine, species_data, phys_props)
    n_per_cell = Fnum * counts

    for i in 1:grid_fine.n_cells
        @test abs(phys_props.n[i, 1] - n_per_cell[i]) < 2*eps()
        @test phys_props.np[i, 1] == counts[i]
    end

    # repeat the test case above
    # but now we use the sorting routine which assumes particle cells have already been set during convection
    # so the coordinates of the particles play no role whatsoever
    # now we split the particles unevenly across 3 cells
    # [0.0,2.0] - 3 particles (6, 7, 8)
    # [2.0, 4.0] - 1 particle (5)
    # [4.0, 6.0] - 0 particles ()
    # [6.0, 8.0] - 4 particles (1, 2, 3, 4)

    # we reset the index
    particles[1].index = [1, 2, 3, 4, 5, 6, 7, 8]
    particles[1].cell = [4, 4, 4, 4, 2, 1, 1, 1]

    for i in 1:4
        particles[1][i].x = SVector{3,Float64}(0.05, 0.0, 0.0)
    end
    for i in 5:5
        particles[1][i].x = SVector{3,Float64}(0.1, 0.0, 0.0)
    end
    for i in 6:8
        particles[1][i].x = SVector{3,Float64}(0.3, 0.0, 0.0)
    end

    sort_particles!(gridsorter_fine, particles[1], pia_fine, 1)
    @test particles[1].index == [6, 7, 8, 5, 1, 2, 3, 4]
    @test gridsorter_fine.occ_lo == 1
    @test gridsorter_fine.occ_hi == 4
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

    phys_props = PhysProps(grid_fine.n_cells, 1)
    compute_props!(particles, pia_fine, species_data, phys_props)
    n_per_cell = Fnum * counts

    for i in 1:grid_fine.n_cells
        @test abs(phys_props.n[i, 1] - n_per_cell[i]) < 2*eps()
        @test phys_props.np[i, 1] == counts[i]
    end
end 
@testset "grid sorting occupancy tracking" begin
    n_cells = 6
    np = 8

    particles = [ParticleVector(np)]
    pia = ParticleIndexerArray(n_cells, 1)

    # the sorting routine taking pre-computed cells is used, so only the cell of each
    # particle matters here, not its position. The per-cell indexing is left exactly as the
    # previous sort produced it, which is the state convection leaves the pia in
    function sort_in_cells!(gridsorter, cells)
        particles[1].index = collect(1:np)
        particles[1].cell = copy(cells)

        pia.n_total[1] = length(cells)
        pia.index_last[1] = length(cells)
        pia.contiguous[1] = true

        sort_particles!(gridsorter, particles[1], pia, 1)
    end

    # the cells in which the pia actually indexes particles
    function occupied(pia)
        return [cell for cell in 1:pia.n_cells if pia.indexer[cell,1].n_group1 > 0]
    end

    # every cell not holding particles has to be indexed as empty
    function test_empty_cells(pia)
        for cell in 1:pia.n_cells
            if pia.indexer[cell,1].n_group1 == 0
                @test pia.indexer[cell,1].start1 == 0
                @test pia.indexer[cell,1].end1 == -1
                @test pia.indexer[cell,1].n_local == 0
                @test pia.indexer[cell,1].start2 == 0
                @test pia.indexer[cell,1].end2 == -1
                @test pia.indexer[cell,1].n_group2 == 0
            end
        end
    end

    gridsorter = GridSortInPlace(n_cells, np)
    @test gridsorter.occ_lo == 1
    @test gridsorter.occ_hi == n_cells

    sort_in_cells!(gridsorter, [2, 2, 3, 3, 4, 4, 5, 5])
    @test gridsorter.occ_lo == 2
    @test gridsorter.occ_hi == 5
    @test occupied(pia) == [2, 3, 4, 5]
    test_empty_cells(pia)

    # the occupancy range shrinks to a single cell: the indexing of the cells occupied
    # during the previous sort has to be cleared even though they are outside of the new range
    sort_in_cells!(gridsorter, [3, 3, 3, 3, 3, 3, 3, 3])
    @test gridsorter.occ_lo == 3
    @test gridsorter.occ_hi == 3
    @test occupied(pia) == [3]
    @test pia.indexer[3,1].start1 == 1
    @test pia.indexer[3,1].end1 == 8
    test_empty_cells(pia)

    # and now it widens in both directions, with an empty cell inside the range
    sort_in_cells!(gridsorter, [1, 1, 2, 2, 4, 4, 6, 6])
    @test gridsorter.occ_lo == 1
    @test gridsorter.occ_hi == 6
    @test occupied(pia) == [1, 2, 4, 6]
    test_empty_cells(pia)

    # sorting an empty pia clears everything and records an empty range
    sort_in_cells!(gridsorter, Int64[])
    @test gridsorter.occ_lo == n_cells + 1
    @test gridsorter.occ_hi == 0
    @test occupied(pia) == []
    test_empty_cells(pia)

    # sorting an empty pia twice in a row is a no-op
    sort_in_cells!(gridsorter, Int64[])
    @test gridsorter.occ_lo == n_cells + 1
    @test gridsorter.occ_hi == 0
    test_empty_cells(pia)

    # and particles can be sorted again afterwards
    sort_in_cells!(gridsorter, [5, 5, 5, 5, 6, 6, 6, 6])
    @test gridsorter.occ_lo == 5
    @test gridsorter.occ_hi == 6
    @test occupied(pia) == [5, 6]
    test_empty_cells(pia)

    # indexing written by other means than the sorting routines is cleared as well
    stale = pia.indexer[1,1]
    stale.n_local = 1
    stale.n_group1 = 1
    stale.start1 = 1
    stale.end1 = 1

    sort_in_cells!(gridsorter, [4, 4, 4, 4, 4, 4, 4, 4])
    @test gridsorter.occ_lo == 4
    @test gridsorter.occ_hi == 4
    @test occupied(pia) == [4]
    test_empty_cells(pia)
end
