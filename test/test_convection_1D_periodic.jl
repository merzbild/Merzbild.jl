@testset "convection 1D on periodic grid" begin
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")

    seed = 1234
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    grid = Grid1DUniform(50.0, 100)

    gridsorter = GridSortInPlace(grid, 15000)
    particles = [ParticleVector(4)]

    pia = ParticleIndexerArray(grid.n_cells, 1)
    pia.n_total[1] = 4
    pia.indexer[1,1].n_local = 4
    pia.indexer[1,1].n_group1 = 4
    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = 4
    pia.contiguous[1] = true

    # will just move: new x_coord = 20.5
    particles[1][1] = Particle(1.0, [-1.25, -1.5, 4.0], [23.0, -8.0, 7.5])

    # will wrap around: new x_cord = 49.0 + 2 * 11.0 = 71.0 -> 71.0 - 50.0 = 21.0
    particles[1][2] = Particle(2.0, [11.0, -3.0, 1.0], [49.0, 6.0, -3.0])

    # will wrap around: 17.0 - 20.0 * 2.0 = -23.0 -> -23.0 + 50.0 = 27.0
    particles[1][3] = Particle(3.0, [-20.0, 0.0, 2.0], [17.0, 1.0, 3.0])

    # will wrap around twice: 1.55 - 49.0*2 = -96.45 -> -96.45 + 50.0*2 = 3.55
    particles[1][4] = Particle(4.0, [-49.0, -20.0, 13.0], [1.55, -1.0, 9.0])

    convect_particles_periodic!(grid, particles[1], pia, 1, 2.0)
    sort_particles!(gridsorter, grid, particles[1], pia, 1)

    @test particles[1].index == [4, 1, 2, 3]

    @test maximum(abs.(particles[1][1].x - [3.55, -1.0, 9.0])) < 3.65e-15 # * eps()
    @test particles[1][1].v == [-49.0, -20.0, 13.0]
    @test particles[1][1].w == 4.0

    @test maximum(abs.(particles[1][2].x - [20.5, -8.0, 7.5])) < 2 * eps()
    @test particles[1][2].v == [-1.25, -1.5, 4.0]
    @test particles[1][2].w == 1.0

    @test maximum(abs.(particles[1][3].x - [21.0, 6.0, -3.0])) < 2 * eps()
    @test particles[1][3].v == [11.0, -3.0, 1.0]
    @test particles[1][3].w == 2.0

    @test maximum(abs.(particles[1][4].x - [27.0, 1.0, 3.0])) < 2 * eps()
    @test particles[1][4].v == [-20.0, 0.0, 2.0]
    @test particles[1][4].w == 3.0

    phys_props = PhysProps(grid.n_cells, 1)
    compute_props!(particles, pia, species_data, phys_props)

    for i in 1:grid.n_cells
        if i == 8
            w = 4.0
        elseif i == 42
            w = 1.0
        elseif i == 43
            w = 2.0
        elseif i == 55
            w = 3.0
        else
            w = 0.0
        end
        @test abs(phys_props.n[i, 1] - w) < eps()
    end

    # non-contiguous indexing
    pia = ParticleIndexerArray(grid.n_cells, 1)
    pia.n_total[1] = 4
    pia.indexer[1,1].n_local = 4
    pia.indexer[1,1].n_group1 = 2
    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = 2
    pia.indexer[1,1].n_group2 = 2
    pia.indexer[1,1].start2 = 4
    pia.indexer[1,1].end2 = 5
    pia.contiguous[1] = false

    particles = [ParticleVector(5)]

    # will just move: new x_coord = 20.5
    particles[1][1] = Particle(1.0, [-1.25, -1.5, 4.0], [23.0, -8.0, 7.5])

    # will wrap around: new x_cord = 49.0 + 2 * 11.0 = 71.0 -> 71.0 - 50.0 = 21.0
    particles[1][2] = Particle(2.0, [11.0, -3.0, 1.0], [49.0, 6.0, -3.0])

    # just to check that this doesn't pop up somewhere
    particles[1][3] = Particle(300.0, [-20.0, 0.0, 2.0], [17.0, 1.0, 3.0])

    # will wrap around: 17.0 - 20.0 * 2.0 = -23.0 -> -23.0 + 50.0 = 27.0
    particles[1][4] = Particle(3.0, [-20.0, 0.0, 2.0], [17.0, 1.0, 3.0])

    # will wrap around twice: 1.55 - 49.0*2 = -96.45 -> -96.45 + 50.0*2 = 3.55
    particles[1][5] = Particle(4.0, [-49.0, -20.0, 13.0], [1.55, -1.0, 9.0])

    convect_particles_periodic!(grid, particles[1], pia, 1, 2.0)
    sort_particles!(gridsorter, grid, particles[1], pia, 1)

    @test particles[1].index[1:4] == [5, 1, 2, 4]

    @test maximum(abs.(particles[1][1].x - [3.55, -1.0, 9.0])) < 3.65e-15 # * eps()
    @test particles[1][1].v == [-49.0, -20.0, 13.0]
    @test particles[1][1].w == 4.0

    @test maximum(abs.(particles[1][2].x - [20.5, -8.0, 7.5])) < 2 * eps()
    @test particles[1][2].v == [-1.25, -1.5, 4.0]
    @test particles[1][2].w == 1.0

    @test maximum(abs.(particles[1][3].x - [21.0, 6.0, -3.0])) < 2 * eps()
    @test particles[1][3].v == [11.0, -3.0, 1.0]
    @test particles[1][3].w == 2.0

    @test maximum(abs.(particles[1][4].x - [27.0, 1.0, 3.0])) < 2 * eps()
    @test particles[1][4].v == [-20.0, 0.0, 2.0]
    @test particles[1][4].w == 3.0

    phys_props = PhysProps(grid.n_cells, 1)
    compute_props!(particles, pia, species_data, phys_props)

    for i in 1:grid.n_cells
        if i == 8
            w = 4.0
        elseif i == 42
            w = 1.0
        elseif i == 43
            w = 2.0
        elseif i == 55
            w = 3.0
        else
            w = 0.0
        end
        @test abs(phys_props.n[i, 1] - w) < eps()
    end

    # computing cells in convection routine
    particles = [ParticleVector(4)]

    pia = ParticleIndexerArray(grid.n_cells, 1)
    pia.n_total[1] = 4
    pia.indexer[1,1].n_local = 4
    pia.indexer[1,1].n_group1 = 4
    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = 4
    pia.contiguous[1] = true

    # will just move: new x_coord = 20.5
    particles[1][1] = Particle(1.0, [-1.25, -1.5, 4.0], [23.0, -8.0, 7.5])

    # will wrap around: new x_cord = 49.0 + 2 * 11.0 = 71.0 -> 71.0 - 50.0 = 21.0
    particles[1][2] = Particle(2.0, [11.0, -3.0, 1.0], [49.0, 6.0, -3.0])

    # will wrap around: 17.0 - 20.0 * 2.0 = -23.0 -> -23.0 + 50.0 = 27.0
    particles[1][3] = Particle(3.0, [-20.0, 0.0, 2.0], [17.0, 1.0, 3.0])

    # will wrap around twice: 1.55 - 49.0*2 = -96.45 -> -96.45 + 50.0*2 = 3.55
    particles[1][4] = Particle(4.0, [-49.0, -20.0, 13.0], [1.55, -1.0, 9.0])

    convect_particles_and_compute_cell_periodic!(grid, particles[1], pia, 1, 2.0)

    # cells computed correctly in routine
    @test particles[1].cell[1:4] == [42, 43, 55, 8]
    sort_particles!(gridsorter, particles[1], pia, 1)

    @test particles[1].index == [4, 1, 2, 3]

    @test maximum(abs.(particles[1][1].x - [3.55, -1.0, 9.0])) < 3.65e-15 # * eps()
    @test particles[1][1].v == [-49.0, -20.0, 13.0]
    @test particles[1][1].w == 4.0

    @test maximum(abs.(particles[1][2].x - [20.5, -8.0, 7.5])) < 2 * eps()
    @test particles[1][2].v == [-1.25, -1.5, 4.0]
    @test particles[1][2].w == 1.0

    @test maximum(abs.(particles[1][3].x - [21.0, 6.0, -3.0])) < 2 * eps()
    @test particles[1][3].v == [11.0, -3.0, 1.0]
    @test particles[1][3].w == 2.0

    @test maximum(abs.(particles[1][4].x - [27.0, 1.0, 3.0])) < 2 * eps()
    @test particles[1][4].v == [-20.0, 0.0, 2.0]
    @test particles[1][4].w == 3.0

    phys_props = PhysProps(grid.n_cells, 1)
    compute_props!(particles, pia, species_data, phys_props)

    for i in 1:grid.n_cells
        if i == 8
            w = 4.0
        elseif i == 42
            w = 1.0
        elseif i == 43
            w = 2.0
        elseif i == 55
            w = 3.0
        else
            w = 0.0
        end
        @test abs(phys_props.n[i, 1] - w) < eps()
    end

    # computing cells in convection routine, non-contiguous indexing
    particles = [ParticleVector(4)]
    pia = ParticleIndexerArray(grid.n_cells, 1)
    pia.n_total[1] = 4
    pia.indexer[1,1].n_local = 4
    pia.indexer[1,1].n_group1 = 2
    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = 2
    pia.indexer[1,1].n_group2 = 2
    pia.indexer[1,1].start2 = 4
    pia.indexer[1,1].end2 = 5
    pia.contiguous[1] = false

    particles = [ParticleVector(5)]

    # will just move: new x_coord = 20.5
    particles[1][1] = Particle(1.0, [-1.25, -1.5, 4.0], [23.0, -8.0, 7.5])

    # will wrap around: new x_cord = 49.0 + 2 * 11.0 = 71.0 -> 71.0 - 50.0 = 21.0
    particles[1][2] = Particle(2.0, [11.0, -3.0, 1.0], [49.0, 6.0, -3.0])

    # just to check that this doesn't pop up somewhere
    particles[1][3] = Particle(300.0, [-20.0, 0.0, 2.0], [17.0, 1.0, 3.0])

    # will wrap around: 17.0 - 20.0 * 2.0 = -23.0 -> -23.0 + 50.0 = 27.0
    particles[1][4] = Particle(3.0, [-20.0, 0.0, 2.0], [17.0, 1.0, 3.0])

    # will wrap around twice: 1.55 - 49.0*2 = -96.45 -> -96.45 + 50.0*2 = 3.55
    particles[1][5] = Particle(4.0, [-49.0, -20.0, 13.0], [1.55, -1.0, 9.0])

    convect_particles_and_compute_cell_periodic!(grid, particles[1], pia, 1, 2.0)

    # cells computed correctly in routine: 3rd cell is 0 as particle is skipped
    @test particles[1].cell[1:5] == [42, 43, 0, 55, 8]
    sort_particles!(gridsorter, particles[1], pia, 1)

    @test particles[1].index[1:4] == [5, 1, 2, 4]

    @test maximum(abs.(particles[1][1].x - [3.55, -1.0, 9.0])) < 3.65e-15 # * eps()
    @test particles[1][1].v == [-49.0, -20.0, 13.0]
    @test particles[1][1].w == 4.0

    @test maximum(abs.(particles[1][2].x - [20.5, -8.0, 7.5])) < 2 * eps()
    @test particles[1][2].v == [-1.25, -1.5, 4.0]
    @test particles[1][2].w == 1.0

    @test maximum(abs.(particles[1][3].x - [21.0, 6.0, -3.0])) < 2 * eps()
    @test particles[1][3].v == [11.0, -3.0, 1.0]
    @test particles[1][3].w == 2.0

    @test maximum(abs.(particles[1][4].x - [27.0, 1.0, 3.0])) < 2 * eps()
    @test particles[1][4].v == [-20.0, 0.0, 2.0]
    @test particles[1][4].w == 3.0

    phys_props = PhysProps(grid.n_cells, 1)
    compute_props!(particles, pia, species_data, phys_props)

    for i in 1:grid.n_cells
        if i == 8
            w = 4.0
        elseif i == 42
            w = 1.0
        elseif i == 43
            w = 2.0
        elseif i == 55
            w = 3.0
        else
            w = 0.0
        end
        @test abs(phys_props.n[i, 1] - w) < eps()
    end
end