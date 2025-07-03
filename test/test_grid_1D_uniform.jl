@testset "1D uniform grid sampling and computes" begin
    Δlarge = 24.0

    seed = 1234
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    grid = Grid1DUniform(4.0, 8)

    @test grid.L == 4
    @test grid.n_cells == 8
    @test grid.Δx == 0.5

    xlo = 0.0
    xhi = 0.0

    for i in 1:8
        @test grid.cells[i].xlo == xlo
        @test grid.cells[i].xhi == xlo + 0.5
        @test grid.cells[i].V == 0.5

        xlo += 0.5
    end

    n_per_cell = 1e10
    ppc = 1000
    Fnum::Float64 = n_per_cell / ppc
    T = 500.0

    particles = [ParticleVector(ppc * grid.n_cells)]

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, "Ar")

    pia = ParticleIndexerArray(grid.n_cells, 1)

    sample_particles_equal_weight!(rng, grid, particles[1], pia, 1,
                                   species_data, ppc, T, Fnum)

    @test pia.n_total[1] == ppc * grid.n_cells
    for i in 1:grid.n_cells
        @test pia.indexer[i, 1].n_local == ppc

        @test pia.indexer[i, 1].start1 == 1 + ppc * (i-1)
        @test pia.indexer[i, 1].end1 == ppc * i
        @test pia.indexer[i, 1].n_group1 == ppc

        @test pia.indexer[i, 1].start2 == 0
        @test pia.indexer[i, 1].end2 == 0
        @test pia.indexer[i, 1].n_group2 == 0

        @test particles[1].cell[(i-1)*ppc + 1] == i
    end

    phys_props::PhysProps = PhysProps(grid.n_cells, 1, [], Tref=1)
    compute_props!(particles, pia, species_data, phys_props)

    @test phys_props.n_species == 1
    @test phys_props.n_cells == grid.n_cells

    for i in 1:grid.n_cells
        @test abs(phys_props.n[i, 1] - n_per_cell) < 2*eps()
        @test phys_props.np[i, 1] == ppc
        @test abs(phys_props.T[i, 1] - T) / T < 0.075

        @test abs((phys_props.v[1,i,1])) < Δlarge
        @test abs((phys_props.v[2,i,1])) < Δlarge
        @test abs((phys_props.v[3,i,1])) < Δlarge
    end

    @test sum(phys_props.np) == ppc * grid.n_cells

    @test Merzbild.get_cell(grid, [0.001, 0.0, 0.0]) == 1
    @test Merzbild.get_cell(grid, [0.4, 0.0, 0.0]) == 1
    @test Merzbild.get_cell(grid, [0.501, 0.0, 0.0]) == 2
    @test Merzbild.get_cell(grid, [3.999, 0.0, 0.0]) == 8

    # test filling domain with particles
    # based on given ndens
    # cell volume is 0.5, so cell will have 500 particles
    particles2 = [ParticleVector(ppc * grid.n_cells)]

    pia.n_total[1] = 0

    ndens = 1e23
    Fnum = 1e20
    sample_particles_equal_weight!(rng, grid, particles2[1], pia, 1,
                                   species_data, ndens, T, Fnum)

    compute_props!(particles2, pia, species_data, phys_props)

    # round-off errors, therefore we round up
    @test sum(phys_props.np) == 0.5 * round(Int64, ndens / Fnum) * grid.n_cells

    # we have 500 ppc instead of 1000, so we relax the constraints a bit
    for i in 1:grid.n_cells
        @test abs(phys_props.T[i, 1] - T) / T < 0.1

        @test abs((phys_props.v[1,i,1])) < Δlarge * 1.5
        @test abs((phys_props.v[2,i,1])) < Δlarge * 1.5
        @test abs((phys_props.v[3,i,1])) < Δlarge * 1.5
    end

    # test PIA construction from grid and species data
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_data(particles_data_path, ["Ar", "He"])

    pia = ParticleIndexerArray(grid, species_data)
    @test length(pia.n_total) == 2
    @test size(pia.indexer) == (8, 2)
end 