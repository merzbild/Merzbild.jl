@testset "malloc tests: electrostatic PIC" begin

    # The Poisson solver is parametrized on the types of the boundary conditions, so a solver taken
    # from a list of solvers with different BCs is abstractly typed; the routines called with such
    # an instance then access its fields generically, which boxes them - on Julia LTS that shows up
    # as spurious allocations. Measuring inside a function gives the solver a concrete type, which
    # is also how it is used in a simulation.
    function check_allocations(poisson_solver, grid, particles, pia, species_data, field_props, Δt)
        for t in 1:2
            clear_charge_density!(field_props)
            deposit_charge!(poisson_solver, grid, particles, pia, species_data, field_props)
            normalize_charge_density!(poisson_solver, field_props)
            solve_poisson!(poisson_solver, field_props)
            accelerate_electric_field_x!(grid, particles[1], pia, 1, species_data, field_props, Δt)
            accelerate_electric_field_x!(grid, particles[2], pia, 1, 2, species_data, field_props, Δt)
            deposit_charge!(grid, particles[1], pia, 1, species_data, field_props)
        end

        for t in 1:2
            bytes = @allocated clear_charge_density!(field_props)
            @test bytes == 0

            bytes = @allocated deposit_charge!(poisson_solver, grid, particles, pia, species_data, field_props)
            @test bytes == 0

            bytes = @allocated normalize_charge_density!(poisson_solver, field_props)
            @test bytes == 0

            bytes = @allocated solve_poisson!(poisson_solver, field_props)
            @test bytes == 0

            bytes = @allocated accelerate_electric_field_x!(grid, particles[1], pia, 1, species_data,
                                                            field_props, Δt)
            @test bytes == 0

            bytes = @allocated accelerate_electric_field_x!(grid, particles[2], pia, 1, 2, species_data,
                                                            field_props, Δt)
            @test bytes == 0

            bytes = @allocated deposit_charge!(grid, particles[1], pia, 1, species_data, field_props)
            @test bytes == 0
        end
    end

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, ["e-", "He+"])

    L = 0.01
    nx = 32
    ppc = 20
    ndens = 1e15
    Δt = 1e-11

    seed = 1234
    rng = StableRNG(seed)

    grid = Grid1DUniform(L, nx)

    n_particles = ppc * nx
    particles = [ParticleVector(n_particles), ParticleVector(n_particles)]
    pia = ParticleIndexerArray(grid.n_cells, 2)
    gridsorter = GridSortInPlace(grid, n_particles)

    Fnum = grid.cells[1].V * ndens / ppc
    sample_particles_equal_weight!(rng, grid, particles[1], pia, 1, species_data, ppc, 11604.0, Fnum)
    sample_particles_equal_weight!(rng, grid, particles[2], pia, 2, species_data, ppc, 300.0, Fnum)

    field_props = ElectrostaticFieldProps(grid)

    for (bc_left, bc_right) in [(PeriodicFieldBC1D(), PeriodicFieldBC1D()),
                                (DirichletFieldBC1D(0.0), DirichletFieldBC1D(50.0)),
                                (NeumannFieldBC1D(0.0), DirichletFieldBC1D(0.0))]
        poisson_solver = PoissonSolver1DUniform(grid, bc_left, bc_right)
        check_allocations(poisson_solver, grid, particles, pia, species_data, field_props, Δt)
    end
end
