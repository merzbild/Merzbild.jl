@testset "charge deposition" begin
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, ["e-", "He+", "Ar"])

    q_e = Merzbild.q_e

    L = 1.0
    nx = 4
    grid = Grid1DUniform(L, nx)
    Δx = grid.Δx

    n_particles = 8
    particles = [ParticleVector(n_particles), ParticleVector(n_particles), ParticleVector(n_particles)]
    pia = ParticleIndexerArray(grid.n_cells, 3)
    gridsorter = GridSortInPlace(grid, n_particles)

    field_props = ElectrostaticFieldProps(grid)

    function set_particles!(species, xs, ws)
        for cell in 1:grid.n_cells
            pia.indexer[cell,species] = ParticleIndexer()
        end

        pia.n_total[species] = length(xs)
        pia.index_last[species] = length(xs)
        pia.indexer[1,species].n_local = length(xs)
        pia.indexer[1,species].n_group1 = length(xs)
        pia.indexer[1,species].start1 = 1
        pia.indexer[1,species].end1 = length(xs)
        pia.contiguous[species] = true

        for (i, (xp, wp)) in enumerate(zip(xs, ws))
            particles[species][i] = Particle(wp, [0.0, 0.0, 0.0], [xp, 0.0, 0.0])
        end

        sort_particles!(gridsorter, grid, particles[species], pia, species)
    end

    # a single electron sitting exactly on node 2 deposits all its charge there
    w = 1e6
    set_particles!(1, [Δx], [w])
    clear_charge_density!(field_props)
    deposit_charge!(grid, particles[1], pia, 1, species_data, field_props)

    @test field_props.charge_density[2] == -q_e * w
    @test field_props.charge_density[1] == 0.0
    @test field_props.charge_density[3] == 0.0
    @test sum(field_props.charge_density) ≈ -q_e * w

    # a single ion at the centre of cell 1 splits its charge equally between nodes 1 and 2
    set_particles!(2, [0.5 * Δx], [w])
    clear_charge_density!(field_props)
    deposit_charge!(grid, particles[2], pia, 2, species_data, field_props)

    @test field_props.charge_density[1] ≈ 0.5 * q_e * w
    @test field_props.charge_density[2] ≈ 0.5 * q_e * w
    @test maximum(abs.(field_props.charge_density[3:end])) == 0.0

    # a neutral species contributes nothing
    set_particles!(3, [0.3 * Δx, 2.7 * Δx], [w, w])
    clear_charge_density!(field_props)
    deposit_charge!(grid, particles[3], pia, 3, species_data, field_props)
    @test maximum(abs.(field_props.charge_density)) == 0.0

    # deposition is additive across species
    x_e = [0.2 * Δx, 1.35 * Δx, 3.9 * Δx]
    w_e = [1e6, 2e6, 3e6]
    x_i = [0.75 * Δx, 2.5 * Δx]
    w_i = [4e6, 5e6]

    set_particles!(1, x_e, w_e)
    set_particles!(2, x_i, w_i)

    clear_charge_density!(field_props)
    deposit_charge!(grid, particles[1], pia, 1, species_data, field_props)
    ρ_electrons = copy(field_props.charge_density)

    deposit_charge!(grid, particles[2], pia, 2, species_data, field_props)
    ρ_both = copy(field_props.charge_density)

    clear_charge_density!(field_props)
    deposit_charge!(grid, particles[2], pia, 2, species_data, field_props)
    @test maximum(abs.(ρ_both .- ρ_electrons .- field_props.charge_density)) < 1e-14 * maximum(abs.(ρ_both))

    # total deposited charge is conserved exactly
    q_total = -q_e * sum(w_e) + q_e * sum(w_i)
    @test abs(sum(ρ_both) - q_total) < 1e-15 * abs(q_total)

    # charge conservation after the non-periodic normalization: the boundary nodes have a Δx/2 dual cell
    poisson_solver = PoissonSolver1DUniform(grid, DirichletFieldBC1D(0.0), DirichletFieldBC1D(0.0))

    deposit_charge!(poisson_solver, grid, particles, pia, species_data, field_props)

    q_integrated = 0.5 * field_props.charge_density[1] + 0.5 * field_props.charge_density[end]
    for j in 2:field_props.n_nodes-1
        q_integrated += field_props.charge_density[j]
    end
    q_integrated *= Δx
    @test abs(q_integrated - q_total) < 1e-14 * abs(q_total)

    # charge conservation after the periodic normalization: every node has a full Δx dual cell
    poisson_solver_periodic = PoissonSolver1DUniform(grid, PeriodicFieldBC1D(), PeriodicFieldBC1D())

    deposit_charge!(poisson_solver_periodic, grid, particles, pia, species_data, field_props)

    q_integrated = 0.0
    for j in 1:grid.n_cells
        q_integrated += field_props.charge_density[j]
    end
    q_integrated *= Δx
    @test abs(q_integrated - q_total) < 1e-14 * abs(q_total)
    @test field_props.charge_density[end] == field_props.charge_density[1]

    # the multi-species deposition clears the previously deposited charge instead of accumulating it
    ρ_deposited = copy(field_props.charge_density)
    deposit_charge!(poisson_solver_periodic, grid, particles, pia, species_data, field_props)
    @test field_props.charge_density == ρ_deposited

    # chunked deposition into per-chunk holders, reduced, reproduces the serial deposition
    for n_chunks in [2, 3]
        cell_chunks = chunks(Vector(1:grid.n_cells); n=n_chunks)
        field_props_chunks = [ElectrostaticFieldProps(grid) for cell_chunk in cell_chunks]

        for (chunk_id, cell_chunk) in enumerate(cell_chunks)
            clear_charge_density!(field_props_chunks[chunk_id])

            for species in 1:pia.n_species
                deposit_charge!(grid, particles[species], pia, species, species_data,
                                field_props_chunks[chunk_id], cell_chunk)
            end
        end

        # each chunk only deposits the cells it owns, so the per-chunk charges add up to the total
        q_chunks = 0.0
        for chunk_id in 1:n_chunks
            q_chunks += sum(field_props_chunks[chunk_id].charge_density)
        end
        @test abs(q_chunks - q_total) < 1e-14 * abs(q_total)

        field_props_reduced = ElectrostaticFieldProps(grid)
        fill!(field_props_reduced.potential, 7.0)
        fill!(field_props_reduced.electric_field, -9.0)
        fill!(field_props_reduced.charge_density, 1.0)

        reduce_field_props!(field_props_reduced, field_props_chunks)

        # the reduction clears the target's charge density and leaves the potential and field alone
        @test minimum(field_props_reduced.potential) == 7.0
        @test minimum(field_props_reduced.electric_field) == -9.0
        @test field_props_reduced.net_charge_density == 0.0

        clear_charge_density!(field_props)
        for species in 1:pia.n_species
            deposit_charge!(grid, particles[species], pia, species, species_data, field_props)
        end

        @test maximum(abs.(field_props_reduced.charge_density .- field_props.charge_density)) <
              1e-14 * maximum(abs.(field_props.charge_density))

        # the normalization runs once, on the reduced result
        normalize_charge_density!(poisson_solver_periodic, field_props_reduced)
        normalize_charge_density!(poisson_solver_periodic, field_props)

        @test maximum(abs.(field_props_reduced.charge_density .- field_props.charge_density)) <
              1e-14 * maximum(abs.(field_props.charge_density))
    end

    # a quasineutral two-species plasma has a vanishing charge density and field
    L_qn = 0.01
    nx_qn = 20
    ppc = 500
    grid_qn = Grid1DUniform(L_qn, nx_qn)

    seed = 1234
    rng = StableRNG(seed)

    ndens = 1e16
    n_particles_qn = ppc * nx_qn
    particles_qn = [ParticleVector(n_particles_qn), ParticleVector(n_particles_qn)]
    pia_qn = ParticleIndexerArray(grid_qn.n_cells, 2)
    Fnum = grid_qn.cells[1].V * ndens / ppc

    sample_particles_equal_weight!(rng, grid_qn, particles_qn[1], pia_qn, 1, species_data, ppc, 5000.0, Fnum)
    sample_particles_equal_weight!(rng, grid_qn, particles_qn[2], pia_qn, 2, species_data, ppc, 300.0, Fnum)

    # the ions are placed at the same positions as the electrons, so the plasma is exactly quasineutral
    for i in 1:pia_qn.n_total[1]
        particles_qn[2][i].x = particles_qn[1][i].x
    end

    field_props_qn = ElectrostaticFieldProps(grid_qn)
    poisson_solver_qn = PoissonSolver1DUniform(grid_qn, PeriodicFieldBC1D(), PeriodicFieldBC1D())

    deposit_charge!(poisson_solver_qn, grid_qn, particles_qn, pia_qn, species_data, field_props_qn)
    solve_poisson!(poisson_solver_qn, field_props_qn)

    ρ_scale = q_e * ndens
    @test maximum(abs.(field_props_qn.charge_density)) < 1e-12 * ρ_scale
    @test maximum(abs.(field_props_qn.electric_field)) < 1e-12 * ρ_scale * L_qn / Merzbild.eps_0
    @test abs(field_props_qn.net_charge_density) < 1e-12 * ρ_scale
end
