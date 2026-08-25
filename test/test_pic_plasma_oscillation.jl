@testset "cold plasma oscillation" begin
    # A cold electron plasma in a periodic box on top of a fixed neutralizing ion background is
    # perturbed sinusoidally in space; the electrons then oscillate at the plasma frequency.
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, ["e-", "He+"])

    L = 0.01
    nx = 32
    ppc = 200
    ndens = 1e14

    ω_p = plasma_frequency(1, species_data, ndens)
    Δt = 0.1 / ω_p
    n_timesteps = 1200

    λ_D = debye_length(ndens, 11604.0)
    @test λ_D > 0.0

    grid = Grid1DUniform(L, nx)

    n_particles = ppc * nx
    particles = [ParticleVector(n_particles), ParticleVector(n_particles)]
    pia = ParticleIndexerArray(grid.n_cells, 2)
    gridsorter = GridSortInPlace(grid, n_particles)

    w_particle = ndens * L / n_particles
    amplitude = 0.01 * L
    k_wave = 2π / L

    for i in 1:n_particles
        x0 = L * (i - 0.5) / n_particles
        x_perturbed = x0 + amplitude * sin(k_wave * x0)

        particles[1][i] = Particle(w_particle, [0.0, 0.0, 0.0], [x_perturbed, 0.0, 0.0])
        particles[2][i] = Particle(w_particle, [0.0, 0.0, 0.0], [x0, 0.0, 0.0])
    end

    for species in 1:2
        pia.n_total[species] = n_particles
        pia.index_last[species] = n_particles
        pia.indexer[1,species].n_local = n_particles
        pia.indexer[1,species].n_group1 = n_particles
        pia.indexer[1,species].start1 = 1
        pia.indexer[1,species].end1 = n_particles
        pia.contiguous[species] = true

        sort_particles!(gridsorter, grid, particles[species], pia, species)
    end

    field_props = ElectrostaticFieldProps(grid)
    poisson_solver = PoissonSolver1DUniform(grid, PeriodicFieldBC1D(), PeriodicFieldBC1D())

    # the box is neutral overall, so no net charge density should be subtracted
    deposit_charge!(poisson_solver, grid, particles, pia, species_data, field_props)
    solve_poisson!(poisson_solver, field_props)
    @test abs(field_props.net_charge_density) < 1e-12 * Merzbild.q_e * ndens

    # leapfrog initialization: a single backward half-kick of the electrons
    accelerate_electric_field_x!(grid, particles[1], pia, 1, species_data, field_props, -0.5 * Δt)

    monitor_node = 1 + nx ÷ 4
    E_monitor = zeros(n_timesteps)
    energy_total = zeros(n_timesteps)

    m_e = species_data[1].mass

    for t in 1:n_timesteps
        deposit_charge!(poisson_solver, grid, particles, pia, species_data, field_props)
        solve_poisson!(poisson_solver, field_props)

        # the ions are a fixed background and are never pushed
        accelerate_electric_field_x!(grid, particles[1], pia, 1, species_data, field_props, Δt)

        convect_particles_periodic!(grid, particles[1], pia, 1, Δt)
        sort_particles!(gridsorter, grid, particles[1], pia, 1)

        E_monitor[t] = field_props.electric_field[monitor_node]

        energy_field = 0.0
        for j in 1:nx
            energy_field += 0.5 * Merzbild.eps_0 * field_props.electric_field[j]^2 * grid.Δx
        end

        energy_kinetic = 0.0
        for i in 1:pia.n_total[1]
            energy_kinetic += 0.5 * m_e * particles[1][i].w * particles[1][i].v[1]^2
        end

        energy_total[t] = energy_field + energy_kinetic
    end

    # measure the oscillation period from the upward zero crossings of the monitored field
    # linear interpolate to get crossing time
    crossings = Float64[]
    for t in 2:n_timesteps
        if E_monitor[t-1] < 0.0 && E_monitor[t] >= 0.0
            push!(crossings, (t - 1 + E_monitor[t-1] / (E_monitor[t-1] - E_monitor[t])) * Δt)
        end
    end

    @test length(crossings) > 15

    period = (crossings[end] - crossings[1]) / (length(crossings) - 1)
    ω_measured = 2π / period

    @test abs(ω_measured / ω_p - 1.0) < 0.02

    # the total energy stays bounded over the whole simulation; a cold plasma with Δx ≫ λ_D
    # heats slowly due to the finite-grid instability, so a few percent of growth is expected
    n_steps_per_period = round(Int64, 2π / (ω_p * Δt))
    energy_scale = maximum(energy_total[1:n_steps_per_period])
    @test abs(maximum(energy_total) / energy_scale - 1.0) < 0.1
    @test minimum(energy_total) > 0.0
end
