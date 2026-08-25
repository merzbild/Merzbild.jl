function run_chunked_plasma_oscillation(n_chunks, n_timesteps; threaded, factorized_exchange)
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, ["e-", "He+"])

    L = 0.01
    nx = 32
    ppc = 200
    ndens = 1e14
    preallocation_margin_multiplier = 1.5

    ω_p = plasma_frequency(1, species_data, ndens)
    Δt = 0.1 / ω_p

    grid = Grid1DUniform(L, nx)

    cell_chunks = chunks(Vector(1:nx); n=n_chunks)
    chunk_of_cell = zeros(Int64, nx)
    for (chunk_id, cell_chunk) in enumerate(cell_chunks)
        for cell in cell_chunk
            chunk_of_cell[cell] = chunk_id
        end
    end

    n_particles = ppc * nx
    n_particles_chunks = [floor(Int64, ppc * length(cell_chunk) * preallocation_margin_multiplier)
                          for cell_chunk in cell_chunks]
    particles_chunks = [[ParticleVector(np), ParticleVector(np)] for np in n_particles_chunks]
    pia_chunks = [ParticleIndexerArray(grid.n_cells, 2) for cell_chunk in cell_chunks]
    gridsorter_chunks = [GridSortInPlace(grid, np) for np in n_particles_chunks]

    chunk_exchanger = ChunkExchanger(cell_chunks, nx)
    exchange_list = generate_1_factorization(n_chunks)

    w_particle = ndens * L / n_particles
    amplitude = 0.01 * L
    k_wave = 2π / L

    np_chunks = [zeros(Int64, 2) for cell_chunk in cell_chunks]
    v_zero = SVector{3, Float64}(0.0, 0.0, 0.0)

    for i in 1:n_particles
        x0 = L * (i - 0.5) / n_particles
        x_perturbed = x0 + amplitude * sin(k_wave * x0)

        for (species, x) in ((1, x_perturbed), (2, x0))
            x_vec = SVector{3, Float64}(x, 0.0, 0.0)
            chunk_id = chunk_of_cell[Merzbild.get_cell(grid, x_vec)]

            np_chunks[chunk_id][species] += 1
            Merzbild.add_particle!(particles_chunks[chunk_id][species], np_chunks[chunk_id][species],
                                   w_particle, v_zero, x_vec)
        end
    end

    for chunk_id in 1:n_chunks
        pia = pia_chunks[chunk_id]

        for species in 1:2
            np = np_chunks[chunk_id][species]

            pia.n_total[species] = np
            pia.index_last[species] = np
            pia.indexer[1,species].n_local = np
            pia.indexer[1,species].n_group1 = np
            pia.indexer[1,species].start1 = 1
            pia.indexer[1,species].end1 = np
            pia.contiguous[species] = true

            sort_particles!(gridsorter_chunks[chunk_id], grid, particles_chunks[chunk_id][species], pia, species)
        end
    end

    field_props = ElectrostaticFieldProps(grid)
    field_props_chunks = [ElectrostaticFieldProps(grid) for cell_chunk in cell_chunks]
    poisson_solver = PoissonSolver1DUniform(grid, PeriodicFieldBC1D(), PeriodicFieldBC1D())

    function run_chunks!(f)
        if threaded
            @threads for chunk_id in 1:n_chunks
                f(chunk_id)
            end
        else
            for chunk_id in 1:n_chunks
                f(chunk_id)
            end
        end
    end

    # each chunk deposits the cells it owns into its own ElectrostaticFieldProps instance, as the
    # nodes shared by cells of different chunks would otherwise be written to concurrently
    function deposit_chunk!(chunk_id)
        clear_charge_density!(field_props_chunks[chunk_id])

        for species in 1:2
            deposit_charge!(grid, particles_chunks[chunk_id][species], pia_chunks[chunk_id], species,
                            species_data, field_props_chunks[chunk_id], cell_chunks[chunk_id])
        end
    end

    function solve_field!()
        run_chunks!(deposit_chunk!)

        reduce_field_props!(field_props, field_props_chunks)
        normalize_charge_density!(poisson_solver, field_props)
        solve_poisson!(poisson_solver, field_props)
    end

    # the gather only reads the field, so all the chunks push into the shared instance
    function push_chunk!(chunk_id, Δt_push)
        for cell in cell_chunks[chunk_id]
            accelerate_electric_field_x!(grid, particles_chunks[chunk_id][1], pia_chunks[chunk_id], cell, 1,
                                         species_data, field_props, Δt_push)
        end
    end

    # the ions are a fixed background and are never pushed, convected or exchanged
    function advance_chunk!(chunk_id)
        push_chunk!(chunk_id, Δt)

        convect_particles_periodic!(grid, particles_chunks[chunk_id][1], pia_chunks[chunk_id], 1, Δt)

        sort_particles!(gridsorter_chunks[chunk_id], grid, particles_chunks[chunk_id][1], pia_chunks[chunk_id], 1)
        update_occupancy_bounds!(chunk_exchanger, gridsorter_chunks[chunk_id], pia_chunks[chunk_id], chunk_id, 1)
    end

    function restore_indexing!(chunk_id)
        sort_particles_after_exchange!(chunk_exchanger, gridsorter_chunks[chunk_id],
                                       particles_chunks[chunk_id][1], pia_chunks[chunk_id],
                                       cell_chunks[chunk_id], 1)
    end

    function exchange!()
        if factorized_exchange
            for exchange_partial in exchange_list
                if threaded
                    @threads for exchange_pair in exchange_partial
                        exchange_particles!(chunk_exchanger, particles_chunks, pia_chunks, cell_chunks,
                                            1, exchange_pair[1], exchange_pair[2])
                    end
                else
                    for exchange_pair in exchange_partial
                        exchange_particles!(chunk_exchanger, particles_chunks, pia_chunks, cell_chunks,
                                            1, exchange_pair[1], exchange_pair[2])
                    end
                end
            end
        else
            exchange_particles!(chunk_exchanger, particles_chunks, pia_chunks, cell_chunks, 1)
        end
    end

    # leapfrog initialization: a single backward half-kick of the electrons
    solve_field!()
    net_charge_density = field_props.net_charge_density
    run_chunks!(chunk_id -> push_chunk!(chunk_id, -0.5 * Δt))

    monitor_node = 1 + nx ÷ 4
    E_monitor = zeros(n_timesteps)
    energy_total = zeros(n_timesteps)

    m_e = species_data[1].mass

    for t in 1:n_timesteps
        solve_field!()

        run_chunks!(advance_chunk!)
        exchange!()
        run_chunks!(restore_indexing!)

        E_monitor[t] = field_props.electric_field[monitor_node]

        energy_field = 0.0
        for j in 1:nx
            energy_field += 0.5 * Merzbild.eps_0 * field_props.electric_field[j]^2 * grid.Δx
        end

        energy_kinetic = 0.0
        for chunk_id in 1:n_chunks
            for i in 1:pia_chunks[chunk_id].n_total[1]
                p = particles_chunks[chunk_id][1][i]
                energy_kinetic += 0.5 * m_e * p.w * p.v[1]^2
            end
        end

        energy_total[t] = energy_field + energy_kinetic
    end

    return (E_monitor=E_monitor, energy_total=energy_total,
            charge_density=copy(field_props.charge_density),
            electric_field=copy(field_props.electric_field),
            net_charge_density=net_charge_density,
            n_total=[pia_chunks[chunk_id].n_total[1] for chunk_id in 1:n_chunks],
            pia_correct=[check_pia_is_correct(pia_chunks[chunk_id], 1) for chunk_id in 1:n_chunks],
            index_correct=[check_unique_index(particles_chunks[chunk_id][1], pia_chunks[chunk_id], 1) for chunk_id in 1:n_chunks],
            buffer_correct=[check_unique_buffer(particles_chunks[chunk_id][1]) for chunk_id in 1:n_chunks],
            n_particles=n_particles, ndens=ndens, ω_p=ω_p, Δt=Δt)
end

@testset "chunked cold plasma oscillation: threaded results identical to serial" begin
    n_chunks = 4
    n_timesteps = 1200

    if Threads.nthreads() == 1
        @info "chunked plasma oscillation threading test running on 1 thread: no concurrency is exercised"
    end

    # a cold electron plasma in a periodic box on top of a fixed neutralizing ion background is
    # perturbed sinusoidally in space and oscillates at the plasma frequency; the electrons are
    # deposited, pushed and convected per cell chunk, and cross chunk boundaries as they oscillate
    serial = run_chunked_plasma_oscillation(n_chunks, n_timesteps; threaded=false, factorized_exchange=false)

    # the simulation contains no randomness, each chunk deposits into its own field props instance
    # and the instances are summed in a fixed order, and sort_particles_after_exchange! restores a
    # canonical cell-sorted ordering every timestep, so neither the thread count nor the order in
    # which the chunks are processed may change a single bit of the answer
    for (threaded, factorized_exchange) in ((true, false), (false, true), (true, true))
        result = run_chunked_plasma_oscillation(n_chunks, n_timesteps;
                                                threaded=threaded, factorized_exchange=factorized_exchange)

        @test result.E_monitor == serial.E_monitor
        @test result.energy_total == serial.energy_total
        @test result.charge_density == serial.charge_density
        @test result.electric_field == serial.electric_field
        @test result.n_total == serial.n_total

        for chunk_id in 1:n_chunks
            @test result.pia_correct[chunk_id] == (1, 0)
            @test result.index_correct[chunk_id] == (1, 0)
            @test result.buffer_correct[chunk_id] == (1, 0)
        end
    end

    # no electron is lost or duplicated by the exchange
    @test sum(serial.n_total) == serial.n_particles

    # the box is neutral overall, so no net charge density should be subtracted
    @test abs(serial.net_charge_density) < 1e-12 * Merzbild.q_e * serial.ndens

    # measure the oscillation period from the upward zero crossings of the monitored field,
    # linearly interpolating to get the crossing time
    crossings = Float64[]
    for t in 2:n_timesteps
        if serial.E_monitor[t-1] < 0.0 && serial.E_monitor[t] >= 0.0
            push!(crossings, (t - 1 + serial.E_monitor[t-1] / (serial.E_monitor[t-1] - serial.E_monitor[t])) * serial.Δt)
        end
    end

    @test length(crossings) > 15

    period = (crossings[end] - crossings[1]) / (length(crossings) - 1)
    ω_measured = 2π / period

    @test abs(ω_measured / serial.ω_p - 1.0) < 0.02

    # the total energy stays bounded over the whole simulation; a cold plasma with Δx ≫ λ_D
    # heats slowly due to the finite-grid instability, so a few percent of growth is expected
    n_steps_per_period = round(Int64, 2π / (serial.ω_p * serial.Δt))
    energy_scale = maximum(serial.energy_total[1:n_steps_per_period])
    @test abs(maximum(serial.energy_total) / energy_scale - 1.0) < 0.1
    @test minimum(serial.energy_total) > 0.0
end
