function run_chunked_deposition(n_chunks; threaded)
    L = 5e-3
    nx = 64
    ppc = 200
    ndens = 1e15
    T_e = 11604.0
    T_i = 300.0

    seed = 1234
    rng_chunks = [StableRNG(seed + i) for i in 0:n_chunks-1]

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, ["e-", "He+"])

    grid = Grid1DUniform(L, nx)

    cell_chunks = chunks(Vector(1:nx); n=n_chunks)

    n_particles_chunks = [2 * ppc * length(cell_chunk) for cell_chunk in cell_chunks]
    particles_chunks = [[ParticleVector(n_particles), ParticleVector(n_particles)]
                        for n_particles in n_particles_chunks]
    pia_chunks = [ParticleIndexerArray(grid.n_cells, 2) for cell_chunk in cell_chunks]

    Fnum = grid.cells[1].V * ndens / ppc

    for (chunk_id, cell_chunk) in enumerate(cell_chunks)
        sample_particles_equal_weight!(rng_chunks[chunk_id], grid, particles_chunks[chunk_id][1],
                                       pia_chunks[chunk_id], 1, species_data, ppc, T_e, Fnum, cell_chunk)
        sample_particles_equal_weight!(rng_chunks[chunk_id], grid, particles_chunks[chunk_id][2],
                                       pia_chunks[chunk_id], 2, species_data, ppc, T_i, Fnum, cell_chunk)
    end

    field_props = ElectrostaticFieldProps(grid)
    field_props_chunks = [ElectrostaticFieldProps(grid) for cell_chunk in cell_chunks]
    poisson_solver = PoissonSolver1DUniform(grid, DirichletFieldBC1D(0.0), NeumannFieldBC1D(0.0))

    function deposit_chunk!(chunk_id)
        clear_charge_density!(field_props_chunks[chunk_id])

        for species in 1:2
            deposit_charge!(grid, particles_chunks[chunk_id][species], pia_chunks[chunk_id], species,
                            species_data, field_props_chunks[chunk_id], cell_chunks[chunk_id])
        end
    end

    if threaded
        @threads for chunk_id in 1:n_chunks
            deposit_chunk!(chunk_id)
        end
    else
        for chunk_id in 1:n_chunks
            deposit_chunk!(chunk_id)
        end
    end

    reduce_field_props!(field_props, field_props_chunks)
    normalize_charge_density!(poisson_solver, field_props)
    solve_poisson!(poisson_solver, field_props)

    return (charge_density=copy(field_props.charge_density),
            potential=copy(field_props.potential),
            electric_field=copy(field_props.electric_field))
end

@testset "chunked charge deposition: threaded results identical to serial" begin
    n_chunks = 4

    if Threads.nthreads() == 1
        @info "chunked charge deposition threading test running on 1 thread: no concurrency is exercised"
    end

    # each chunk deposits the cells it owns into its own ElectrostaticFieldProps instance and the
    # instances are summed in a fixed order afterwards, so neither the thread count nor the order in
    # which the chunks are processed may change a single bit of the answer
    serial = run_chunked_deposition(n_chunks; threaded=false)
    threaded = run_chunked_deposition(n_chunks; threaded=true)

    @test threaded.charge_density == serial.charge_density
    @test threaded.potential == serial.potential
    @test threaded.electric_field == serial.electric_field

    # the deposited charge is not identically zero, so the comparison above is not trivial
    @test maximum(abs.(serial.charge_density)) > 0.0
    @test maximum(abs.(serial.electric_field)) > 0.0
end
