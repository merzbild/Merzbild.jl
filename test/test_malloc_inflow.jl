
@testset "malloc tests: inflow velocity sampling" begin

    T = 300.0
    ndens = 1e20
    v_wall = 500.0

    seed = 1234
    rng = StableRNG(seed)

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")

    normal_left = SVector{3,Float64}(1.0, 0.0, 0.0)
    normal_right = SVector{3,Float64}(-1.0, 0.0, 0.0)

    bc_no_drift = InflowBC1D(1, species_data, T, SVector{3,Float64}(0.0, v_wall, 0.0), normal_left, ndens)
    bc_envelope4 = InflowBC1D(1, species_data, T, SVector{3,Float64}(v_wall, 0.0, 0.0), normal_left, ndens)
    bc_envelope4_right = InflowBC1D(1, species_data, T, SVector{3,Float64}(-v_wall, 0.0, 0.0), normal_right, ndens)

    @test bc_no_drift.sampler === InflowNoDrift()
    @test bc_envelope4.sampler === InflowEnvelope4()

    # sampling many particles in a loop, as an emission routine would
    function sample_many(rng, inflowbc, n)
        v = SVector{3,Float64}(0.0, 0.0, 0.0)
        for _ in 1:n
            v += sample_inflow_velocity(rng, inflowbc)
        end
        return v
    end

    # The BCs in the tuple have different types, so a BC taken from it is abstractly typed and the
    # returned velocity gets boxed - on Julia LTS that shows up as spurious allocations. Measuring
    # inside a function gives the BC a concrete type, which is also how it is used in a simulation.
    function check_sampling_allocations(rng, inflowbc, n)
        sample_many(rng, inflowbc, 100)

        # the acceptance-rejection loop must not allocate over many samples either
        bytes = @allocated sample_many(rng, inflowbc, n)
        @test bytes == 0
    end

    function check_flux_allocations(species_data, T, v, normal, ndens)
        compute_inflow_number_flux(1, species_data, T, v, normal, ndens)
        bytes = @allocated compute_inflow_number_flux(1, species_data, T, v, normal, ndens)
        @test bytes == 0
    end

    for bc in (bc_no_drift, bc_envelope4, bc_envelope4_right)
        check_sampling_allocations(rng, bc, 10000)
    end

    check_flux_allocations(species_data, T, SVector{3,Float64}(v_wall, 0.0, 0.0), normal_left, ndens)
end

@testset "malloc tests: emission across inflow boundaries" begin

    T = 300.0
    L = 5e-4
    ndens = 5e22
    nx = 50
    ppc = 100
    Δt = 2.59e-9

    seed = 1234
    rng = StableRNG(seed)

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")

    grid = Grid1DUniform(L, nx)

    # both reservoirs are in equilibrium with the domain, so the number of particles stays close to
    # the initial one, and the particle vector does not have to be resized
    bc_list = (InflowBC1D(1, species_data, T, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], ndens),
               InflowBC1D(1, species_data, T, [0.0, 0.0, 0.0], [-1.0, 0.0, 0.0], ndens))

    n_particles = 4 * ppc * nx
    particles = ParticleVector(n_particles)
    pia = ParticleIndexerArray(grid.n_cells, 1)
    gridsorter = GridSortInPlace(grid, n_particles)
    surf_props = SurfProps(pia, grid)

    Fnum = grid.cells[1].V * ndens / ppc
    sample_particles_equal_weight!(rng, grid, particles, pia, 1, species_data, ndens, T, Fnum)

    for t in 1:2
        convect_particles!(rng, grid, bc_list, particles, pia, 1, species_data, Δt)
        emit_particles!(rng, grid, bc_list, particles, pia, 1, species_data, Fnum, Δt)
        sort_particles!(gridsorter, grid, particles, pia, 1)

        convect_particles!(rng, grid, bc_list, particles, pia, 1, species_data, surf_props, Δt)
        emit_particles!(rng, grid, bc_list, particles, pia, 1, species_data, surf_props, Fnum, Δt)
        sort_particles!(gridsorter, grid, particles, pia, 1)

        convect_particles_and_compute_cell!(rng, grid, bc_list, particles, pia, 1, species_data, Δt)
        emit_particles!(rng, grid, bc_list, particles, pia, 1, species_data, Fnum, Δt, 1:nx)
        sort_particles!(gridsorter, particles, pia, 1)
    end

    n_particles_allocated = length(particles)

    for t in 1:2
        bytes = @allocated convect_particles!(rng, grid, bc_list, particles, pia, 1, species_data, Δt)
        @test bytes == 0
        bytes = @allocated emit_particles!(rng, grid, bc_list, particles, pia, 1, species_data, Fnum, Δt)
        @test bytes == 0
        bytes = @allocated sort_particles!(gridsorter, grid, particles, pia, 1)
        @test bytes == 0

        bytes = @allocated convect_particles!(rng, grid, bc_list, particles, pia, 1, species_data, surf_props, Δt)
        @test bytes == 0
        bytes = @allocated emit_particles!(rng, grid, bc_list, particles, pia, 1, species_data, surf_props, Fnum, Δt)
        @test bytes == 0
        bytes = @allocated sort_particles!(gridsorter, grid, particles, pia, 1)
        @test bytes == 0

        bytes = @allocated convect_particles_and_compute_cell!(rng, grid, bc_list, particles, pia, 1, species_data, Δt)
        @test bytes == 0
        bytes = @allocated emit_particles!(rng, grid, bc_list, particles, pia, 1, species_data, Fnum, Δt, 1:nx)
        @test bytes == 0
        bytes = @allocated sort_particles!(gridsorter, particles, pia, 1)
        @test bytes == 0
    end

    @test length(particles) == n_particles_allocated
    @test check_pia_is_correct(pia, 1) == (true, 0)
    @test check_unique_index(particles, pia, 1) == (true, 0)
    @test check_unique_buffer(particles) == (true, 0)
end
