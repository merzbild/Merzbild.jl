@testset "emission across 1D inflow boundaries" begin
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")
    mass = species_data[1].mass

    T = 300.0
    ndens = 1e20
    v_thermal = Merzbild.compute_thermal_velocity(mass, T)

    normal_left = SVector{3,Float64}(1.0, 0.0, 0.0)
    normal_right = SVector{3,Float64}(-1.0, 0.0, 0.0)
    v_zero = SVector{3,Float64}(0.0, 0.0, 0.0)

    bc_left = InflowBC1D(1, species_data, T, v_zero, normal_left, ndens)
    bc_right = InflowBC1D(1, species_data, T, v_zero, normal_right, ndens)

    # the particles emitted over a timestep stay far from the opposite wall
    L = 1.0
    nx = 50
    grid = Grid1DUniform(L, nx)
    Δt = 1e-3 * L / v_thermal

    # expected number of particles emitted across the left wall into an initially empty collisionless domain
    # with no particles coming back from the right wall, after a time t:
    # integral of (n/2) erfc((x/t - u) / v_thermal) over [0, L]
    G(s) = s * erfc(s) - exp(-s^2) / sqrt(π)
    function n_particles_free_flight(bc, L, t, Fnum)
        s0 = -bc.a
        s1 = L / (t * bc.v_thermal) - bc.a
        return bc.ndens / 2 * t * bc.v_thermal * (G(s1) - G(s0)) / Fnum
    end

    @testset "number of particles and residual flight time" begin
        λ = 50.0
        Fnum = bc_left.ndens_flux * Δt / λ
        n_steps = 2000

        bc_list = (bc_left, bc_right)
        particles = ParticleVector(ceil(Int64, 2.2 * λ * n_steps))
        pia = ParticleIndexerArray(grid.n_cells, 1)
        rng = StableRNG(1234)

        counts = zeros(Int64, n_steps)
        for k in 1:n_steps
            counts[k] = emit_particles!(rng, grid, bc_list, particles, pia, 1, species_data, Fnum, Δt)
        end
        n_emitted = sum(counts)
        @test pia.n_total[1] == n_emitted

        # Poisson-distributed number of particles entering across each of the two boundaries
        mean_count = n_emitted / n_steps
        var_count = sum((counts .- mean_count) .^ 2) / (n_steps - 1)
        @test abs(mean_count - 2 * λ) < 5 * sqrt(2 * λ / n_steps)
        @test abs(var_count / (2 * λ) - 1.0) < 5 * sqrt(2 / n_steps)

        # all particles in left part are coming L->R
        # all particles in right part are coming R->L
        # and didn't travel too far (based on tau)
        τ = zeros(n_emitted)
        inflowing = true
        for i in 1:n_emitted
            p = particles[i]
            if p.x[1] < L / 2
                inflowing &= p.v[1] > 0.0
                τ[i] = (p.x[1] - grid.min_x) / p.v[1]
            else
                inflowing &= p.v[1] < 0.0
                τ[i] = (p.x[1] - grid.max_x) / p.v[1]
            end
        end
        @test inflowing
        @test all(0.0 .< τ .<= Δt * (1.0 + 1e-9))

        # the residual flight times are uniformly distributed over the timestep
        sort!(τ)
        ks = 0.0
        for i in 1:n_emitted
            F = τ[i] / Δt
            ks = max(ks, abs(F - i / n_emitted), abs(F - (i - 1) / n_emitted))
        end
        @test ks < 1.95 / sqrt(n_emitted)

        @test all(particles.cell[i] == Merzbild.get_cell(grid, particles[i].x) for i in 1:n_emitted)
        @test all(particles[i].w == Fnum for i in 1:n_emitted)
        @test check_pia_is_correct(pia, 1) == (true, 0)
        @test check_unique_index(particles, pia, 1) == (true, 0)
        @test check_unique_buffer(particles) == (true, 0)
    end

    @testset "collisionless filling of an empty domain" begin
        L_fill = 1.0
        grid_fill = Grid1DUniform(L_fill, 20)
        Δt_fill = 0.1 * grid_fill.Δx / v_thermal
        n_max = 40000

        for a_target in [0.0, 2.0]
            bc_in = InflowBC1D(1, species_data, T, SVector{3,Float64}(a_target * v_thermal, 0.0, 0.0), normal_left, ndens)
            bc_list = (bc_in, OutflowBC())

            # about 20000 particles in the domain at steady state
            Fnum = bc_in.ndens * (1.0 + erf(bc_in.a)) / 2 * L_fill / 20000
            n_steps = a_target == 0.0 ? 400 : 1000  # 2 L / v_thermal and 5 L / v_thermal

            particles = ParticleVector(n_max)
            pia = ParticleIndexerArray(grid_fill.n_cells, 1)
            gridsorter = GridSortInPlace(grid_fill, n_max)
            surf_props = SurfProps(pia, grid_fill)
            rng = StableRNG(1234)

            flux_emitted = 0.0
            flux_outlet = 0.0
            n_avg_outlet = 0
            for k in 1:n_steps
                convect_particles!(rng, grid_fill, bc_list, particles, pia, 1, species_data, surf_props, Δt_fill)
                emit_particles!(rng, grid_fill, bc_list, particles, pia, 1, species_data, surf_props, Fnum, Δt_fill)
                sort_particles!(gridsorter, grid_fill, particles, pia, 1)

                flux_emitted -= surf_props.flux_reflected[1, 1]
                if k * Δt_fill >= 3 * L_fill / v_thermal
                    flux_outlet += surf_props.flux_incident[2, 1]
                    n_avg_outlet += 1
                end
            end

            n_expected = n_particles_free_flight(bc_in, L_fill, n_steps * Δt_fill, Fnum)
            @test abs(pia.n_total[1] / n_expected - 1.0) < 0.04

            # the emitted particles are recorded as leaving the inlet
            flux_mass = bc_in.ndens_flux * mass
            @test abs(flux_emitted / n_steps / flux_mass - 1.0) < 0.04

            if a_target > 0.0
                # nearly all of the particles are fast enough to have crossed the domain
                @test abs(flux_outlet / n_avg_outlet / flux_mass - 1.0) < 0.04
            end

            @test check_pia_is_correct(pia, 1) == (true, 0)
            @test check_unique_index(particles, pia, 1) == (true, 0)
            @test check_unique_buffer(particles) == (true, 0)
        end
    end

    @testset "equilibrium between two reservoirs" begin
        grid_eq = Grid1DUniform(1.0, 2)
        Δt_eq = 0.2 * grid_eq.Δx / v_thermal
        bc_list = (bc_left, bc_right)

        # about 10000 particles in the domain at steady state
        Fnum = ndens * grid_eq.L / 10000
        n_steps = 2000  # 200 L / v_thermal
        n_avg = 1000

        n_max = 20000
        particles = ParticleVector(n_max)
        pia = ParticleIndexerArray(grid_eq.n_cells, 1)
        gridsorter = GridSortInPlace(grid_eq, n_max)
        rng = StableRNG(1234)

        function velocity_moments(particles, n_p)
            v_sum = SVector{3,Float64}(0.0, 0.0, 0.0)
            v_sq_sum = 0.0
            for i in 1:n_p
                v = particles[i].v
                v_sum += v
                v_sq_sum += v[1]^2 + v[2]^2 + v[3]^2
            end
            return v_sum / n_p, v_sq_sum / n_p
        end

        T_avg = 0.0
        u_avg = 0.0
        for k in 1:n_steps
            convect_particles!(rng, grid_eq, bc_list, particles, pia, 1, species_data, Δt_eq)
            emit_particles!(rng, grid_eq, bc_list, particles, pia, 1, species_data, Fnum, Δt_eq)
            sort_particles!(gridsorter, grid_eq, particles, pia, 1)

            if k > n_steps - n_avg
                v_mean, v_sq_mean = velocity_moments(particles, pia.n_total[1])
                T_avg += mass * (v_sq_mean - sum(v_mean .^ 2)) / (3 * k_B) / n_avg
                u_avg += v_mean[1] / n_avg
            end
        end

        n_expected = 2 * n_particles_free_flight(bc_left, grid_eq.L, n_steps * Δt_eq, Fnum)
        @test abs(pia.n_total[1] / n_expected - 1.0) < 0.05
        @test abs(T_avg / T - 1.0) < 0.015
        @test abs(u_avg) < 0.02 * v_thermal
        @test check_pia_is_correct(pia, 1) == (true, 0)
        @test check_unique_index(particles, pia, 1) == (true, 0)
        @test check_unique_buffer(particles) == (true, 0)
    end

    @testset "appending to an indexer with populated 2nd groups" begin
        # this mimics effect of variable-weight collisions
        # that would populate group2 of indexer in some cells
        n_max = 20000
        particles = ParticleVector(n_max)
        pia = ParticleIndexerArray(grid.n_cells, 1)
        gridsorter = GridSortInPlace(grid, n_max)
        rng = StableRNG(1234)

        Fnum = ndens * grid.cells[1].V / 100
        sample_particles_equal_weight!(rng, grid, particles, pia, 1, species_data, ndens, T, Fnum)
        sort_particles!(gridsorter, grid, particles, pia, 1)

        # particles created in cells 3 and 7 are held in their 2nd groups;
        # this is the same sequence as the particle splitting in ntc!
        for (cell, n_new) in ((3, 4), (7, 2))
            for _ in 1:n_new
                Merzbild.update_buffer_index_new_particle!(particles, pia, cell, 1)
                p = particles[pia.index_last[1]]
                p.w = Fnum
                p.v = SVector{3,Float64}(0.0, 0.0, 0.0)
                p.x = SVector{3,Float64}(grid.cells[cell].xlo + 0.5 * grid.Δx, 0.0, 0.0)
            end
        end
        @test check_pia_is_correct(pia, 1) == (true, 0)

        n_total_before = pia.n_total[1]
        w_total_before = sum(particles[i].w for i in 1:n_total_before)

        Fnum_emit = bc_left.ndens_flux * Δt / 50
        n_emitted = emit_particles!(rng, grid, (bc_left, FullySpecularBC1D()), particles, pia, 1, species_data,
                                    Fnum_emit, Δt)

        @test n_emitted > 0
        @test pia.n_total[1] == n_total_before + n_emitted

        # a 2nd group of cell 1 cannot follow that of cell 7, so the particles are appended to the latter
        @test pia.indexer[1, 1].n_group2 == 0
        @test pia.indexer[7, 1].n_group2 == 2 + n_emitted
        @test check_pia_is_correct(pia, 1) == (true, 0)
        @test check_unique_index(particles, pia, 1) == (true, 0)
        @test check_unique_buffer(particles) == (true, 0)

        sort_particles!(gridsorter, grid, particles, pia, 1)

        @test check_pia_is_correct(pia, 1) == (true, 0)
        @test check_unique_index(particles, pia, 1) == (true, 0)
        @test check_unique_buffer(particles) == (true, 0)
        @test pia.n_total[1] == n_total_before + n_emitted
        @test sum(particles[i].w for i in 1:pia.n_total[1]) ≈ w_total_before + n_emitted * Fnum_emit
    end

    @testset "cell indices of the emitted particles" begin
        # the 3-argument (cell-based) sort has to agree with the 4-argument (position-based) one;
        # this relies on convect_particles_and_compute_cell! and convect_particles! consuming the RNG identically,
        # so that both runs see the same particles and differ only in how the cell indices are obtained
        Fnum = ndens * grid.cells[1].V / 100
        Δt_cells = 100 * Δt
        Fnum_emit = bc_left.ndens_flux * Δt_cells / 50
        bc_list = (bc_left, OutflowBC())

        particles_cells = ParticleVector(20000)
        pia_cells = ParticleIndexerArray(grid.n_cells, 1)
        particles_pos = ParticleVector(20000)
        pia_pos = ParticleIndexerArray(grid.n_cells, 1)
        gridsorter = GridSortInPlace(grid, 20000)

        rng_cells = StableRNG(1234)
        sample_particles_equal_weight!(rng_cells, grid, particles_cells, pia_cells, 1, species_data, ndens, T, Fnum)
        sort_particles!(gridsorter, grid, particles_cells, pia_cells, 1)

        rng_pos = StableRNG(1234)
        sample_particles_equal_weight!(rng_pos, grid, particles_pos, pia_pos, 1, species_data, ndens, T, Fnum)
        sort_particles!(gridsorter, grid, particles_pos, pia_pos, 1)

        for k in 1:5
            convect_particles_and_compute_cell!(rng_cells, grid, bc_list, particles_cells, pia_cells, 1, species_data,
                                                Δt_cells)
            emit_particles!(rng_cells, grid, bc_list, particles_cells, pia_cells, 1, species_data, Fnum_emit, Δt_cells)
            sort_particles!(gridsorter, particles_cells, pia_cells, 1)

            convect_particles!(rng_pos, grid, bc_list, particles_pos, pia_pos, 1, species_data, Δt_cells)
            emit_particles!(rng_pos, grid, bc_list, particles_pos, pia_pos, 1, species_data, Fnum_emit, Δt_cells)
            sort_particles!(gridsorter, grid, particles_pos, pia_pos, 1)
        end

        n_p = pia_cells.n_total[1]
        @test n_p == pia_pos.n_total[1]
        @test particles_cells.index[1:n_p] == particles_pos.index[1:n_p]
        @test check_pia_is_correct(pia_cells, 1) == (true, 0)
        @test check_unique_index(particles_cells, pia_cells, 1) == (true, 0)
        @test check_unique_buffer(particles_cells) == (true, 0)
    end

    @testset "particles leaving the domain within the timestep" begin
        grid_short = Grid1DUniform(1e-3, 2)
        Δt_long = 10 * grid_short.L / v_thermal
        λ = 10000.0
        Fnum = bc_left.ndens_flux * Δt_long / λ

        # the particle vector has to be resized during the emission
        particles = ParticleVector(100)
        pia = ParticleIndexerArray(grid_short.n_cells, 1)
        rng = StableRNG(1234)

        n_emitted = emit_particles!(rng, grid_short, (bc_left, OutflowBC()), particles, pia, 1, species_data, Fnum, Δt_long)

        # a particle with normal velocity v entering with residual flight time τ ~ U(0, Δt) remains in the domain
        # if v τ < L, so the number of retained particles is Poisson-distributed with mean λ E[min(1, L / (v Δt))]
        # over the flux-weighted distribution 2 v exp(-v^2 / v_thermal^2) / v_thermal^2
        s = grid_short.L / (Δt_long * v_thermal)
        n_expected = λ * (1.0 - exp(-s^2) + sqrt(π) * s * erfc(s))
        @test abs(n_emitted - n_expected) < 5 * sqrt(n_expected)
        @test pia.n_total[1] == n_emitted
        @test length(particles) == pia.n_total[1] + particles.nbuffer
        @test all(0.0 < particles[i].x[1] < grid_short.L for i in 1:n_emitted)
        @test check_pia_is_correct(pia, 1) == (true, 0)
        @test check_unique_index(particles, pia, 1) == (true, 0)
        @test check_unique_buffer(particles) == (true, 0)
    end

    @testset "chunks of cells and boundary normals" begin
        Fnum = bc_left.ndens_flux * Δt / 50
        bc_list = (bc_left, bc_right)

        particles = ParticleVector(1000)
        pia = ParticleIndexerArray(grid.n_cells, 1)
        rng = StableRNG(1234)

        @test emit_particles!(rng, grid, bc_list, particles, pia, 1, species_data, Fnum, Δt, 2:nx-1) == 0
        @test pia.n_total[1] == 0

        n_left = emit_particles!(rng, grid, bc_list, particles, pia, 1, species_data, Fnum, Δt, 1:10)
        @test n_left > 0
        @test all(particles[i].v[1] > 0.0 for i in 1:n_left)

        n_right = emit_particles!(rng, grid, bc_list, particles, pia, 1, species_data, Fnum, Δt, 11:nx)
        @test n_right > 0
        @test all(particles[i].v[1] < 0.0 for i in n_left+1:n_left+n_right)

        @test emit_particles!(rng, grid, (FullySpecularBC1D(), OutflowBC()), particles, pia, 1, species_data,
                              Fnum, Δt) == 0

        # an inflow boundary condition placed at the wrong wall
        # current version does not throw
        # @test_throws ArgumentError emit_particles!(rng, grid, (bc_right, OutflowBC()), particles, pia, 1,
        #                                            species_data, Fnum, Δt)
        # @test_throws ArgumentError emit_particles!(rng, grid, (OutflowBC(), bc_left), particles, pia, 1,
        #                                            species_data, Fnum, Δt)
    end

    @testset "surface properties of the emitted particles" begin
        Fnum = bc_left.ndens_flux * Δt / 50
        bc_list = (bc_left, FullySpecularBC1D())

        particles = ParticleVector(1000)
        pia = ParticleIndexerArray(grid.n_cells, 1)
        surf_props = SurfProps(pia, grid)
        rng = StableRNG(1234)

        n_emitted = emit_particles!(rng, grid, bc_list, particles, pia, 1, species_data, surf_props, Fnum, Δt)

        @test n_emitted > 0
        @test isapprox(surf_props.flux_reflected[1, 1], -n_emitted * Fnum * mass / Δt, rtol=1e-12)
        @test surf_props.flux_incident[1, 1] == 0.0
        @test surf_props.np[1, 1] == 0
        @test surf_props.normal_pressure[1, 1] > 0.0
        @test surf_props.kinetic_energy_flux[1, 1] < 0.0

        # the emitted particles do not reach the opposite wall
        @test surf_props.np[2, 1] == 0
        @test surf_props.flux_incident[2, 1] == 0.0
        @test surf_props.flux_reflected[2, 1] == 0.0

        # recording the surface properties does not change the emitted particles
        particles_ref = ParticleVector(1000)
        pia_ref = ParticleIndexerArray(grid.n_cells, 1)
        rng_ref = StableRNG(1234)
        n_emitted_ref = emit_particles!(rng_ref, grid, bc_list, particles_ref, pia_ref, 1, species_data, Fnum, Δt)

        @test n_emitted_ref == n_emitted
        @test all(particles[i].w == Fnum for i in 1:n_emitted)
        @test all(particles[i].v == particles_ref[i].v for i in 1:n_emitted)
        @test all(particles[i].x == particles_ref[i].x for i in 1:n_emitted)

        # in a short domain, emitted particles reach the opposite (specular) wall within the timestep
        grid_short = Grid1DUniform(1e-3, 2)
        Δt_long = 10 * grid_short.L / v_thermal
        Fnum_short = bc_left.ndens_flux * Δt_long / 1000

        particles_short = ParticleVector(2000)
        pia_short = ParticleIndexerArray(grid_short.n_cells, 1)
        surf_props_short = SurfProps(pia_short, grid_short)

        emit_particles!(rng, grid_short, bc_list, particles_short, pia_short, 1, species_data, surf_props_short,
                        Fnum_short, Δt_long)

        # each impact is recorded with the scaled weight of the particle
        @test surf_props_short.np[2, 1] > 0
        @test isapprox(surf_props_short.flux_incident[2, 1],
                       surf_props_short.np[2, 1] * Fnum_short * mass / Δt_long, rtol=1e-12)
        @test isapprox(surf_props_short.flux_reflected[2, 1], -surf_props_short.flux_incident[2, 1], rtol=1e-12)
        @test surf_props_short.normal_pressure[2, 1] > 0.0

        # particles coming back to the inflow wall leave the domain there
        @test surf_props_short.np[1, 1] > 0
        @test all(particles_short[i].w == Fnum_short for i in 1:pia_short.n_total[1])
    end
end
