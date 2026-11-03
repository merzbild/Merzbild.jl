@testset "inflow sampling" begin

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, "Ar")

    mass = species_data[1].mass
    T = 300.0
    v_thermal = Merzbild.compute_thermal_velocity(mass, T)
    ndens = 1e20

    normal_left = SVector{3,Float64}(1.0, 0.0, 0.0)
    normal_right = SVector{3,Float64}(-1.0, 0.0, 0.0)

    m_a(a) = exp(-a^2) + a * sqrt(π) * (1.0 + erf(a))  # eqn 3.2 in Garcia Wagner

    # Garcia and Wagner, JCP 217 (2006), section 4
    F_a(z, a) = (exp(-z^2) + a * sqrt(π) * (1.0 + erf(z))) / m_a(a)

    # Garcia and Wagner, appendix A, with the paper's T equal to v_thermal^2 / 2
    mean_v_normal(a) = (v_thermal / 2) * (sqrt(π) * (1.0 + erf(a)) * (1.0 + 2 * a^2) + 2 * a * exp(-a^2)) / m_a(a)

    make_bc(a_target, normal) = InflowBC1D(1, species_data, T,
                                           SVector{3,Float64}(a_target * v_thermal * normal[1], 0.0, 0.0),
                                           normal, ndens)

    @testset "constructor and sampler selection" begin
        bc_no_drift = make_bc(0.0, normal_left)
        bc_drift = make_bc(1.0, normal_left)

        @test bc_no_drift.sampler === InflowNoDrift()
        @test bc_drift.sampler === InflowEnvelope4()

        @test bc_no_drift.a == 0.0
        @test abs(bc_drift.a - 1.0) < 1e-12

        # a streaming velocity tangential to the surface gives a zero speed ratio
        bc_tangential = InflowBC1D(1, species_data, T, SVector{3,Float64}(0.0, 500.0, -300.0), normal_left, ndens)
        @test bc_tangential.sampler === InflowNoDrift()

        # the right wall has an inward normal of (-1, 0, 0), so an inflow drifts in -x
        bc_right = make_bc(1.5, normal_right)
        @test bc_right.sampler === InflowEnvelope4()
        @test abs(bc_right.a - 1.5) < 1e-12

        # streaming velocity pointing out of the domain
        @test_throws ArgumentError InflowBC1D(1, species_data, T, SVector{3,Float64}(-500.0, 0.0, 0.0),
                                              normal_left, ndens)
        # explicitly requesting a sampler that does not apply
        @test_throws ArgumentError InflowBC1D(InflowNoDrift(), 1, species_data, T,
                                              SVector{3,Float64}(500.0, 0.0, 0.0), normal_left, ndens)
        @test_throws ArgumentError InflowBC1D(InflowEnvelope4(), 1, species_data, T,
                                              SVector{3,Float64}(-500.0, 0.0, 0.0), normal_left, ndens)

        # Envelope 4 requires a positive speed ratio, otherwise the acceptance-rejection loop never terminates
        @test_throws ArgumentError InflowBC1D(InflowEnvelope4(), 1, species_data, T,
                                              SVector{3,Float64}(0.0, 0.0, 0.0), normal_left, ndens)
        @test_throws ArgumentError InflowBC1D(InflowEnvelope4(), 1, species_data, T,
                                              SVector{3,Float64}(0.0, 500.0, -300.0), normal_left, ndens)

        # a zero temperature gives a non-finite speed ratio
        @test_throws ArgumentError InflowBC1D(1, species_data, 0.0, SVector{3,Float64}(500.0, 0.0, 0.0),
                                              normal_left, ndens)
        @test_throws ArgumentError InflowBC1D(1, species_data, 0.0, SVector{3,Float64}(0.0, 0.0, 0.0),
                                              normal_left, ndens)

        # the sampler assumes an inward normal of (1, 0, 0) or (-1, 0, 0)
        for normal in (SVector{3,Float64}(0.0, 1.0, 0.0), SVector{3,Float64}(2.0, 0.0, 0.0),
                       SVector{3,Float64}(1.0, 0.0, 1e-3))
            @test_throws ArgumentError InflowBC1D(1, species_data, T, SVector{3,Float64}(0.0, 0.0, 0.0), normal, ndens)
            @test_throws ArgumentError InflowBC1D(InflowEnvelope4(), 1, species_data, T,
                                                  SVector{3,Float64}(500.0, 0.0, 0.0), normal, ndens)
        end

        # v and normal are converted to SVector{3,Float64} before use
        bc_from_vectors = InflowBC1D(1, species_data, T, [v_thermal, 0.0, 0.0], [1.0, 0.0, 0.0], ndens)
        @test abs(bc_from_vectors.a - 1.0) < 1e-12

        # the field-wise constructor would bypass validation
        @test_throws MethodError InflowBC1D(InflowEnvelope4(), T, SVector{3,Float64}(0.0, 0.0, 0.0), normal_left,
                                            ndens, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    end

    @testset "number flux" begin
        # reservoir at rest: effusion flux n * sqrt(kT / (2 pi m))
        bc = make_bc(0.0, normal_left)
        @test abs(bc.ndens_flux / (ndens * sqrt(Merzbild.k_B * T / (Merzbild.twopi * mass))) - 1.0) < 1e-14

        # the exported utility and the pre-computed field must agree
        @test compute_inflow_number_flux(1, species_data, T, SVector{3,Float64}(0.0, 0.0, 0.0),
                                         normal_left, ndens) == bc.ndens_flux

        # flux across a surface is unaffected by a purely tangential streaming velocity
        @test abs(compute_inflow_number_flux(1, species_data, T, SVector{3,Float64}(0.0, 1000.0, 500.0),
                                             normal_left, ndens) / bc.ndens_flux - 1.0) < 1e-14

        # compare against numerical quadrature of int_{v_n > 0} v_n * M(v) dv
        for a_target in [0.0, 0.5, 2.0]
            v_drift = a_target * v_thermal
            f(vn) = vn * exp(-(vn - v_drift)^2 / v_thermal^2)

            n_quad = 400000
            h = (v_drift + 12 * v_thermal) / n_quad
            integral = 0.5 * (f(0.0) + f(v_drift + 12 * v_thermal))
            for i in 1:n_quad-1
                integral += f(i * h)
            end
            quad_flux = ndens * integral * h / (v_thermal * sqrt(π))

            flux = compute_inflow_number_flux(1, species_data, T, SVector{3,Float64}(v_drift, 0.0, 0.0),
                                              normal_left, ndens)
            @test abs(flux / quad_flux - 1.0) < 1e-8
        end

        # the flux increases monotonically with the speed ratio
        fluxes = [compute_inflow_number_flux(1, species_data, T,
                                             SVector{3,Float64}(a_target * v_thermal, 0.0, 0.0),
                                             normal_left, ndens) for a_target in 0.0:0.5:5.0]
        @test all(diff(fluxes) .> 0.0)
    end

    n_samples = 200000

    # two-sided Kolmogorov-Smirnov critical value at the 99.9% level
    ks_threshold = 1.95 / sqrt(n_samples)

    @testset "distribution of z" begin
        for a_target in [0.0, 0.25, 1.0, 3.0, 5.0]
            bc = make_bc(a_target, normal_left)
            a = bc.a

            rng::Xoshiro = Xoshiro(1234)
            z = Vector{Float64}(undef, n_samples)
            for i in 1:n_samples
                z[i] = Merzbild.sample_inflow_z(rng, bc.sampler, bc)
            end

            # the distribution has support (-Inf, a), which guarantees an inflowing normal velocity
            @test all(z .< a)

            sort!(z)
            ks = 0.0
            for i in 1:n_samples
                F = F_a(z[i], a)
                ks = max(ks, abs(F - i / n_samples), abs(F - (i - 1) / n_samples))
            end
            @test ks < ks_threshold
        end
    end

    @testset "Envelope 4 in the zero speed ratio limit" begin
        # Envelope 4 at a vanishing speed ratio must reproduce the a = 0 distribution
        bc = InflowBC1D(InflowEnvelope4(), 1, species_data, T,
                        SVector{3,Float64}(1e-9 * v_thermal, 0.0, 0.0), normal_left, ndens)

        rng::Xoshiro = Xoshiro(5678)
        z = Vector{Float64}(undef, n_samples)
        for i in 1:n_samples
            z[i] = Merzbild.sample_inflow_z(rng, bc.sampler, bc)
        end

        sort!(z)
        ks = 0.0
        for i in 1:n_samples
            F = F_a(z[i], 0.0)
            ks = max(ks, abs(F - i / n_samples), abs(F - (i - 1) / n_samples))
        end
        @test ks < ks_threshold
    end

    @testset "sampled velocity moments" begin
        v_tangential = SVector{3,Float64}(0.0, 200.0, -100.0)

        for (a_target, normal) in zip([0.0, 0.75, 2.0, 1.25], [normal_left, normal_left, normal_left, normal_right])
            v_drift = SVector{3,Float64}(a_target * v_thermal * normal[1], v_tangential[2], v_tangential[3])
            bc = InflowBC1D(1, species_data, T, v_drift, normal, ndens)

            rng::Xoshiro = Xoshiro(9012)
            sum_vn = 0.0
            sum_vy = 0.0
            sum_vz = 0.0
            sum_vy_sq = 0.0
            sum_vz_sq = 0.0
            all_inflowing = true

            for _ in 1:n_samples
                v = sample_inflow_velocity(rng, bc)
                vn = v[1] * normal[1]
                all_inflowing &= vn > 0.0

                sum_vn += vn
                sum_vy += v[2]
                sum_vz += v[3]
                sum_vy_sq += (v[2] - v_tangential[2])^2
                sum_vz_sq += (v[3] - v_tangential[3])^2
            end

            # every sampled particle moves into the domain
            @test all_inflowing

            @test abs((sum_vn / n_samples) / mean_v_normal(bc.a) - 1.0) < 1e-2

            # the tangential components are Gaussians centred on the tangential streaming velocity
            # with a variance of v_thermal^2 / 2
            @test abs(sum_vy / n_samples - v_tangential[2]) < 0.02 * v_thermal
            @test abs(sum_vz / n_samples - v_tangential[3]) < 0.02 * v_thermal
            @test abs((sum_vy_sq / n_samples) / (0.5 * v_thermal^2) - 1.0) < 1e-2
            @test abs((sum_vz_sq / n_samples) / (0.5 * v_thermal^2) - 1.0) < 1e-2
        end
    end

    @testset "inflow boundaries delete outgoing particles" begin
        bc_inflow = make_bc(1.0, normal_left)

        @test Merzbild.has_particle_deleting_bc((bc_inflow, FullySpecularBC1D()))
        @test Merzbild.bc_deletes_particles((bc_inflow, FullySpecularBC1D()), 1)
        @test !Merzbild.bc_deletes_particles((bc_inflow, FullySpecularBC1D()), 2)

        # the no-drift sampler parameterization must be recognized as well
        @test Merzbild.has_particle_deleting_bc((FullySpecularBC1D(), make_bc(0.0, normal_right)))

        # boundary condition tuples without inflow/outflow boundaries are unaffected
        @test !Merzbild.has_particle_deleting_bc((FullySpecularBC1D(), FullySpecularBC1D()))

        # applying a particle-deleting boundary condition leaves the particle and the remaining time unchanged
        rng_bc = Xoshiro(1234)
        v0 = SVector{3,Float64}(-150.0, 20.0, -5.0)
        x0 = SVector{3,Float64}(0.0, 1.0, 2.0)
        for (bc_id, bc_list) in ((1, (bc_inflow, FullySpecularBC1D())), (2, (FullySpecularBC1D(), OutflowBC())))
            normal = bc_id == 1 ? normal_left : normal_right

            particle = Particle(2.5, v0, x0)
            @test Merzbild.apply_bc!(rng_bc, particle, bc_list[bc_id], normal, 0.75) == 0.75
            @test particle.v == v0
            @test particle.x == x0
            @test particle.w == 2.5

            @test Merzbild.apply_bc_dispatched!(rng_bc, particle, bc_list, bc_id, normal, 0.75) == 0.75
            @test particle.v == v0
            @test particle.x == x0
            @test particle.w == 2.5
        end
    end
end
