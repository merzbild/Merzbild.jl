@testset "phase box sampling" begin

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, "Ar")

    seed = 1234
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)


    n_dens = 1e20
    Δabs = 15.0
    Δrel_small = 3e-2
    Δrel_xsmall = 1e-13
    n_particles = 30000

    for (v0, T0) in zip([[0.0, 0.0, 0.0], [20.0, -10.0, 30.0], [3000.0, 2000.0, -1000.0]], [273.0, 1000.0, 500.0])
        particles::Vector{ParticleVector} = [ParticleVector(n_particles)]

        pia = ParticleIndexerArray(0)

        sample_particles_phase_box_weighted!(rng, particles[1], pia, 1, 1, n_particles, species_data[1].mass, T0, n_dens,
        0.0, 0.5, 0.0, 1.0, 0.0, 2.0; v_mult=3.5, vx0=v0[1], vy0=v0[2], vz0=v0[3])
        
        phys_props::PhysProps = PhysProps(1, 1, [4, 6, 8], Tref=T0)
        compute_props_with_total_moments!(particles, pia, species_data, phys_props)
        @test abs((phys_props.n[1,1] - n_dens) / n_dens) < Δrel_xsmall
        @test abs((phys_props.v[1,1,1] - v0[1])) < Δabs
        @test abs((phys_props.v[2,1,1] - v0[2])) < Δabs
        @test abs((phys_props.v[3,1,1] - v0[3])) < Δabs
        @test abs((phys_props.T[1,1] - T0) / T0) < Δrel_small
        @test abs(phys_props.moments[1,1,1] .- 1.0) < 5e-2  # test 4th moment
        @test abs(phys_props.moments[2,1,1] .- 1.0) < 5e-2  # test 6th moment
        @test abs(phys_props.moments[3,1,1] .- 1.0) < 5e-2  # test 8th moment

        @inbounds pia.indexer[1, 1].start2 == 0
        @inbounds pia.indexer[1, 1].end2 == -1
        @inbounds pia.indexer[1, 1].n_group2 == 0

        oob = false

        for i in 1:n_particles
            if particles[1][i].x[1] < 0.0 || particles[1][i].x[1] > 0.5
                oob = true
                break
            end
            if particles[1][i].x[2] < 0.0 || particles[1][i].x[2] > 1.0
                oob = true
                break
            end
            if particles[1][i].x[3] < 0.0 || particles[1][i].x[3] > 2.0
                oob = true
                break
            end
        end
        @test oob == false
    end
end