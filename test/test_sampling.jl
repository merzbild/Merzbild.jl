@testset "sampling" begin

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, "Ar")

    seed = 1234
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    n_particles = 20000

    n_dens = 1e20
    Δlarge = 10.0
    Δsmall = 1e-2

    Fnum::Float64 = n_dens / n_particles

    for (v0, T0) in zip([[0.0, 0.0, 0.0], [20.0, -10.0, 30.0], [3000.0, 2000.0, -1000.0]], [273.0, 1000.0, 500.0])
        particles::Vector{ParticleVector} = [ParticleVector(n_particles)]

        pia = ParticleIndexerArray(0)

        # particles, pia, cell, species,
        # nparticles, m, T, Fnum, xlo, xhi, ylo, yhi, zlo, zhi;
        # distribution=:Maxwellian, vx0=0.0, vy0=0.0, vz0=0.0)

        sample_particles_equal_weight!(rng, particles[1], pia, 1, 1, n_particles,
                                       species_data[1].mass, T0, Fnum, 0.0, 0.5, 0.0, 1.0, 0.0, 2.0,
                                       vx0=v0[1], vy0=v0[2], vz0=v0[3])

        phys_props::PhysProps = PhysProps(1, 1)
        compute_props!(particles, pia, species_data, phys_props)
        @test abs((phys_props.n[1,1] - n_dens) / n_dens) < eps()
        @test abs((phys_props.v[1,1,1] - v0[1])) < Δlarge
        @test abs((phys_props.v[2,1,1] - v0[2])) < Δlarge
        @test abs((phys_props.v[3,1,1] - v0[3])) < Δlarge
        @test abs((phys_props.T[1,1] - T0) / T0) < Δsmall

        mlist = [4,6,8]
        mscaling = zeros(3)
        mvals_list = zeros(3)
        compute_moment_scaling!(mscaling, mlist, 1, species_data, T0)
        compute_moments!(mvals_list, mscaling, mlist, particles, pia, 1, 1, species_data, phys_props)

        @test abs(mvals_list[1] - 1.0) < 0.05  # test 4th moment
        @test abs(mvals_list[2] - 1.0) < 0.05  # test 6th moment
        @test abs(mvals_list[3] - 1.0) < 0.12  # test 8th moment

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