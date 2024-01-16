@testset "sampling" begin

    species_list::Vector{Species} = load_species_list("data/particles.toml", "Ar")

    seed = 1234
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    n_particles = 20000

    n_dens = 1e20
    Δlarge = 10.0
    Δsmall = 1e-2

    Fnum::Float64 = n_dens / n_particles

    for (v0, T0) in zip([[0.0, 0.0, 0.0], [20.0, -10.0, 30.0], [3000.0, 2000.0, -1000.0]], [273.0, 1000.0, 500.0])
        particles::Vector{Vector{Particle}} = [Vector{Particle}(undef, n_particles)]

        sample_particles_equal_weight!(rng, particles[1], n_particles, T0, species_list[1].mass, Fnum, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
        vx0=v0[1], vy0=v0[2], vz0=v0[3])

        particle_indexer::Array{ParticleIndexer, 2} = Array{ParticleIndexer, 2}(undef, 1, 1)
        particle_indexer[1,1] = create_particle_indexer(n_particles)

        phys_props::PhysProps = create_props(1, 1, [], Tref=1)
        phys_props_no_moments::PhysProps = create_props(1, 1, [], Tref=1)
        compute_props!(phys_props, particle_indexer, particles, species_list)
        @test abs((phys_props.n[1,1] - n_dens) / n_dens) < eps()
        @test abs((phys_props.v[1,1,1] - v0[1])) < Δlarge
        @test abs((phys_props.v[2,1,1] - v0[2])) < Δlarge
        @test abs((phys_props.v[3,1,1] - v0[3])) < Δlarge
        @test abs((phys_props.T[1,1] - T0) / T0) < Δsmall
    end
end