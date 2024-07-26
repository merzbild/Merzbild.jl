@testset "grid_sampling" begin

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_list::Vector{Species} = load_species_list(particles_data_path, "Ar")

    seed = 1234
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    n_particles = 50000

    n_dens = 1e20
    Δabs = 0.5
    Δrel_small = 1e-2
    Δrel_xsmall = 1e-13
    nv = 20

    Fnum::Float64 = n_dens / n_particles

    for (v0, T0) in zip([[0.0, 0.0, 0.0], [20.0, -10.0, 30.0], [3000.0, 2000.0, -1000.0]], [273.0, 1000.0, 500.0])
        particles::Vector{Vector{Particle}} = [Vector{Particle}(undef, n_particles)]

        n_sampled = sample_maxwellian_on_grid!(rng, particles[1], nv, species_list[1].mass, T0, n_dens,
        0.0, 1.0, 0.0, 1.0, 0.0, 1.0; v_mult=3.5, cutoff_mult=8.0, noise=0.0, v_offset=v0)
        
        pia = create_particle_indexer_array(n_sampled)

        phys_props::PhysProps = create_props(1, 1, [4, 6, 8], Tref=T0)
        compute_props!(particles, pia, species_list, phys_props)
        @test abs((phys_props.n[1,1] - n_dens) / n_dens) < Δrel_xsmall
        @test abs((phys_props.v[1,1,1] - v0[1])) < Δabs
        @test abs((phys_props.v[2,1,1] - v0[2])) < Δabs
        @test abs((phys_props.v[3,1,1] - v0[3])) < Δabs
        @test abs((phys_props.T[1,1] - T0) / T0) < Δrel_small
        @test abs(phys_props.moments[1,1,1] .- 1.0) < 1e-4  # test 4th moment
        @test abs(phys_props.moments[2,1,1] .- 1.0) < 2e-4  # test 6th moment
        @test abs(phys_props.moments[3,1,1] .- 1.0) < 5e-4  # test 8th moment
    end
end