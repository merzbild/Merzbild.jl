@testset "merging_grid_merging" begin

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_list::Vector{Species} = load_species_list(particles_data_path, "Ar")

    seed = 1234
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    Nx = 2
    Ny = 2
    Nz = 2

    phys_props::PhysProps = create_props(1, 1, [], Tref=1)

    Δabs = 2.5
    Δrel_xsmall = 5e-13

    n_particles = 5000
    particles::Vector{Vector{Particle}} = [Vector{Particle}(undef, n_particles)]
    T0 = 300.0
    Fnum = 1e8
    vx0 = 2000.0
    vy0 = 500.0
    vz0 = -400.0
    sample_particles_equal_weight!(rng, particles[1], n_particles, species_list[1].mass, T0, Fnum, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
    vx0=vx0, vy0=vy0, vz0=vz0)

    pia = create_particle_indexer_array(n_particles)

    compute_props!(particles, pia, species_list, phys_props)

    n0_computed = phys_props.n[1,1]
    np0_computed = phys_props.np[1,1]
    T0_computed = phys_props.T[1,1]
    v0_computed = phys_props.v[:,1,1]

    mg = create_merging_grid(Nx, Ny, Nz, 3.5)

    merge_grid_based!(mg, particles[1], pia, 1, 1, species_list, phys_props)

    @test pia.n_total[1] < n_particles
    @test pia.n_total[1] == pia.indexer[1,1].n_local

    compute_props!(particles, pia, species_list, phys_props)
    @test pia.n_total[1] == phys_props.np[1,1]
    @test abs(phys_props.n[1,1] - n0_computed) < eps()
    @test abs(phys_props.T[1,1] - T0_computed) < 3e-12
    @test abs(v0_computed[1] - phys_props.v[1,1,1]) < 3e-12
    @test abs(v0_computed[2] - phys_props.v[2,1,1]) < 3e-12
    @test abs(v0_computed[3] - phys_props.v[3,1,1]) < 3e-12
end