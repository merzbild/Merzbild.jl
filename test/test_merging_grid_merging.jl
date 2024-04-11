@testset "merging_grid_merging" begin

    species_list::Vector{Species} = load_species_list("data/particles.toml", "Ar")

    seed = 1234
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    Nx = 2
    Ny = 2
    Nz = 2

    phys_props::PhysProps = create_props(1, 1, [], Tref=1)

    Δabs = 2.5
    Δrel_xsmall = 5e-13

    n_particles = 40000
    particles::Vector{Vector{Particle}} = [Vector{Particle}(undef, n_particles)]
    T0 = 300.0
    Fnum = 1e20
    vx0 = 2000.0
    vy0 = 500.0
    vz0 = -400.0
    sample_particles_equal_weight!(rng, particles[1], n_particles, T0, species_list[1].mass, Fnum, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
    vx0=vx0, vy0=vy0, vz0=vz0)

    particle_indexer = create_particle_indexer_array(n_particles)

    compute_props!(phys_props, particle_indexer, particles, species_list)
    mg = create_merging_grid(Nx, Ny, Nz, 3.5)

    Merzbild.compute_velocity_extent!(1, 1, mg, phys_props, species_list)
    merge_grid_based!(1, 1, mg, phys_props, species_list, particles, particle_indexer)
end