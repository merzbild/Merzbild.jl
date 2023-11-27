include("../src/merging.jl")

using .Merging
using Random

species_list = load_species_list("data/particles.toml", "Ar")

rng = Xoshiro(1234)

n_particles = 5000

particles = Vector{Particle}(undef, n_particles)
sample_particles_equal_weight!(rng, particles, n_particles, 500.0, species_list[1].mass, 1e10, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)

particle_indexer = Array{ParticleIndexer, 2}(undef, 1, 1)
particle_indexer[1,1] = create_particle_indexer(n_particles)

phys_props = create_props(1, 1, [])
compute_props!(phys_props, particle_indexer, particles, 1, 1, species_list)
println(phys_props.n)
println(phys_props.v)
println(phys_props.T)