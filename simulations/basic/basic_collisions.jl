# Test that sampling from a Maxwellian, doing 1 timestep of collisions and computing physical properties gives reasonable results

include("../../src/merzbild.jl")

using ..Merzbild
using Random

species_list = load_species_list("data/particles.toml", "Ar")
interaction_data = load_interaction_data("data/vhs.toml", species_list)

println([species.name for species in species_list])
println(interaction_data)
rng = Xoshiro(1234)

n_particles = 5000
T0 = 500.0
Fnum = 1e10

particles = [Vector{Particle}(undef, n_particles)]
sample_particles_equal_weight!(rng, particles[1], n_particles, T0, species_list[1].mass, Fnum, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)

particle_indexer = Array{ParticleIndexer, 2}(undef, 1, 1)
particle_indexer[1,1] = create_particle_indexer(n_particles)

phys_props = create_props(1, 1, [])
compute_props!(phys_props, particle_indexer, particles, species_list)
println(phys_props.n)
println(phys_props.v)
println(phys_props.T)

# ds = create_netcdf_phys_props("test_collisions.nc", phys_props, species_list)
# write_netcdf_phys_props(ds, phys_props, 1)

collision_factors = create_collision_factors()
collision_data = create_collision_data()

collision_factors.sigma_g_w_max = estimate_sigma_g_w_max(interaction_data[1,1], species_list[1], T0, Fnum)

Δt = 1e-6
V = 1.0
ntc!(rng, collision_factors, particle_indexer[1,1], collision_data, interaction_data[1,1], particles[1],
    Δt, V)

print(collision_factors)

compute_props!(phys_props, particle_indexer, particles, species_list)
println(phys_props.n)
println(phys_props.v)
println(phys_props.T)
# close(ds)