# Test that sampling two different species from a Maxwellian, doing multiple timesteps of collisions and computing physical properties gives reasonable results

include("../../src/merging.jl")

using ..Merging
using Random

species_list = load_species_list("data/particles.toml", ["Ar", "He"])
interaction_data = load_interaction_data("data/vhs.toml", species_list)

println([species.name for species in species_list])
println(interaction_data)
rng = Xoshiro(1234)

# n_t = 1000

n_particles_Ar = 5000
n_particles_He = 1000
T0_Ar = 1000.0
T0_He = 300.0
Fnum = 1e11

particles = [Vector{Particle}(undef, n_particles_Ar), Vector{Particle}(undef, n_particles_He)]
sample_particles_equal_weight!(rng, particles[1], n_particles_Ar, T0_Ar, species_list[1].mass, Fnum, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
sample_particles_equal_weight!(rng, particles[2], n_particles_He, T0_He, species_list[2].mass, Fnum, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)

# particle_indexer_array[species,cell]
particle_indexer = Array{ParticleIndexer, 2}(undef, 2, 1)
particle_indexer[1,1] = create_particle_indexer(n_particles_Ar)
particle_indexer[2,1] = create_particle_indexer(n_particles_He)

phys_props = create_props(1, 2, [0, 2, 4, 6], Tref=T0_Ar)
compute_props!(phys_props, particle_indexer, particles, species_list)
println(phys_props.n)
println(phys_props.v)
println(phys_props.T)

# ds = create_netcdf_phys_props("test_maxwellian.nc", phys_props, species_list)
# write_netcdf_phys_props(ds, phys_props, 0)

# collision_factors = create_collision_factors()
# collision_data = create_collision_data()

# collision_factors.sigma_g_w_max = estimate_sigma_g_w_max(interaction_data[1,1], species_list[1], T0, Fnum)

# Δt = 1e-6
# V = 1.0

# for ts in 1:n_t
#     ntc_single_species!(rng, collision_factors, particle_indexer[1,1], collision_data, interaction_data[1,1], particles,
#         Δt, V)

#     print(collision_factors)
    
#     compute_props!(phys_props, particle_indexer, particles, species_list)
#     write_netcdf_phys_props(ds, phys_props, ts)
# end
# # println(phys_props.n)
# # println(phys_props.v)
# # println(phys_props.T)
# close(ds)