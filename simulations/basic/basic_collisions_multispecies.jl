# Test that sampling two different species from a Maxwellian, doing multiple timesteps of collisions and computing physical properties gives reasonable results

# SPARTA: 12k particles 13.29s user 0.08s system 98% cpu 13.539 total
# this code: 6.90s user 1.78s system 127% cpu 6.829 total

include("../../src/merging.jl")

using ..Merging
using Random

seed = 1234
Random.seed!(seed)
rng = Xoshiro(seed)

species_list = load_species_list("data/particles.toml", ["Ar", "He"])
interaction_data = load_interaction_data("data/vhs.toml", species_list)
n_species = length(species_list)

println([species.name for species in species_list])
println(interaction_data)

n_t = 5000

n_particles_Ar = 2000
n_particles_He = 10000
T0_Ar = 300.0
T0_He = 2000.0
T0_list = [T0_Ar, T0_He]
Fnum = 5e11

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

ds = create_netcdf_phys_props("test_multi_species.nc", phys_props, species_list)
write_netcdf_phys_props(ds, phys_props, 0)

collision_factors = create_collision_factors(n_species)
collision_data = create_collision_data()

estimate_sigma_g_w_max!(collision_factors, interaction_data, species_list, T0_list, Fnum)

Δt = 2.5e-3
V = 1.0

for ts in 1:n_t
 
    for s2 in 1:n_species
        for s1 in s2:n_species
            if (s1 == s2)
                ntc!(rng, collision_factors[s1,s1], particle_indexer[s1,1], collision_data, interaction_data[s1,s1], particles[s1],
                Δt, V)
            else
                ntc!(rng, collision_factors[s1,s2], particle_indexer[s1,1], particle_indexer[s2,1],
                collision_data, interaction_data[s1,s2], particles[s1], particles[s2],
                Δt, V)
            end
        end
    end

    compute_props!(phys_props, particle_indexer, particles, species_list)
    write_netcdf_phys_props(ds, phys_props, ts)
end
println(phys_props.n)
println(phys_props.v)
println(phys_props.T)

close(ds)