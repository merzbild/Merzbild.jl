# Test that sampling from a Maxwellian, doing multiple timesteps of collisions and computing physical properties gives reasonable results

include("../../src/merzbild.jl")

using ..Merzbild
using Random

seed = 1234
Random.seed!(seed)
rng = Xoshiro(seed)

species_list = load_species_list("data/particles.toml", "Ar")
interaction_data = load_interaction_data("data/vhs.toml", species_list)

println([species.name for species in species_list])
println(interaction_data)

n_t = 12

n_particles = 25000
nv = 20
T0 = 500.0
n_dens = 1e23

v0 = [0.0, 0.0, 0.0]

particles = [Vector{Particle}(undef, n_particles)]
n_sampled = sample_maxwellian_on_grid!(rng, particles[1], nv, T0, species_list[1].mass, n_dens,
                                       0.0, 1.0, 0.0, 1.0, 0.0, 1.0;
                                       v_mult=3.5, cutoff_mult=3.5, noise=0.0, v_offset=v0)

particle_indexer = Array{ParticleIndexer, 2}(undef, 1, 1)
particle_indexer[1,1] = create_particle_indexer(n_sampled)

phys_props = create_props(1, 1, [0, 2, 4, 6], Tref=T0)
compute_props!(phys_props, particle_indexer, particles, species_list)
println(phys_props.n)
println(phys_props.v)
println(phys_props.T)

# ds = create_netcdf_phys_props("test_maxwellian.nc", phys_props, species_list)
# write_netcdf_phys_props(ds, phys_props, 0)

collision_factors = create_collision_factors()
collision_data = create_collision_data()

collision_factors.sigma_g_w_max = estimate_sigma_g_w_max(interaction_data[1,1], species_list[1], T0, n_dens/n_sampled)

Δt = 1e-8
V = 1.0

for ts in 1:n_t
    ntc!(rng, collision_factors, particle_indexer[1,1], collision_data, interaction_data[1,1], particles[1],
        Δt, V)

    println(collision_factors)
    
    compute_props!(phys_props, particle_indexer, particles, species_list)
    # write_netcdf_phys_props(ds, phys_props, ts)
    println("------------")
    println(phys_props.np)
    println(phys_props.n)
    println(phys_props.v)
    println(phys_props.T)
end
# println(phys_props.n)
# println(phys_props.v)
# println(phys_props.T)
# close(ds)