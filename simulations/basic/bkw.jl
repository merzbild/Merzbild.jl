# BKW test case
# 10k particles, dt_scaled = 0.025, n_t = 500 julia --project=. simulations/basic/bkw.jl  4.49s user 1.71s system 134% cpu 4.595 total
# 20k particles, julia --project=. simulations/basic/bkw.jl  4.73s user 1.75s system 133% cpu 4.874 total
# 50k particles, no print julia --project=. simulations/basic/bkw.jl  5.54s user 1.65s system 134% cpu 5.357 total
# 200k particles, 10 moments, no print julia --project=. simulations/basic/bkw.jl  13.29s user 1.78s system 114% cpu 13.211 total

include("../../src/merging.jl")

using ..Merging
using Random

const seed = 1234
Random.seed!(seed)
rng::Xoshiro = Xoshiro(seed)

const species_list::Vector{Species} = load_species_list("data/particles.toml", "Ar")
const interaction_data::Array{Interaction, 2} = load_interaction_data("data/pseudo_maxwell.toml", species_list)

println([species.name for species in species_list])
println(interaction_data)

# Important!
# The time scaling in the analytical solution is different
# Tref = 273.0
# mref = 66.3e-27 
# mcd = mref / 2.0
# dref = 4.11e-10
# nref = 1e23
# Lref = 1.0 / (nref * constants.pi * dref**2)
# vref = ((2 * constants.k * Tref) / mref)**0.5
# time_ref = Lref / vref

# kappa_mult = constants.pi * dref**2 * (mcd / (2 * constants.k * tref))**(-0.5) / gamma(5/2 - 1.0)
# ttt_bkw = 1 / (4 * constants.pi * n * kappa_mult)
# magic_factor = time_ref / ttt_bkw / (4 * constants.pi)
# print(magic_factor)  # approximately 1.59577 for Argon, 1.5963 for N
    
# def analytic(time, N):    
#     C = 1. - 0.4 * np.exp(-time * magic_factor / 6)
#     kk = N // 2
#     return C**(kk - 1) * (kk - (kk - 1) * C)

const dt_scaled = 0.025
const n_t = 500
const n_particles = 500000

const T0::Float64 = 273.0
const sigma_ref = π * (interaction_data[1,1].vhs_d^2)
const n_dens = 1e23

const vref = sqrt(2 * k_B * T0 / species_list[1].mass)
const Lref = 1.0 / (n_dens * sigma_ref)
const tref = Lref / vref

const Fnum::Float64 = n_dens / n_particles

particles::Vector{Vector{Particle}} = [Vector{Particle}(undef, n_particles)]
sample_particles_equal_weight!(rng, particles[1], n_particles, T0, species_list[1].mass, Fnum,
                               0.0, 1.0, 0.0, 1.0, 0.0, 1.0; distribution=:BKW)

particle_indexer::Array{ParticleIndexer, 2} = Array{ParticleIndexer, 2}(undef, 1, 1)
particle_indexer[1,1] = create_particle_indexer(n_particles)

phys_props::PhysProps = create_props(1, 1, [4, 6, 8, 10], Tref=T0)
compute_props!(phys_props, particle_indexer, particles, species_list)
println(phys_props.n)
println(phys_props.v)
println(phys_props.T)

ds = create_netcdf_phys_props("bkw.nc", phys_props, species_list)
write_netcdf_phys_props(ds, phys_props, 0)

collision_factors::CollisionFactors = create_collision_factors()
collision_data::CollisionData = create_collision_data()

collision_factors.sigma_g_w_max = estimate_sigma_g_w_max(interaction_data[1,1], species_list[1], T0, Fnum)

const Δt::Float64 = dt_scaled * tref
const V::Float64 = 1.0

@time for ts in 1:n_t
    ntc!(rng, collision_factors, particle_indexer[1,1], collision_data, interaction_data[1,1], particles[1],
        Δt, V)

    # println(collision_factors)
    
    compute_props!(phys_props, particle_indexer, particles, species_list)
    write_netcdf_phys_props(ds, phys_props, ts)
end
close(ds)