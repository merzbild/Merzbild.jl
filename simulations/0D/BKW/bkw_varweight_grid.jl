# BKW with grid merging, also with @code_warntype to check for type instabilities
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
include("../../../src/Merzbild.jl")

using ..Merzbild
using Random
using InteractiveUtils

function run(seed)
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    species_data::Vector{Species} = load_species_data("data/particles.toml", "Ar")
    interaction_data::Array{Interaction, 2} = load_interaction_data("data/pseudo_maxwell.toml", species_data)

    dt_scaled = 0.025
    n_t = 500

    nv = 40
    np_base = 40^3  # some initial guess on # of particle in simulation

    threshold = 10000
    Nmerging = 16  # ~8000 after merging

    mg = GridN2Merge(Nmerging, Nmerging, Nmerging, 3.5)

    T0::Float64 = 273.0
    sigma_ref = π * (interaction_data[1,1].vhs_d^2)
    n_dens = 1e23

    vref = sqrt(2 * k_B * T0 / species_data[1].mass)
    Lref = 1.0 / (n_dens * sigma_ref)
    tref = Lref / vref
    moments_list = [4, 6, 8, 10]

    particles::Vector{ParticleVector} = [ParticleVector(np_base)]

    vdf0 = (vx, vy, vz) -> bkw(vx, vy, vz, species_data[1].mass, T0, 0.0)

    n_sampled = sample_on_grid!(rng, vdf0, particles[1], nv, species_data[1].mass, T0, n_dens,
                                0.0, 1.0, 0.0, 1.0, 0.0, 1.0;
                                v_mult=3.5, cutoff_mult=3.5, noise=0.0, v_offset=[0.0, 0.0, 0.0])
    # println(n_sampled)

    pia = ParticleIndexerArray(n_sampled)

    phys_props::PhysProps = PhysProps(1, 1, moments_list, Tref=T0)
    compute_props!(particles, pia, species_data, phys_props)

    ds = NCDataHolder("scratch/data/tmp_bkw_grid.nc", species_data, phys_props)
    write_netcdf_phys_props(ds, phys_props, 0)

    collision_factors::CollisionFactors = CollisionFactors()
    collision_data::CollisionData = CollisionData()

    Fnum = n_dens/n_sampled
    collision_factors.sigma_g_w_max = estimate_sigma_g_w_max(interaction_data[1,1], species_data[1], T0, Fnum)

    Δt::Float64 = dt_scaled * tref
    V::Float64 = 1.0

    for ts in 1:n_t
        ntc!(rng, collision_factors, collision_data, interaction_data, particles[1], pia, 1, 1, Δt, V)

        if phys_props.np[1,1] > threshold
            merge_grid_based!(rng, mg, particles[1], pia, 1, 1, species_data, phys_props)
        end
        
        compute_props!(particles, pia, species_data, phys_props)
        write_netcdf_phys_props(ds, phys_props, ts)
    end
    close_netcdf(ds)
end
@code_warntype run(1234)