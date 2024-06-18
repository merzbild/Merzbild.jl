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
include("../../src/merzbild.jl")

using ..Merzbild
using Random
using InteractiveUtils
using TimerOutputs

function run(seed)
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    reset_timer!()
    species_list::Vector{Species} = load_species_list("data/particles.toml", "Ar")
    interaction_data::Array{Interaction, 2} = load_interaction_data("data/pseudo_maxwell.toml", species_list)

    dt_scaled = 0.025
    n_t = 500

    nv = 30
    np_base = 40^3  # some initial guess on # of particle in simulation
    n_moms = 6

    mim = []
    for i in 1:n_moms
        append!(mim, compute_multi_index_moments(i))
    end
    # append!(mim, [[4, 0, 0], [0, 4, 0], [0, 0, 4]])
    append!(mim, [[7, 0, 0], [0, 7, 0], [0, 0, 7]])
    append!(mim, [[8, 0, 0], [0, 8, 0], [0, 0, 8]])
    # append!(mim, [[10, 0, 0], [0, 10, 0], [0, 0, 10]])
    # append!(mim, )
    # append!(mim, )
    println(length(mim))

    threshold = 200
    ntarget = 300
    ntarget_octree = length(mim) + 16

    # nbins * length(mim) ~ ntarget
    # Nmerging = 16  # ~8000 after merging

    @timeit "NNLSinit" mnnls = create_nnls_merging(mim, threshold)
    @timeit "OctreeInit" ocm = create_merging_octree(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel, max_Nbins=6000)

    T0::Float64 = 273.0
    sigma_ref = π * (interaction_data[1,1].vhs_d^2)
    n_dens = 1e23

    vref = sqrt(2 * k_B * T0 / species_list[1].mass)
    Lref = 1.0 / (n_dens * sigma_ref)
    tref = Lref / vref
    moments_list = [4, 6, 8, 10]

    particles::Vector{Vector{Particle}} = [Vector{Particle}(undef, np_base)]

    vdf0 = (vx, vy, vz) -> bkw(T0, 0.0, species_list[1].mass, vx, vy, vz)

    n_sampled = sample_on_grid!(rng, vdf0, particles[1], nv, T0, species_list[1].mass, n_dens,
                                0.0, 1.0, 0.0, 1.0, 0.0, 1.0;
                                v_mult=3.5, cutoff_mult=3.5, noise=0.0, v_offset=[0.0, 0.0, 0.0])
    # println(n_sampled)

    pia = create_particle_indexer_array(n_sampled)

    phys_props::PhysProps = create_props(1, 1, moments_list, Tref=T0)
    compute_props!(phys_props, pia, particles, species_list)
    
    # if phys_props.np[1,1] > threshold
    #     @timeit "t=0 NNLSmerge" merge_nnls_based!(rng, mnnls, vref, particles, 1, 1, pia, 0)
    #     compute_props!(phys_props, pia, particles, species_list)
    # end

    ds = create_netcdf_phys_props("test.nc",phys_props, species_list)
    write_netcdf_phys_props(ds, phys_props, 0)

    collision_factors::CollisionFactors = create_collision_factors()
    collision_data::CollisionData = create_collision_data()

    Fnum = n_dens/n_sampled
    collision_factors.sigma_g_w_max = estimate_sigma_g_w_max(interaction_data[1,1], species_list[1], T0, Fnum)

    Δt::Float64 = dt_scaled * tref
    V::Float64 = 1.0

    firstm = true
    nnls_success_flag = 1

    for ts in 1:n_t
        if ts % 10 == 0
            println(ts)
        end

        ntc!(1, 1, rng, collision_factors, pia, collision_data, interaction_data[1,1], particles[1],
            Δt, V)

        if phys_props.np[1,1] > threshold
            if firstm
                @timeit "NNLSmerge: 1st time" nnls_success_flag = merge_nnls_based!(rng, mnnls, vref, particles, 1, 1, pia, 10)
                firstm = false
            else
                @timeit "NNLSmerge" nnls_success_flag = merge_nnls_based!(rng, mnnls, vref, particles, 1, 1, pia, 10)
            end

            if nnls_success_flag == -1
                println("Resorting to octree merging")
                @timeit "Octreemerge" merge_octree_N2_based!(1, 1, ocm, particles, pia, ntarget_octree)
            end
        end
        
        compute_props!(phys_props, pia, particles, species_list)
        write_netcdf_phys_props(ds, phys_props, ts)
    end
    close(ds)
    print_timer()
end
# @code_warntype run(1234)
run(1234)

run(1234)   # 30^3, 400 threshold, 6 moments, 0.965364 seconds (528.50 k allocations: 53.510 MiB, 1.01% gc time)
# 0.968882 seconds (528.28 k allocations: 52.935 MiB, 0.82% gc time)
# 0.940580 seconds (527.77 k allocations: 52.930 MiB, 0.61% gc time)