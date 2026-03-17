# BKW with NNLS merging + timing outputs
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
using TimerOutputs

"""
    run(seed, n_full_up_to_total, threshold, ntarget_octree)

Run the BKW simulation with NNLS merging.

Positional arguments:
* `seed`: random seed value
* `n_full_up_to_total`: all mixed moments up to this order are conserved during merging
* `threshold`: threshold number of particles after which they are merged
* `ntarget_octree`: in case NNLS merging fails and backup NNLS merging fails, resorting
to octree merging with this target number of particles
"""
function run(seed, n_full_up_to_total, threshold, ntarget_octree)
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    # reset_timer!()
    species_data::Vector{Species} = load_species_data("data/particles.toml", "Ar")
    interaction_data::Array{Interaction, 2} = load_interaction_data("data/pseudo_maxwell.toml", species_data)

    dt_scaled = 0.025
    n_t = 600

    nv = 36
    np_base = 40^3  # some initial guess on # of particle in simulation

    mim = []
    n_moms = n_full_up_to_total
    for i in 1:n_moms
        append!(mim, compute_multi_index_moments(i))
    end

    mim_backup = []
    n_moms_backup = n_full_up_to_total - 1
    for i in 1:n_moms_backup
        append!(mim_backup, compute_multi_index_moments(i))
    end

    mnnls = NNLSMerge(mim, threshold)
    mnnls_backup = NNLSMerge(mim_backup, threshold)
    ocm = OctreeN2Merge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel, max_Nbins=6000)

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
                                v_mult=4.0, cutoff_mult=4.0, noise=0.0, v_offset=[0.0, 0.0, 0.0])

    pia = ParticleIndexerArray(n_sampled)

    phys_props::PhysProps = PhysProps(1, 1, moments_list, Tref=T0)
    compute_props_with_total_moments!(particles, pia, species_data, phys_props)

    ds = NCDataHolder("scratch/data/bkw_nnls_$(n_full_up_to_total)full_$(threshold)_$(seed).nc", species_data, phys_props)

    write_netcdf(ds, phys_props, 0)

    if phys_props.np[1,1] > threshold
        nnls_success_flag = merge_nnls_based!(rng, mnnls, particles[1], pia, 1, 1;
                                                                            vref=vref, scaling=:variance, centered_at_mean=false, v_multipliers=[], iteration_mult=4)

        if nnls_success_flag == -1
            nnls_success_flag = merge_nnls_based!(rng, mnnls_backup, particles[1], pia, 1, 1;
                                                                    vref=vref, scaling=:variance, centered_at_mean=false, v_multipliers=[], iteration_mult=4)
            
            if nnls_success_flag == -1
                merge_octree_N2_based!(rng, ocm, particles[1], pia, 1, 1, ntarget_octree)
            end
        end
    end

    collision_factors::CollisionFactors = CollisionFactors()
    collision_data::CollisionData = CollisionData()

    # average Fnum is: 
    # Fnum = n_dens/n_sampled
    Fnum = n_dens/pia.indexer[1,1].n_local
    collision_factors.sigma_g_w_max = estimate_sigma_g_w_max(interaction_data[1,1], species_data[1], T0, Fnum)

    Δt::Float64 = dt_scaled * tref
    V::Float64 = 1.0

    nnls_success_flag = 1

    # npp = 0.0

    for ts in 1:n_t
        # if ts % 100 == 0
        #     println(ts)
        # end

        ntc!(rng, collision_factors, collision_data, interaction_data, particles[1], pia, 1, 1, Δt, V)

        # npp += pia.indexer[1,1].n_local
        if pia.indexer[1,1].n_local > threshold
            nnls_success_flag = merge_nnls_based!(rng, mnnls, particles[1], pia, 1, 1;
                                                                      vref=vref, scaling=:variance, centered_at_mean=false, v_multipliers=[], iteration_mult=4)
            if nnls_success_flag == -1
                nnls_success_flag = merge_nnls_based!(rng, mnnls_backup, particles[1], pia, 1, 1;
                                                                      vref=vref, scaling=:variance, centered_at_mean=false, v_multipliers=[], iteration_mult=4)
                
                if nnls_success_flag == -1
                    merge_octree_N2_based!(rng, ocm, particles[1], pia, 1, 1, ntarget_octree)
                end
            end
        end

        squash_pia!(particles, pia, 1)
        
        compute_props_with_total_moments!(particles, pia, species_data, phys_props)
        write_netcdf(ds, phys_props, ts)
    end
    close_netcdf(ds)
end

# run with merging conserving all mixed moments up to order 6; merging is called when particle count exceeds 100
# octree merging with 85 target particles is used when all else fails
run(0, 6, 100, 85)


#  # Uncomment set-up below to run over the parameter sets used for "Moment-preserving particle merging via non-negative least squares"
#  # each element of params is a 4-element vector: first element is the n_full_up_to_total parameter (maximum order of moments conserved)
#  # second is the threshold number of particles
#  # third is the target number of particles for the backup octree merging
#  # fourth is the number of ensembles that are run with different random seeds 

# const params = [[4, 50, 36, 10800], [5, 66, 55, 4800], [6, 100, 85, 2160], [7, 142, 120, 1040], [8, 200, 164, 560], [9, 264, 220, 400]]

# for paramset in params
#     println("$paramset with $(paramset[4]) seeds")
#     for seed in 1:paramset[4]
#         if seed % 100 == 0
#             println(seed, ", ", paramset)
#         end
#         run(seed, paramset[1], paramset[2], paramset[3])
#     end
# end
