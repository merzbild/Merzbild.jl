

# include("../../../src/Merzbild.jl")

using Merzbild
using Random
using TimerOutputs

"""
    run(seed, E_Tn, n_t,
             n_full_up_to_total, threshold_electrons, np_target_electrons_octree;
             adds=0, rate_preserving=:off, do_es=true)
             
Simulate electrons accelerated by a constant electric field collding and ionizing neutrals, with or without event splitting.
Electrons are merged away using NNLS merging, the ions are merged away (very coarsely),
neutrals are merged using octree merging.
Requires external data (not included in the repository), available from LXCat!
Data is written to `scratch/data` by default, adjust `fname` if necessary.
**Warning**: this simulation produces very large amounts of data when running over multiple random
seed values and parameter values.

Positional arguments:
* `seed`: random seed value
* `E_Tn`: electric field value in Townsend
* `n_t`: number of timesteps to run for
* `n_full_up_to_total`: all mixed moments up to this order are conserved during merging
* `threshold_electrons`: threshold number of electrons after which they are merged
* `np_target_electrons_octree`: in case NNLS merging fails and backup NNLS merging fails, resorting
to octree merging with this target number of electrons

Keyword arguments:
* `adds`: value added to `seed` to obtain the random seed used in the simulation, if not 0, will be also appended
to end of output filename (use for ensemble simulations)
* `rate_preserving`: if `:off`, standard NNLS merging is used, if `:approximate`, approximate rate-preserving NNLS merging is used,
if `:exact`, exact rate-preserving NNLS merging is used
* `do_es`: if `true`, event splitting is used for electron-neutral collisions
"""
function run(seed, E_Tn, n_t,
             n_full_up_to_total, threshold_electrons, np_target_electrons_octree,
             cross_section_filepath;
             adds=0, rate_preserving=:off, do_es=true)
    
            
    T0 = 300.0
    T0_e = Merzbild.eV * 2.0  # T_e(t=0) = 2eV
    n_dens_neutrals = 1e23
    n_dens_e = 1e-7 * n_dens_neutrals
    n_dens_ions = n_dens_e

    E_field = E_Tn * n_dens_neutrals * 1e-21 # convert Tn to V/m

    Δt = 5e-14
    V = 1.0

    threshold_neutrals = 200
    threshold_ion = 100

    np_target_neutrals = 100
    Nmerging_ions = 1  # 8-16 after merging

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

    if length(mim) + 2 >= threshold_electrons
        throw(ErrorException("# of post-merge particles larger than threshold (NNLS)"))
    end
    if np_target_electrons_octree >= threshold_electrons
        throw(ErrorException("# of post-merge particles larger than threshold (octree)"))
    end

    Efield_int = round(Int64, E_Tn)

    fname =  "scratch/data/"

    addstr = ""

    if do_es
        addstr = "_es"
    end
    
    if rate_preserving == :approximate
        fname = fname * "ionization_Ar_$(Efield_int)Tn_NNLSrate_approx"
    elseif rate_preserving == :exact
        fname = fname * "ionization_Ar_$(Efield_int)Tn_NNLSrate_exact"
    else
        fname = fname * "ionization_Ar_$(Efield_int)Tn_NNLS"
    end

    if adds == 0
        fname = fname * "_$(n_full_up_to_total)full_$(threshold_electrons)$(addstr).nc"
    else
        fname = fname * "_$(n_full_up_to_total)full_$(threshold_electrons)$(addstr)_seed$(adds).nc"
    end

    println(fname)
    Random.seed!(seed+adds)
    rng::Xoshiro = Xoshiro(seed+adds)

    reset_timer!()

    species_data::Vector{Species} = load_species_data("data/particles.toml", ["Ar", "Ar+", "e-"])
    interaction_data::Array{Interaction, 2} = load_interaction_data_with_dummy("data/vhs.toml", species_data)

    n_e_interactions = load_electron_neutral_interactions(species_data, cross_section_filepath,
                                                          Dict("Ar" => "IST-Lisbon"),
                                                          Dict("Ar" => ScatteringIsotropic),
                                                          Dict("Ar" => ElectronEnergySplitEqual))

    ref_cs_elastic = maximum(n_e_interactions.elastic[1].data.sigma)
    ref_cs_ionization = maximum(n_e_interactions.ionization[1].data.sigma)

    println(ref_cs_elastic, ", ", ref_cs_ionization)
    
    @timeit "NNLSinit" mnnls = NNLSMerge(mim, threshold_electrons, rate_preserving=false)
    @timeit "NNLSinit" mnnls_backup = NNLSMerge(mim_backup, threshold_electrons, rate_preserving=(rate_preserving != :off))
    @timeit "NNLSinit RP" mnnls_rp = NNLSMerge(mim, threshold_electrons, rate_preserving=true)

    println("Total conservation eqns: $(length(mnnls.rhs_vector))")

    n_e_cs = create_computed_crosssections(n_e_interactions)

    index_neutral = 1
    index_ion = 2
    index_electron = 3

    nv_heavy = 20  # init neutrals and ions on a coarser grid
    nv_electrons = 40
    np_base_heavy = nv_heavy^3  # some initial guess on # of particles in simulation
    np_base_electrons = nv_electrons^3  # some initial guess on # of particles in simulation

    oc = OctreeN2Merge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel, max_Nbins=6000)
    oc_electrons = OctreeN2Merge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel, max_Nbins=6000)
    mg_ions = GridN2Merge(Nmerging_ions, Nmerging_ions, Nmerging_ions, 3.5)

    particles = [ParticleVector(np_base_heavy),
                 ParticleVector(np_base_heavy),
                 ParticleVector(np_base_electrons)]
    n_sampled = [0, 0, 0]

    for (index, (n_dens, T0, nv)) in enumerate(zip([n_dens_neutrals, n_dens_ions, n_dens_e],
                                                   [T0, T0, T0_e],
                                                   [nv_heavy, nv_heavy, nv_electrons]))
        n_sampled[index] = sample_maxwellian_on_grid!(rng, particles[index], nv, species_data[index].mass, T0, n_dens,
                                                      0.0, 1.0, 0.0, 1.0, 0.0, 1.0;
                                                      v_mult=3.5, cutoff_mult=8.0, noise=0.0, v_offset=[0.0, 0.0, 0.0])
    end

    pia = ParticleIndexerArray(n_sampled)

    phys_props::PhysProps = PhysProps(1, 3, [], Tref=T0)
    compute_props!(particles, pia, species_data, phys_props)

    ds = NCDataHolder(fname, species_data, phys_props)
    write_netcdf(ds, phys_props, 0)

    collision_factors = create_collision_factors_array(3)
    collision_data = CollisionData()

    vref = sqrt(2 * k_B * 3 * 11605.0 / species_data[1].mass)

    if pia.n_total[1] > threshold_neutrals
        @timeit "merge n t=0" merge_octree_N2_based!(rng, oc, particles[1], pia, 1, 1, np_target_neutrals)
    end

    if pia.n_total[2] > threshold_ion
        @timeit "merge i t=0" merge_grid_based!(rng, mg_ions, particles[2], pia, 1, 2, species_data, phys_props)
    end

    if pia.n_total[3] > threshold_electrons
        @timeit "merge e t=0" nnls_success_flag = merge_nnls_based!(rng, mnnls, particles[3], pia, 1, 3;
                                                                    vref=vref, scaling=:variance, centered_at_mean=false, v_multipliers=[], iteration_mult=5)
        
        if nnls_success_flag == -1
            println("resorting to octree for electrons at t=0")
            merge_octree_N2_based!(rng, oc_electrons, particles[3], pia, 1, 3, np_target_electrons_octree)
        end
    end

    squash_pia!(particles, pia)

    # neutral-neutral
    Fnum_neutral_mean = n_dens_neutrals / pia.indexer[1,1].n_local
    collision_factors[1,1].sigma_g_w_max = estimate_sigma_g_w_max(interaction_data[1,1],
                                                                  species_data[1], T0,
                                                                  Fnum_neutral_mean)
                             
    s1 = index_neutral
    s2 = index_electron
    s3 = index_ion
    # neutral-electron
    estimate_sigma_g_w_max_ntc_n_e!(rng, collision_factors[s1,s2,1], collision_data, interaction_data,
                                    n_e_interactions, n_e_cs,
                                    particles[s1], particles[s2], pia, 1, s1, s2, Δt, V, min_coll=15, n_loops=6)

    firstm = true
    nnls_success_flag = 1

    for ts in 1:n_t

        if ts % 20000 == 0
            println(ts)
        end
        # collide neutrals and neutrals
        @timeit "coll n-n" ntc!(rng, collision_factors[s1,s1,1], collision_data, interaction_data, particles[s1], pia, 1, s1, Δt, V)

        if do_es
            @timeit "coll n-e ES" ntc_n_e_es!(rng, collision_factors[s1,s2,1], collision_data, interaction_data,
                                            n_e_interactions, n_e_cs,
                                            particles[s1], particles[s2], particles[s3], pia, 1, s1, s2, s3, Δt, V)
        else
            @timeit "coll n-e" ntc_n_e!(rng, collision_factors[s1,s2,1], collision_data, interaction_data,
                                        n_e_interactions, n_e_cs,
                                        particles[s1], particles[s2], particles[s3], pia, 1, s1, s2, s3, Δt, V)
        end

        if pia.n_total[1] > threshold_neutrals
            @timeit "merge n" merge_octree_N2_based!(rng, oc, particles[1], pia, 1, 1, np_target_neutrals)
        end

        if pia.n_total[2] > threshold_ion
            @timeit "merge i" merge_grid_based!(rng, mg_ions, particles[2], pia, 1, 2, species_data, phys_props)
        end

        vref = sqrt(2 * k_B * phys_props.T[1,index_electron] / species_data[index_electron].mass)

        if pia.n_total[3] > threshold_electrons
            # println("wat")
            if rate_preserving == :approximate
                @timeit "NNLSmergeARP e" nnls_success_flag = merge_nnls_based_rate_preserving!(rng, mnnls_rp, 
                interaction_data, n_e_interactions, n_e_cs,
                particles[3], pia, 1, 3, index_neutral,
                ref_cs_elastic, ref_cs_ionization; centered_at_mean=false, v_multipliers=[], vref=vref, scaling=:variance, iteration_mult=5)
            elseif rate_preserving == :exact
                @timeit "NNLSmergeERP e" nnls_success_flag = merge_nnls_based_rate_preserving!(rng, mnnls_rp, 
                interaction_data, n_e_interactions, n_e_cs,
                particles[3], particles[1], pia, 1, 3, index_neutral,
                ref_cs_elastic, ref_cs_ionization; vref=vref, scaling=:variance, iteration_mult=5)
            else
                @timeit "NNLSmerge e" nnls_success_flag = merge_nnls_based!(rng, mnnls, particles[3], pia, 1, 3;
                                                                          centered_at_mean=false, v_multipliers=[], vref=vref, scaling=:variance, iteration_mult=5)
            end
            
            if nnls_success_flag == -1

                
                @timeit "NNLSmerge bup e" nnls_success_flag = merge_nnls_based!(rng, mnnls_backup, particles[3], pia, 1, 3;
                                                      centered_at_mean=false, v_multipliers=[], vref=vref, scaling=:variance, iteration_mult=5) 

                if nnls_success_flag == -1
                    println("Resorting to octree merging")
                    @timeit "Octreemerge e" merge_octree_N2_based!(rng, oc_electrons, particles[3], pia, 1, 3, np_target_electrons_octree)
                end
            end
        end

        @timeit "acc e" accelerate_constant_field_x!(particles[index_electron],
                                     pia, 1, index_electron, species_data,
                                     E_field, Δt)
        
        @timeit "props" compute_props!(particles, pia, species_data, phys_props)
        @timeit "I/O" Merzbild.write_netcdf(ds, phys_props, ts, sync_freq=5000)

    end
    print_timer()
    close_netcdf(ds)
end

# paramset is a list with 3 elements: number of velocity moments conserved, threshold number of particles,
# target number of particles for backup octree merging (in case octree fails)
# if NNLS with param[1] velocity moments fails, a backup NNLS is called with param[1]-1 moments, if that fails - octree is called
paramset = [5, 69, 58]
external_E_field_Tn = 400.0

# path to IST-Lisbon cross-section data
ist_lisbon_filepath = "../../Data/cross_sections/Ar_IST_Lisbon.xml"

# set this to a smaller value (i.e. 10000) if you just want to check that the file runs
n_t = 500000
for do_event_splitting in [false, true]  # try out different collision schemes
    run(1234, external_E_field_Tn, n_t, paramset[1], paramset[2], paramset[3],
        ist_lisbon_filepath; adds=0, rate_preserving=:off, do_es=do_event_splitting)
    run(1234, external_E_field_Tn, n_t, paramset[1], paramset[2], paramset[3],
        ist_lisbon_filepath; adds=0, rate_preserving=:approximate, do_es=do_event_splitting)
    run(1234, external_E_field_Tn, n_t, paramset[1], paramset[2], paramset[3],
        ist_lisbon_filepath; adds=0, rate_preserving=:exact, do_es=do_event_splitting)
end

#  # Uncomment set-up below to run over the parameter sets used for "Moment-preserving particle merging via non-negative least squares"
#  # the 4th value in each parameter list is the number of ensembles that are run with different random seeds 
#  # these simulations produce large amounts of data - a single 400 Tn simulation file is 200+ MB in size,
#  # a single 100 Tn simulation file is 2+ GB in size

# params = [[4, 41, 38, 64], [5, 62, 58, 64], [6, 95, 88, 64], [7, 131, 122, 16], [8, 178, 166, 16], [9, 236, 220, 16]]  # 1.075

#  # we need more timesteps for the weaker field
# for (n_t, external_E_field_Tn) in zip([500000, 5000000], [400.0, 100.0])

#     # iterate over parameters 
#     for paramset in params
#         for sadd in 0:paramset[4]-1
#             run(1234, external_E_field_Tn, n_t,
#                 paramset[1], paramset[2], paramset[3], adds=sadd, rate_preserving=:off, do_es=true)
#         end
#     end

#     for paramset in params
#         for sadd in 0:paramset[4]-1
#             run(1234, external_E_field_Tn, n_t,
#                 paramset[1], paramset[2], paramset[3], adds=sadd, rate_preserving=:approximate, do_es=true)
#         end
#     end

#     for paramset in params
#         for sadd in 0:paramset[4]-1
#             run(1234, external_E_field_Tn, n_t,
#                 paramset[1], paramset[2], paramset[3], adds=sadd, rate_preserving=:exact, do_es=true)
#         end
#     end
# end
