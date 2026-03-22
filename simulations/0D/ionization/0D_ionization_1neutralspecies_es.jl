# Simulate electrons accelerated by a constant electric field collding and ionizing neutrals
# With or without event splitting
# The ions are merged away (very coarsely)
# Electrons and neutrals are merged using octree merging
# Requires external data (not included in the repository), available from LXCat!

include("../../../src/Merzbild.jl")


using ..Merzbild
using Random
using TimerOutputs

"""
    run(seed, E_Tn, n_t, threshold_electrons, np_target_electrons;
        merging_bin_split=OctreeBinMidSplit, adds=0, do_es=true)
             
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
* `threshold_electrons`: threshold number of electrons after which they are merged
* `np_target_electrons`: in case NNLS merging fails and backup NNLS merging fails, resorting
to octree merging with this target number of electrons

Keyword arguments:
* `merging_bin_split`: how merging bins are split (`OctreeBinMidSplit` is the only recommended option)
* `adds`: value added to `seed` to obtain the random seed used in the simulation, if not 0, will be also appended
to end of output filename (use for ensemble simulations)
* `do_es`: if `true`, event splitting is used for electron-neutral collisions
"""
function run(seed, E_Tn, n_t, threshold_electrons, np_target_electrons,
             cross_section_filepath; merging_bin_split=OctreeBinMidSplit, adds=0, do_es=true)
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

    Efield_int = round(Int64, E_Tn)

    fname = "scratch/data/"
    addstr = ""

    if do_es
        addstr = "_es"
    end

    fname = fname * "ionization_Ar_$(Efield_int)Tn_octree"
    if merging_bin_split == OctreeBinMidSplit
        fname = fname * "_mid"
    elseif merging_bin_split == OctreeBinMeanSplit
        fname = "_mean"
    elseif merging_bin_split == OctreeBinMedianSplit
        fname = "_median"
    end

    fname = fname * "_$(threshold_electrons)_to_$(np_target_electrons)"

    if adds == 0
        fname = fname * "$(addstr).nc"
    else
        fname = fname * "$(addstr)_seed$(adds).nc"
    end

    Random.seed!(seed+adds)
    rng::Xoshiro = Xoshiro(seed+adds)

    reset_timer!()

    species_data::Vector{Species} = load_species_data("data/particles.toml", ["Ar", "Ar+", "e-"])
    interaction_data::Array{Interaction, 2} = load_interaction_data_with_dummy("data/vhs.toml", species_data)


    n_e_interactions = load_electron_neutral_interactions(species_data, cross_section_filepath,
                                                          Dict("Ar" => "IST-Lisbon"),
                                                          Dict("Ar" => ScatteringIsotropic),
                                                          Dict("Ar" => ElectronEnergySplitEqual))

    n_e_cs = create_computed_crosssections(n_e_interactions)

    index_neutral = 1
    index_ion = 2
    index_electron = 3

    nv_heavy = 20  # init neutrals and ions on a coarser grid
    nv_electrons = 40
    np_base_heavy = nv_heavy^3  # some initial guess on # of particles in simulation
    np_base_electrons = nv_electrons^3  # some initial guess on # of particles in simulation

    oc = OctreeN2Merge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel, max_Nbins=8000)
    oc_electrons = OctreeN2Merge(merging_bin_split; init_bin_bounds=OctreeInitBinMinMaxVel, max_Nbins=8000)
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

    if pia.n_total[1] > threshold_neutrals
        @timeit "merge n" merge_octree_N2_based!(rng, oc, particles[1], pia, 1, 1, np_target_neutrals)
    end

    if pia.n_total[2] > threshold_ion
        @timeit "merge i" merge_grid_based!(rng, mg_ions, particles[2], pia, 1, 2, species_data, phys_props)
    end

    if pia.n_total[3] > threshold_electrons
        @timeit "merge e" merge_octree_N2_based!(rng, oc_electrons, particles[3], pia, 1, 3, np_target_electrons)
    end

    squash_pia!(particles, pia)

    ds = NCDataHolder(fname, ["v", "moments"], species_data, phys_props)
    write_netcdf(ds, phys_props, 0)

    collision_factors = create_collision_factors_array(3)
    collision_data = CollisionData()

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

        if pia.n_total[3] > threshold_electrons
            @timeit "merge e" merge_octree_N2_based!(rng, oc_electrons, particles[3], pia, 1, 3, np_target_electrons)
        end

        @timeit "acc e" accelerate_constant_field_x!(particles[index_electron],
                                     pia, 1, index_electron, species_data,
                                     E_field, Δt)
        
        @timeit "props" compute_props!(particles, pia, species_data, phys_props)
        @timeit "I/O" write_netcdf(ds, phys_props, ts, sync_freq=1000)

    end
    print_timer()
    close_netcdf(ds)
end

# paramset is a list with 2 elements: threshold number of particles,
# target number of particles for merging
paramset = [45, 38]
external_E_field_Tn = 400.0
n_t = 500000

# path to IST-Lisbon cross-section data
ist_lisbon_filepath = "../../Data/cross_sections/Ar_IST_Lisbon.xml"

for do_event_splitting in [false, true]  # try out different collision schemes
    run(1234, external_E_field_Tn, n_t, paramset[1], paramset[2], ist_lisbon_filepath; merging_bin_split=OctreeBinMidSplit, adds=0, do_es=do_event_splitting)
end


#  # Uncomment set-up below to run over the parameter sets used for "Moment-preserving particle merging via non-negative least squares"
#  # the 3rd value in each parameter list is the number of ensembles that are run with different random seeds 

# params = [[41, 38, 64], [62, 58, 64], [95, 88, 64], [131, 122, 16], [178, 166, 16], [236, 220, 16]] 

#  # we need more timesteps for the weaker field
# for (n_t, external_E_field_Tn) in zip([500000, 5000000], [400.0, 100.0])

#     # iterate over parameters 
#     for paramset in params
#         for sadd in 0:params[3]-1
#             run(1234, external_E_field_Tn, n_t,
#                 paramset[1], paramset[2], ist_lisbon_filepath; merging_bin_split=OctreeBinMidSplit, adds=sadd, do_es=true)
#         end
#     end

#     for paramset in params
#         for sadd in 0:params[3]-1
#             run(1234, external_E_field_Tn, n_t,
#                 paramset[1], paramset[2], ist_lisbon_filepath; merging_bin_split=OctreeBinMidSplit, adds=sadd, do_es=true)
#         end
#     end

#     for paramset in params
#         for sadd in 0:params[3]-1
#             run(1234, external_E_field_Tn, n_t,
#                 paramset[1], paramset[2], ist_lisbon_filepath; merging_bin_split=OctreeBinMidSplit, adds=sadd, do_es=true)
#         end
#     end
# end
