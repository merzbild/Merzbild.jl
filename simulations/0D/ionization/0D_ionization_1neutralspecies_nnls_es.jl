# Simulate electrons accelerated by a constant electric field collding and ionizing neutrals
# With or without event splitting
# Electrons are merged away using NNLS merging
# The ions are merged away (very coarsely)
# Neutrals are merged using octree merging

include("../../../src/Merzbild.jl")


using ..Merzbild
using Random
using TimerOutputs

function run(seed, E_Tn, n_t,
             n_up_to_total, n_full_up_to_total, threshold_electrons, np_target_electrons_octree;
             adds=0, rate_preserving=true, do_es=true)
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

    if n_full_up_to_total > n_up_to_total
        throw(ErrorException("n_full_up_to_total cannot be larger than n_up_to_total"))
    end

    mim = []
    n_moms = n_full_up_to_total
    for i in 1:n_moms
        append!(mim, compute_multi_index_moments(i))
    end

    for i in n_full_up_to_total+1:n_up_to_total
        append!(mim, [[i, 0, 0], [0, i, 0], [0, 0, i]])
    end

    if length(mim) + 2 >= threshold_electrons
        throw(ErrorException("# of post-merge particles larger than threshold (NNLS)"))
    end
    if np_target_electrons_octree >= threshold_electrons
        throw(ErrorException("# of post-merge particles larger than threshold (octree)"))
    end

    Efield_int = round(Int64, E_Tn)

    fname =  "../../Data/PIC_DSMC/NNLS_merging/0D_ionization/Ar_e/NNLS_macc"
    
    if do_es
        fname = fname * "_es/"
    else
        fname = fname * "/"
    end

    if rate_preserving
        if adds == 0
            fname = fname * "Ar_$(Efield_int)Tn_NNLSrate_$(n_full_up_to_total)full_upto$(n_up_to_total)_$(threshold_electrons).nc"
        else
            fname = fname * "Ar_$(Efield_int)Tn_NNLSrate_$(n_full_up_to_total)full_upto$(n_up_to_total)_$(threshold_electrons)_seed$(adds).nc"
        end
    else
        if adds == 0
            fname = fname * "Ar_$(Efield_int)Tn_NNLS_$(n_full_up_to_total)full_upto$(n_up_to_total)_$(threshold_electrons).nc"
        else
            fname = fname * "Ar_$(Efield_int)Tn_NNLS_$(n_full_up_to_total)full_upto$(n_up_to_total)_$(threshold_electrons)_seed$(adds).nc"
        end
    end

    Random.seed!(seed+adds)
    rng::Xoshiro = Xoshiro(seed+adds)

    reset_timer!()

    species_list::Vector{Species} = load_species_list("data/particles.toml", ["Ar", "Ar+", "e-"])
    interaction_data::Array{Interaction, 2} = load_interaction_data_with_dummy("data/vhs.toml", species_list)


    n_e_interactions = load_electron_neutral_interactions(species_list, "../../Data/cross_sections/Ar_IST_Lisbon.xml",
                                                          Dict("Ar" => "IST-Lisbon"),
                                                          Dict("Ar" => ScatteringIsotropic),
                                                          Dict("Ar" => ElectronEnergySplitEqual))

    ref_cs_elastic = maximum(n_e_interactions.elastic[1].data.sigma)
    ref_cs_ionization = maximum(n_e_interactions.ionization[1].data.sigma)

    println(ref_cs_elastic, ", ", ref_cs_ionization)
    
    @timeit "NNLSinit" mnnls = create_nnls_merging(mim, threshold_electrons)
    @timeit "NNLSinit RP" mnnls_rp = create_nnls_merging_rate_conserving(mim, threshold_electrons)

    println("Total conservation eqns: $(length(mnnls.rhs_vector))")

    n_e_cs = create_computed_crosssections(n_e_interactions)

    index_neutral = 1
    index_ion = 2
    index_electron = 3


    nv_heavy = 20  # init neutrals and ions on a coarser grid
    nv_electrons = 40
    np_base_heavy = nv_heavy^3  # some initial guess on # of particles in simulation
    np_base_electrons = nv_electrons^3  # some initial guess on # of particles in simulation

    oc = create_merging_octree(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel, max_Nbins=6000)
    oc_electrons = create_merging_octree(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel, max_Nbins=6000)
    mg_ions = create_merging_grid(Nmerging_ions, Nmerging_ions, Nmerging_ions, 3.5)

   
    particles::Vector{Vector{Particle}} = [Vector{Particle}(undef, np_base_heavy),
                                           Vector{Particle}(undef, np_base_heavy),
                                           Vector{Particle}(undef, np_base_electrons)]
    n_sampled = [0, 0, 0]

    for (index, (n_dens, T0, nv)) in enumerate(zip([n_dens_neutrals, n_dens_ions, n_dens_e],
                                                   [T0, T0, T0_e],
                                                   [nv_heavy, nv_heavy, nv_electrons]))
        n_sampled[index] = sample_maxwellian_on_grid!(rng, particles[index], nv, species_list[index].mass, T0, n_dens,
                                                      0.0, 1.0, 0.0, 1.0, 0.0, 1.0;
                                                      v_mult=3.5, cutoff_mult=8.0, noise=0.0, v_offset=[0.0, 0.0, 0.0])
    end

    pia = create_particle_indexer_array(n_sampled)

    phys_props::PhysProps = create_props(1, 3, [], Tref=T0)
    compute_props!(particles, pia, species_list, phys_props)

    ds = create_netcdf_phys_props(fname, species_list, phys_props)
    write_netcdf_phys_props(ds, phys_props, 0)

    collision_factors = create_collision_factors(3)
    collision_data = create_collision_data()

    # neutral-neutral
    Fnum_neutral_mean = n_dens_neutrals / n_sampled[1]
    collision_factors[1,1].sigma_g_w_max = estimate_sigma_g_w_max(interaction_data[1,1],
                                                                  species_list[1], T0,
                                                                  Fnum_neutral_mean)
                             
    s1 = index_neutral
    s2 = index_electron
    s3 = index_ion
    # neutral-electron
    estimate_sigma_g_w_max_ntc_n_e!(rng, collision_factors[s1,s2], collision_data, interaction_data,
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
            @timeit "merge n" merge_octree_N2_based!(oc, particles[1], pia, 1, 1, np_target_neutrals)
        end

        if pia.n_total[2] > threshold_ion
            @timeit "merge i" merge_grid_based!(mg_ions, particles[2], pia, 1, 2, species_list, phys_props)
        end

        vref = sqrt(2 * k_B * phys_props.T[1,index_electron] / species_list[index_electron].mass)

        if pia.n_total[3] > threshold_electrons
            if firstm
                if rate_preserving
                    @timeit "NNLSmerge: 1st time" nnls_success_flag = merge_nnls_based_rate_preserving!(rng, mnnls_rp, 
                    interaction_data, n_e_interactions, n_e_cs,
                    particles[3], pia, 1, 3, index_neutral,
                    vref, ref_cs_elastic, ref_cs_ionization)
                else
                    @timeit "NNLSmerge: 1st time" nnls_success_flag = merge_nnls_based!(rng, mnnls, particles[3], pia, 1, 3, vref)
                end
                firstm = false
            else
                if rate_preserving
                    @timeit "NNLSmerge" nnls_success_flag = merge_nnls_based_rate_preserving!(rng, mnnls_rp, 
                    interaction_data, n_e_interactions, n_e_cs,
                    particles[3], pia, 1, 3, index_neutral,
                    vref, ref_cs_elastic, ref_cs_ionization)
                else
                    @timeit "NNLSmerge" nnls_success_flag = merge_nnls_based!(rng, mnnls, particles[3], pia, 1, 3, vref)
                end
                
            end

            if nnls_success_flag == -1
                println("Resorting to octree merging")
                @timeit "Octreemerge e" merge_octree_N2_based!(oc_electrons, particles[3], pia, 1, 3, np_target_electrons_octree)
            end
        end


        @timeit "acc e" accelerate_constant_field_x!(particles[index_electron],
                                     pia, 1, index_electron, species_list,
                                     E_field, Δt)
        
        @timeit "props" compute_props!(particles, pia, species_list, phys_props)
        @timeit "I/O" Merzbild.write_netcdf_phys_props(ds, phys_props, ts, sync_freq=5000)

    end
    print_timer()
    close_netcdf(ds)
end

const paramset = [8, 6, 150, 100]
run(1234, 400.0, 500000, paramset[1], paramset[2], paramset[3], paramset[4], adds=0, rate_preserving=false, do_es=false)
run(1234, 400.0, 500000, paramset[1], paramset[2], paramset[3], paramset[4], adds=0, rate_preserving=false, do_es=true)

run(1234, 400.0, 500000, paramset[1], paramset[2], paramset[3], paramset[4], adds=0, rate_preserving=true, do_es=false)
run(1234, 400.0, 500000, paramset[1], paramset[2], paramset[3], paramset[4], adds=0, rate_preserving=true, do_es=true)
# const params = [[10, 8, 250, 170], [9, 7, 200, 140], [8, 6, 150, 100], [7, 5, 100, 70], [6, 4, 75, 50], [5, 3, 50, 30]]
# const seeds = [1, 2, 3, 4, 5, 6, 7]
# for sadd in seeds
#     for paramset in params
#         run(1234, 400.0, 500000, paramset[1], paramset[2], paramset[3], paramset[4], adds=sadd)
#     end
# end