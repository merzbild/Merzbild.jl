@testset "ionization cross-section indexing" begin
    # test that internal indexing of species works correctly

    # first we test with small T0_e
    E_Tn = 100 # field strength in Tn

    T0 = 300.0
    T0_e = Merzbild.eV * 5.0
    n_dens_neutrals = 1e23
    n_dens_e = 1e-5 * n_dens_neutrals
    n_dens_ions = n_dens_e

    E_field = E_Tn * n_dens_neutrals * 1e-21 # convert Tn to V/m

    Δt = 5e-14
    V = 1.0

    merging_bin_split = OctreeBinMidSplit

    n_t = 100

    threshold_electrons = 150
    np_target_electrons = 120

    threshold_neutrals = 100
    threshold_ion = 50

    np_target_neutrals = 100
    Nmerging_ions = 1  # 8-16 after merging

    seed = 123
    Random.seed!(seed)
    rng = StableRNG(seed)


    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    interaction_data_path = joinpath(@__DIR__, "..", "data", "vhs.toml")
    cs_data_path = joinpath(@__DIR__, "..", "data", "test_neutral_electron_data.xml")

    species_data::Vector{Species} = load_species_data(particles_data_path, ["Ar", "Ar+", "Xe+", "Xe", "e-"])
    interaction_data::Array{Interaction, 2} = load_interaction_data_with_dummy(interaction_data_path, species_data)

    n_e_interactions = load_electron_neutral_interactions(species_data, cs_data_path,
                                                          Dict("Ar" => "ConstantDB",
                                                               "Xe" => "ZeroDB"),
                                                          Dict("Ar" => ScatteringIsotropic,
                                                               "Xe" => ScatteringIsotropic),
                                                          Dict("Ar" => ElectronEnergySplitEqual,
                                                               "Xe" => ElectronEnergySplitEqual))

    n_e_cs = create_computed_crosssections(n_e_interactions)

    @test n_e_interactions.n_neutrals == 2
    @test length(n_e_interactions.neutral_indexer) == 5
    @test n_e_interactions.neutral_indexer[1] == 1
    @test n_e_interactions.neutral_indexer[4] == 2

    for i in [2,3,5]
        @test n_e_interactions.neutral_indexer[i] == -1
    end

    nv_heavy = 15  # init neutrals and ions on a coarser grid
    nv_electrons = 20
    np_base_heavy = nv_heavy^3  # some initial guess on # of particles in simulation
    np_base_electrons = nv_electrons^3  # some initial guess on # of particles in simulation

    oc = OctreeN2Merge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel, max_Nbins=6000)
    oc_electrons = OctreeN2Merge(merging_bin_split; init_bin_bounds=OctreeInitBinMinMaxVel, max_Nbins=6000)
    mg_ions = GridN2Merge(Nmerging_ions, Nmerging_ions, Nmerging_ions, 3.5)

    particles = [ParticleVector(np_base_heavy),
                 ParticleVector(np_base_heavy),
                 ParticleVector(np_base_heavy),
                 ParticleVector(np_base_heavy),
                 ParticleVector(np_base_electrons)]
    n_sampled = [0, 0, 0, 0, 0]

    for (index, (n_dens, T0, nv)) in enumerate(zip([n_dens_neutrals, n_dens_ions, n_dens_ions, n_dens_neutrals, n_dens_e],
                                                   [T0, T0, T0, T0, T0_e],
                                                   [nv_heavy, nv_heavy, nv_heavy, nv_heavy, nv_electrons]))
        n_sampled[index] = sample_maxwellian_on_grid!(rng, particles[index], nv, species_data[index].mass, T0, n_dens,
                                                      0.0, 1.0, 0.0, 1.0, 0.0, 1.0;
                                                      v_mult=3.5, cutoff_mult=8.0, noise=0.0, v_offset=[0.0, 0.0, 0.0])
    end

    pia = ParticleIndexerArray(n_sampled)

    phys_props::PhysProps = PhysProps(1, 5, [], Tref=T0)

    collision_factors = create_collision_factors_array(5)
    collision_data = CollisionData()
    compute_props!(particles, pia, species_data, phys_props)

    if pia.n_total[1] > threshold_neutrals
        merge_octree_N2_based!(rng, oc, particles[1], pia, 1, 1, np_target_neutrals)
    end

    if pia.n_total[2] > threshold_ion
        merge_grid_based!(rng, mg_ions, particles[2], pia, 1, 2, species_data, phys_props)
    end

    if pia.n_total[3] > threshold_ion
        merge_grid_based!(rng, mg_ions, particles[3], pia, 1, 3, species_data, phys_props)
    end

    if pia.n_total[4] > threshold_ion
        merge_grid_based!(rng, mg_ions, particles[4], pia, 1, 4, species_data, phys_props)
    end

    if pia.n_total[5] > threshold_electrons
        merge_octree_N2_based!(rng, oc_electrons, particles[5], pia, 1, 5, np_target_electrons)
    end

    # neutral-neutral
    Fnum_neutral_mean = n_dens_neutrals / np_target_neutrals
    collision_factors[1,1,1].sigma_g_w_max = estimate_sigma_g_w_max(interaction_data[1,1],
                                                                    species_data[1], T0,
                                                                    Fnum_neutral_mean)
    
    collision_factors[4,4,1].sigma_g_w_max = estimate_sigma_g_w_max(interaction_data[4,4],
                                                                    species_data[4], T0,
                                                                    Fnum_neutral_mean)
    
    for s1 in [1,4]
        estimate_sigma_g_w_max_ntc_n_e!(rng, collision_factors[s1,5,1], collision_data, interaction_data,
                                        n_e_interactions, n_e_cs,
                                        particles[s1], particles[5], pia, 1, s1, 5, Δt, V, min_coll=2, n_loops=1)

        cf_e_n = collision_factors[s1,5,1].sigma_g_w_max
        @test cf_e_n > 0

        # test that more samples lead to larger estimate
        estimate_sigma_g_w_max_ntc_n_e!(rng, collision_factors[s1,5,1], collision_data, interaction_data,
                                        n_e_interactions, n_e_cs,
                                        particles[s1], particles[5], pia, 1, s1, 5, Δt, V, min_coll=15, n_loops=5)
        @test collision_factors[s1,5,1].sigma_g_w_max > cf_e_n
    end
end