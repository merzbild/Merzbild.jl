@testset "full ionization simulation (no ES)" begin

    # Simulate electrons accelerated by a constant electric field colliding and ionizing neutrals
    # Without event splitting
    # The ions are merged away (very coarsely)
    # Electrons and neutrals are merged using octree merging
    # Uses a very simple constant cross-section model
    # running this for 300k timesteps and a T0_e of Merzbild.eV * 2.0
    # with sampling on a finer grid (nv_heavy = 20, nv_electrons = 40)
    # and computing rate based
    # on densities for 100_000 < nt < 290_000
    # gives a rate of 4.0565095021838734e-16
    # reference bolsig value (bolsig version from 07/2024) using same data
    # (input file bolsig_input_const_Ar.dat and cross-section file bolsig_const_Ar_cs.txt)
    # gives a rate of 3.970707E-16

    # first we test with small T0_e
    E_Tn = 100 # field strength in Tn

    T0 = 300.0
    T0_e = Merzbild.eV * 2.0  # T_e(t=0) = 2eV
    n_dens_neutrals = 1e23
    n_dens_e = 1e-7 * n_dens_neutrals
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

    species_data::Vector{Species} = load_species_data(particles_data_path, ["Ar", "Ar+", "e-"])
    interaction_data::Array{Interaction, 2} = load_interaction_data_with_dummy(interaction_data_path, species_data)

    n_e_interactions = load_electron_neutral_interactions(species_data, cs_data_path,
                                                          Dict("Ar" => "ConstantDB"),
                                                          Dict("Ar" => ScatteringIsotropic),
                                                          Dict("Ar" => ElectronEnergySplitEqual))

    n_e_cs = create_computed_crosssections(n_e_interactions)

    index_neutral = 1
    index_ion = 2
    index_electron = 3


    nv_heavy = 15  # init neutrals and ions on a coarser grid
    nv_electrons = 20
    np_base_heavy = nv_heavy^3  # some initial guess on # of particles in simulation
    np_base_electrons = nv_electrons^3  # some initial guess on # of particles in simulation

    oc = OctreeN2Merge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel, max_Nbins=6000)
    oc_electrons = OctreeN2Merge(merging_bin_split; init_bin_bounds=OctreeInitBinMinMaxVel, max_Nbins=6000)
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

    sol_path = joinpath(@__DIR__, "data", "tmp_ion_no_es.nc")

    # exclude v and moments
    ds = NCDataHolder(sol_path, ["v", "moments"], species_data, phys_props)

    compute_props!(particles, pia, species_data, phys_props)
    write_netcdf(ds, phys_props, 0)

    collision_factors = create_collision_factors_array(3)
    collision_data = CollisionData()

    if pia.n_total[1] > threshold_neutrals
        merge_octree_N2_based!(rng, oc, particles[1], pia, 1, 1, np_target_neutrals)
    end

    if pia.n_total[2] > threshold_ion
        merge_grid_based!(rng, mg_ions, particles[2], pia, 1, 2, species_data, phys_props)
    end

    if pia.n_total[3] > threshold_electrons
        merge_octree_N2_based!(rng, oc_electrons, particles[3], pia, 1, 3, np_target_electrons)
    end

    # neutral-neutral
    Fnum_neutral_mean = n_dens_neutrals / np_target_neutrals
    collision_factors[1,1,1].sigma_g_w_max = estimate_sigma_g_w_max(interaction_data[1,1],
                                                                  species_data[1], T0,
                                                                  Fnum_neutral_mean)
                             
    s1 = index_neutral
    s2 = index_electron
    s3 = index_ion
    # neutral-electron

    @test collision_factors[s1,s2,1].sigma_g_w_max == 0
    
    estimate_sigma_g_w_max_ntc_n_e!(rng, collision_factors[s1,s2,1], collision_data, interaction_data,
                                    n_e_interactions, n_e_cs,
                                    particles[s1], particles[s2], pia, 1, s1, s2, Δt, V, min_coll=2, n_loops=1)

    cf_e_n = collision_factors[s1,s2,1].sigma_g_w_max
    @test cf_e_n > 0

    # test that more samples lead to larger estimate
    estimate_sigma_g_w_max_ntc_n_e!(rng, collision_factors[s1,s2,1], collision_data, interaction_data,
                                    n_e_interactions, n_e_cs,
                                    particles[s1], particles[s2], pia, 1, s1, s2, Δt, V, min_coll=15, n_loops=5)
    @test collision_factors[s1,s2,1].sigma_g_w_max > cf_e_n

    for ts in 1:n_t
        # collide neutrals and neutrals
        ntc!(rng, collision_factors[s1,s1,1], collision_data, interaction_data, particles[s1], pia, 1, s1, Δt, V)

        ntc_n_e!(rng, collision_factors[s1,s2,1], collision_data, interaction_data,
                 n_e_interactions, n_e_cs,
                 particles[s1], particles[s2], particles[s3], pia, 1, s1, s2, s3, Δt, V)


        if pia.n_total[1] > threshold_neutrals
            merge_octree_N2_based!(rng, oc, particles[1], pia, 1, 1, np_target_neutrals)
        end

        if pia.n_total[2] > threshold_ion
            merge_grid_based!(rng, mg_ions, particles[2], pia, 1, 2, species_data, phys_props)
        end

        if pia.n_total[3] > threshold_electrons
            merge_octree_N2_based!(rng, oc_electrons, particles[3], pia, 1, 3, np_target_electrons)
        end

        accelerate_constant_field_x!(particles[index_electron],
                                     pia, 1, index_electron, species_data,
                                     E_field, Δt)
        
        compute_props!(particles, pia, species_data, phys_props)
        write_netcdf(ds, phys_props, ts, sync_freq=1000)
    end
    close_netcdf(ds)

    sol = NCDataset(sol_path, "r")
 
    ndens_tot = n_dens_e + n_dens_ions + n_dens_neutrals

    # test overall number density conservation, 1.5e-14 is relative error due to initialization
    ndens_conservation = true

    for t in 1:n_t
        if abs(sum(sol["ndens"][1, :, t]) - ndens_tot) / ndens_tot > 1.5e-14
            ndens_conservation = false
        end
    end

    @test ndens_conservation == true

    # test charge neutrality
    charge_neutrality = true

    for t in 1:n_t
        if abs(sol["ndens"][1, 2, t] - sol["ndens"][1, 3, t]) / sol["ndens"][1, 3, t] > 1.5e-14
            charge_neutrality = false
        end
    end

    @test charge_neutrality == true

    # test compared to ref solution
    ref_sol_path = joinpath(@__DIR__, "data", "ionization_Ar_no_es.nc")
    ref_sol = NCDataset(ref_sol_path, "r")

    @test maximum(abs.(ref_sol["ndens"][1, :, 1:n_t] .- sol["ndens"][1, :, 1:n_t])) < 4 * eps()
    @test maximum(abs.(ref_sol["T"][1, :, 1:n_t] .- sol["T"][1, :, 1:n_t])) < 7.5e-12  # not zero
    # because in v0.7.9 some octree merging computations replaced 1/w with 1 * inv_w

    close(sol)
    rm(sol_path)

    # now we test with higher T0_e and E_tn and should get ionization
    E_Tn = 400
    T0_e = Merzbild.eV * 4.5  # T_e(t=0) = 4.5eV

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


    sol_path = joinpath(@__DIR__, "data", "tmp_ion_no_es_with_ion.nc")
    # exclude v and moments
    ds = NCDataHolder(sol_path, ["v", "moments"], species_data, phys_props)

    compute_props!(particles, pia, species_data, phys_props)
    write_netcdf(ds, phys_props, 0)

    collision_factors = create_collision_factors_array(3)
    collision_data = CollisionData()

    if pia.n_total[1] > threshold_neutrals
        merge_octree_N2_based!(rng, oc, particles[1], pia, 1, 1, np_target_neutrals)
    end

    if pia.n_total[2] > threshold_ion
        merge_grid_based!(rng, mg_ions, particles[2], pia, 1, 2, species_data, phys_props)
    end

    if pia.n_total[3] > threshold_electrons
        merge_octree_N2_based!(rng, oc_electrons, particles[3], pia, 1, 3, np_target_electrons)
    end

    # neutral-neutral
    Fnum_neutral_mean = n_dens_neutrals / np_target_neutrals
    collision_factors[1,1,1].sigma_g_w_max = estimate_sigma_g_w_max(interaction_data[1,1],
                                                                  species_data[1], T0,
                                                                  Fnum_neutral_mean)
                             
    s1 = index_neutral
    s2 = index_electron
    s3 = index_ion
    # neutral-electron

    @test collision_factors[s1,s2,1].sigma_g_w_max == 0
    
    estimate_sigma_g_w_max_ntc_n_e!(rng, collision_factors[s1,s2,1], collision_data, interaction_data,
                                    n_e_interactions, n_e_cs,
                                    particles[s1], particles[s2], pia, 1, s1, s2, Δt, V, min_coll=2, n_loops=1)

    cf_e_n = collision_factors[s1,s2,1].sigma_g_w_max
    @test cf_e_n > 0

    # test that more samples lead to larger estimate
    estimate_sigma_g_w_max_ntc_n_e!(rng, collision_factors[s1,s2,1], collision_data, interaction_data,
                                    n_e_interactions, n_e_cs,
                                    particles[s1], particles[s2], pia, 1, s1, s2, Δt, V, min_coll=15, n_loops=5)
    @test collision_factors[s1,s2,1].sigma_g_w_max > cf_e_n

    for ts in 1:n_t
        # collide neutrals and neutrals
        ntc!(rng, collision_factors[s1,s1,1], collision_data, interaction_data, particles[s1], pia, 1, s1, Δt, V)

        ntc_n_e!(rng, collision_factors[s1,s2,1], collision_data, interaction_data,
                 n_e_interactions, n_e_cs,
                 particles[s1], particles[s2], particles[s3], pia, 1, s1, s2, s3, Δt, V)


        if pia.n_total[1] > threshold_neutrals
            merge_octree_N2_based!(rng, oc, particles[1], pia, 1, 1, np_target_neutrals)
        end

        if pia.n_total[2] > threshold_ion
            merge_grid_based!(rng, mg_ions, particles[2], pia, 1, 2, species_data, phys_props)
        end

        if pia.n_total[3] > threshold_electrons
            merge_octree_N2_based!(rng, oc_electrons, particles[3], pia, 1, 3, np_target_electrons)
        end

        accelerate_constant_field_x!(particles[index_electron],
                                     pia, 1, index_electron, species_data,
                                     E_field, Δt)
        
        compute_props!(particles, pia, species_data, phys_props)
        write_netcdf(ds, phys_props, ts, sync_freq=1000)
    end
    close_netcdf(ds)

    sol = NCDataset(sol_path, "r")
 
    ndens_tot = n_dens_e + n_dens_ions + n_dens_neutrals

    # test overall number density conservation, we don't have it due to ionization
    ndens_conservation = true

    for t in 1:n_t
        if abs(sum(sol["ndens"][1, :, t]) - ndens_tot) / ndens_tot > 1.5e-14
            ndens_conservation = false
        end
    end

    @test ndens_conservation == false

    # test that we actually get ionization by looking at ion density at end of simulation
    @test sol["ndens"][1, 2, n_t] > n_dens_ions

    # test charge neutrality
    charge_neutrality = true

    for t in 1:n_t
        if abs(sol["ndens"][1, 2, t] - sol["ndens"][1, 3, t]) / sol["ndens"][1, 3, t] > 1.5e-14
            charge_neutrality = false
        end
    end

    @test charge_neutrality == true

    # test compared to ref solution
    ref_sol_path = joinpath(@__DIR__, "data", "ionization_Ar_no_es_v2.nc")
    ref_sol = NCDataset(ref_sol_path, "r")

    @test maximum(abs.(ref_sol["ndens"][1, :, 1:n_t] .- sol["ndens"][1, :, 1:n_t])) < 4 * eps()
    @test maximum(abs.(ref_sol["T"][1, :, 1:n_t] .- sol["T"][1, :, 1:n_t])) < 7.5e-14

    close(sol)
    rm(sol_path)
end