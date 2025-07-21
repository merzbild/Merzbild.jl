@testset "2species elastic variable weight" begin

    seed = 1234
    rng = StableRNG(seed)

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, ["Ar", "He"])

    interaction_data_path = joinpath(@__DIR__, "..", "data", "vhs.toml")
    interaction_data::Array{Interaction, 2} = load_interaction_data(interaction_data_path, species_data)
    n_species = length(species_data)

    oc = OctreeN2Merge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel, max_Nbins=6000)

    n_t = 800

    T0_Ar = 3000.0
    T0_He = 360.0  # equilibrium T = 600.0K
    T0_list = [T0_Ar, T0_He]

    n_Ar = 2e15
    n_He = 2e16

    Fnum_Ar = 5e11
    Fnum_He = 5e12

    n_particles_Ar = round(Int32, n_Ar / Fnum_Ar)
    n_particles_He = round(Int32, n_He / Fnum_He)

    threshold_Ar = round(Int32, n_particles_Ar * 1.2)
    threshold_He = round(Int32, n_particles_He * 1.2)

    T_eq = (n_Ar * T0_Ar + n_He * T0_He) / (n_Ar + n_He)

    particles::Vector{ParticleVector} = [ParticleVector(n_particles_Ar), ParticleVector(n_particles_He)]

    pia = ParticleIndexerArray([0, 0])

    sample_particles_equal_weight!(rng, particles[1], pia, 1, 1, 
                                   n_particles_Ar, species_data[1].mass, T0_Ar, Fnum_Ar, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
    sample_particles_equal_weight!(rng, particles[2], pia, 1, 2, 
                                   n_particles_He, species_data[2].mass, T0_He, Fnum_He, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)

    phys_props::PhysProps = PhysProps(1, 2, [], Tref=T0_Ar)
    compute_props!(particles, pia, species_data, phys_props)
    
    sol_path = joinpath(@__DIR__, "data", "tmp_2species_varweight_elastic.nc")
    ds = NCDataHolder(sol_path, species_data, phys_props)
    write_netcdf_phys_props(ds, phys_props, 0)

    collision_factors::Array{CollisionFactors, 3} = create_collision_factors_array(n_species)
    collision_data::CollisionData = CollisionData()

    estimate_sigma_g_w_max!(collision_factors, interaction_data, species_data, T0_list, [Fnum_Ar, Fnum_He])

    Δt = 2.5e-3
    V = 1.0

    for ts in 1:n_t
        for s2 in 1:n_species
            for s1 in s2:n_species
                if (s1 == s2)
                    ntc!(rng, collision_factors[s1,s1,1], collision_data, interaction_data, particles[s1], pia, 1, s1, Δt, V)
                else
                    ntc!(rng, collision_factors[s1,s2,1], collision_data, interaction_data, particles[s1], particles[s2],
                         pia, 1, s1, s2, Δt, V)
                end
            end
        end

        if pia.indexer[1,1].n_local > threshold_Ar
            merge_octree_N2_based!(rng, oc, particles[1], pia, 1, 1, n_particles_Ar)
        end

        if pia.indexer[1,2].n_local > threshold_He
            merge_octree_N2_based!(rng, oc, particles[2], pia, 1, 2, n_particles_He)
        end

        compute_props!(particles, pia, species_data, phys_props)
        write_netcdf_phys_props(ds, phys_props, ts)
    end

    close_netcdf(ds)

    ref_sol_path = joinpath(@__DIR__, "data", "2species_varweight_octree_seed1234.nc")
    ref_sol = NCDataset(ref_sol_path, "r")
    sol = NCDataset(sol_path, "r")

    ref_T = ref_sol["T"]
    sol_T = sol["T"]

    @test maximum(abs.(sol["ndens"][1, 1, :] .- n_Ar)) / n_Ar < 2e-15
    @test maximum(abs.(sol["ndens"][1, 2, :] .- n_He)) / n_He < 6e-15

    for species in 1:2
        diff = abs.(ref_T[1, species, :] - sol_T[1, species, :])
        @test maximum(diff) < 9.3e-13
    end

    for species in 1:2
        @test abs(sol_T[1, species, end] - T_eq) / T_eq < 0.055
    end

    close(sol)
    close(ref_sol)

    rm(sol_path)
end