@testset "2species elastic" begin

    seed = 1234
    rng = StableRNG(seed)

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, ["Ar", "He"])

    interaction_data_path = joinpath(@__DIR__, "..", "data", "vhs.toml")
    interaction_data::Array{Interaction, 2} = load_interaction_data(interaction_data_path, species_data)
    n_species = length(species_data)

    n_t = 800

    n_particles_Ar = 400
    n_particles_He = 4000
    T0_Ar = 3000.0
    T0_He = 360.0  # equilibrium T = 600.0K
    T0_list = [T0_Ar, T0_He]
    Fnum = 5e12

    n_Ar = Fnum * n_particles_Ar
    n_He = Fnum * n_particles_He

    T_eq = (n_Ar * T0_Ar + n_He * T0_He) / (n_Ar + n_He)

    particles::Vector{ParticleVector} = [ParticleVector(n_particles_Ar), ParticleVector(n_particles_He)]

    pia = ParticleIndexerArray([0, 0])

    sample_particles_equal_weight!(rng, particles[1], pia, 1, 1, 
                                   n_particles_Ar, species_data[1].mass, T0_Ar, Fnum, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
    sample_particles_equal_weight!(rng, particles[2], pia, 1, 2, 
                                   n_particles_He, species_data[2].mass, T0_He, Fnum, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)


    phys_props::PhysProps = PhysProps(1, 2, [], Tref=T0_Ar)
    compute_props!(particles, pia, species_data, phys_props)
    
    sol_path = joinpath(@__DIR__, "data", "tmp_2species_elastic.nc")
    ds = NCDataHolder(sol_path, species_data, phys_props)
    write_netcdf_phys_props(ds, phys_props, 0)

    collision_factors::Array{CollisionFactors, 3} = create_collision_factors_array(n_species)
    collision_data::CollisionData = CollisionData()

    estimate_sigma_g_w_max!(collision_factors, interaction_data, species_data, T0_list, Fnum)

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

        compute_props!(particles, pia, species_data, phys_props)
        write_netcdf_phys_props(ds, phys_props, ts)
    end

    close_netcdf(ds)

    ref_sol_path = joinpath(@__DIR__, "data", "2species_seed1234.nc")
    ref_sol = NCDataset(ref_sol_path, "r")
    sol = NCDataset(sol_path, "r")

    ref_T = ref_sol["T"]
    sol_T = sol["T"]

    @test maximum(abs.(sol["ndens"][1, 1, :] .- n_Ar)) / n_Ar < eps()
    @test maximum(abs.(sol["ndens"][1, 2, :] .- n_He)) / n_He < eps()

    @test maximum(sol["np"][1,1,:]) == n_particles_Ar
    @test minimum(sol["np"][1,1,:]) == n_particles_Ar

    @test maximum(sol["np"][1,2,:]) == n_particles_He
    @test minimum(sol["np"][1,2,:]) == n_particles_He

    for species in 1:2
        diff = abs.(ref_T[1, species, :] - sol_T[1, species, :])
        @test maximum(diff) < 1.2e-13
    end

    for species in 1:2
        @test abs(sol_T[1, species, end] - T_eq) / T_eq < 0.12
    end

    close(sol)
    close(ref_sol)

    rm(sol_path)
end