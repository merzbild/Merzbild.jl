@testset "2species elastic" begin


    seed = 1234
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    species_list::Vector{Species} = load_species_list("data/particles.toml", ["Ar", "He"])
    interaction_data::Array{Interaction, 2} = load_interaction_data("data/vhs.toml", species_list)
    n_species = length(species_list)

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

    particles::Vector{Vector{Particle}} = [Vector{Particle}(undef, n_particles_Ar), Vector{Particle}(undef, n_particles_He)]
    sample_particles_equal_weight!(rng, particles[1], n_particles_Ar, T0_Ar, species_list[1].mass, Fnum, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
    sample_particles_equal_weight!(rng, particles[2], n_particles_He, T0_He, species_list[2].mass, Fnum, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)

    particle_indexer::Array{ParticleIndexer, 2} = Array{ParticleIndexer, 2}(undef, 2, 1)
    particle_indexer[1,1] = create_particle_indexer(n_particles_Ar)
    particle_indexer[2,1] = create_particle_indexer(n_particles_He)

    phys_props::PhysProps = create_props(1, 2, [], Tref=T0_Ar)
    compute_props!(phys_props, particle_indexer, particles, species_list)
    
    ds = create_netcdf_phys_props("test/data/tmp_2species_elastic.nc", phys_props, species_list)
    write_netcdf_phys_props(ds, phys_props, 0)

    collision_factors::Array{CollisionFactors, 2} = create_collision_factors(n_species)
    collision_data::CollisionData = create_collision_data()

    estimate_sigma_g_w_max!(collision_factors, interaction_data, species_list, T0_list, Fnum)

    Δt = 2.5e-3
    V = 1.0

    for ts in 1:n_t
        for s2 in 1:n_species
            for s1 in s2:n_species
                if (s1 == s2)
                    ntc!(rng, collision_factors[s1,s1], particle_indexer[s1,1], collision_data, interaction_data[s1,s1], particles[s1],
                    Δt, V)
                else
                    ntc!(rng, collision_factors[s1,s2], particle_indexer[s1,1], particle_indexer[s2,1],
                    collision_data, interaction_data[s1,s2], particles[s1], particles[s2],
                    Δt, V)
                end
            end
        end

        # compute_props!(phys_props, particle_indexer, particles, species_list)
        compute_props_sorted_without_moments!(phys_props, particle_indexer, particles, species_list)
        write_netcdf_phys_props(ds, phys_props, ts)
    end

    close(ds)


    ref_sol = NCDataset("test/data/2species_seed1234.nc", "r")
    sol = NCDataset("test/data/tmp_2species_elastic.nc", "r")

    ref_T = ref_sol["T"]
    sol_T = sol["T"]

    @test maximum(abs.(ref_sol["ndens"][1, 1, :] .- n_Ar)) / n_Ar < eps()
    @test maximum(abs.(ref_sol["ndens"][1, 2, :] .- n_He)) / n_He < eps()


    @test maximum(sol["np"][1,1,:]) == n_particles_Ar
    @test minimum(sol["np"][1,1,:]) == n_particles_Ar

    @test maximum(sol["np"][1,2,:]) == n_particles_He
    @test minimum(sol["np"][1,2,:]) == n_particles_He

    for species in 1:2
        diff = abs.(ref_T[1, species, :] - sol_T[1, species, :])
        @test maximum(diff) < eps()
    end

    for species in 1:2
        @test abs(ref_T[1, species, end] - T_eq) / T_eq < 0.1
    end

    close(sol)
    close(ref_sol)

    rm("test/data/tmp_2species_elastic.nc")
end