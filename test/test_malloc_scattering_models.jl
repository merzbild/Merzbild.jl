@testset "malloc: collisions with VHS/VSS/hard sphere scattering models" begin
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, ["Ar", "He"])
    n_species = length(species_data)

    n_particles = [800, 800]
    T0_list = [3000.0, 360.0]
    Fnum = 5e12
    V = 1.0

    # the hard sphere cross-section is larger at thermal velocities, so a smaller timestep is used
    # to keep the number of particles created by SWPM bounded
    Δt_list = Dict("vhs.toml" => 2.5e-4, "vss.toml" => 2.5e-4, "hard_sphere.toml" => 2.5e-5)

    for interaction_file in ["vhs.toml", "vss.toml", "hard_sphere.toml"]
        rng = StableRNG(1234)

        interaction_data::Array{Interaction, 2} = load_interaction_data(joinpath(@__DIR__, "..", "data", interaction_file),
                                                                        species_data)

        # the particle vectors are over-allocated so that the particles created by SWPM
        # do not lead to a resize! (which would allocate for reasons unrelated to the collisions)
        particles = [ParticleVector(50000), ParticleVector(50000)]
        pia = ParticleIndexerArray([0, 0])
        for s in 1:n_species
            sample_particles_equal_weight!(rng, particles[s], pia, 1, s, n_particles[s], species_data[s].mass,
                                           T0_list[s], Fnum, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
        end

        collision_factors::Array{CollisionFactors, 3} = create_collision_factors_array(n_species)
        collision_factors_swpm::Array{CollisionFactorsSWPM, 3} = create_collision_factors_swpm_array(n_species)
        collision_data::CollisionData = CollisionData()

        estimate_sigma_g_w_max!(collision_factors, interaction_data, species_data, T0_list, Fnum)
        estimate_sigma_g_max!(collision_factors_swpm, interaction_data, species_data, T0_list)

        Δt = Δt_list[interaction_file]

        # warm-up so that everything is compiled before the allocations are measured
        for _ in 1:3
            ntc!(rng, collision_factors[1,1,1], collision_data, interaction_data, particles[1], pia, 1, 1, Δt, V)
            ntc!(rng, collision_factors[1,2,1], collision_data, interaction_data, particles[1], particles[2],
                 pia, 1, 1, 2, Δt, V)
            ntc_equal_weight!(rng, collision_factors[2,2,1], collision_data, interaction_data, particles[2],
                              pia, 1, 2, Δt, V)
            ntc_equal_weight!(rng, collision_factors[1,2,1], collision_data, interaction_data, particles[1], particles[2],
                              pia, 1, 1, 2, Δt, V)
            swpm!(rng, collision_factors_swpm[1,1,1], collision_data, interaction_data, particles[1], pia, 1, 1, 0.5, Δt, V)
        end

        for _ in 1:5
            bytes = @allocated ntc!(rng, collision_factors[1,1,1], collision_data, interaction_data, particles[1],
                                    pia, 1, 1, Δt, V)
            @test bytes == 0

            bytes = @allocated ntc!(rng, collision_factors[1,2,1], collision_data, interaction_data,
                                    particles[1], particles[2], pia, 1, 1, 2, Δt, V)
            @test bytes == 0

            bytes = @allocated ntc_equal_weight!(rng, collision_factors[2,2,1], collision_data, interaction_data,
                                                 particles[2], pia, 1, 2, Δt, V)
            @test bytes == 0

            bytes = @allocated ntc_equal_weight!(rng, collision_factors[1,2,1], collision_data, interaction_data,
                                                 particles[1], particles[2], pia, 1, 1, 2, Δt, V)
            @test bytes == 0

            bytes = @allocated swpm!(rng, collision_factors_swpm[1,1,1], collision_data, interaction_data,
                                     particles[1], pia, 1, 1, 0.5, Δt, V)
            @test bytes == 0
        end
    end
end
