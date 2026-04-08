@testset "collide_2particles_vhs_equal_weight!" begin
    using StaticArrays
    using Random
    using StableRNGs

    # Load species and interaction data
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, "Ar")

    interaction_data_path = joinpath(@__DIR__, "..", "data", "vhs.toml")
    interaction_data = load_interaction_data(interaction_data_path, species_data)

    # Initialize random number generator
    seed = 1234
    rng = StableRNG(seed)

    # Create collision data and factors
    collision_data = CollisionData()
    collision_factors = CollisionFactors()

    # Create particle indexer array
    pia = ParticleIndexerArray(1, 1)  # 1 cell, 1 species
    pia.n_total[1] = 2  # Start with 2 particles

    # Initialize particles
    particles = ParticleVector(2)
    Merzbild.update_particle_buffer_new_particle!(particles, 1)
    Merzbild.update_particle_buffer_new_particle!(particles, 2)

    # Set up particles for testing
    particles[1] = Particle(1.0, SVector{3}(1.0, 0.0, 0.0), SVector{3}(0.0, 0.0, 0.0))
    particles[2] = Particle(1.0, SVector{3}(-1.0, 0.0, 0.0), SVector{3}(0.0, 0.0, 0.0))

    # Test case 1: Equal weight particles
    @testset "Equal weight particles" begin
        # Reset particles to equal weight
        particles[1].w = 1.0
        particles[2].w = 1.0
        
        # Store initial velocities and weights for verification
        v1_initial = particles[1].v
        v2_initial = particles[2].v
        w1_initial = particles[1].w
        w2_initial = particles[2].w

        # compute relative velocity and center-of-mass velocity
        Merzbild.compute_g!(collision_data, particles[1], particles[2])
        Merzbild.compute_com!(collision_data, interaction_data[1,1], particles[1], particles[2])
        
        # Compute initial kinetic energy
        m = species_data[1].mass
        ke_initial = 0.5 * m * (sum(v1_initial.^2) + sum(v2_initial.^2))
        
        # Perform collision
        Merzbild.collide_2particles_vhs_equal_weight!(rng, collision_data, collision_factors, 
                                                    interaction_data[1,1], particles[1], particles[2],
                                                    particles, particles, pia, 1, 1, 1)
        
        # Check that collision was performed
        @test collision_factors.n_coll_performed > 0
        
        # Check energy conservation
        ke_final = 0.5 * m * (sum(particles[1].v.^2) + sum(particles[2].v.^2))
        @test isapprox(ke_initial, ke_final, atol=1e-14)
        
        # Check that no particles were split (equal weight)
        @test collision_factors.n_eq_w_coll_performed > 0
        @test collision_factors.n_eq_w_coll_performed == collision_factors.n_coll_performed
        @test pia.n_total[1] == 2  # No new particles created
        
        # Check that weights remain unchanged
        @test particles[1].w == w1_initial
        @test particles[2].w == w2_initial
    end

    # Test case 2: Unequal weight particles (no splitting should occur)
    @testset "Unequal weight particles (w1>w2)" begin
        # Reset particles with unequal weights
        particles[1].w = 2.0
        particles[2].w = 1.0
        
        # Store initial state
        v1_initial = particles[1].v
        v2_initial = particles[2].v
        w1_initial = particles[1].w
        w2_initial = particles[2].w
        
        # Compute initial kinetic energy (weighted)
        m = species_data[1].mass
        ke_initial = 0.5 * m * (w1_initial * sum(v1_initial.^2) + w2_initial * sum(v2_initial.^2))

        # compute relative velocity and center-of-mass velocity
        Merzbild.compute_g!(collision_data, particles[1], particles[2])
        Merzbild.compute_com!(collision_data, interaction_data[1,1], particles[1], particles[2])
        
        # Perform collision
        Merzbild.collide_2particles_vhs_equal_weight!(rng, collision_data, collision_factors, 
                                                    interaction_data[1,1], particles[1], particles[2],
                                                    particles, particles, pia, 1, 1, 1)
        
        # Check that collision was performed
        @test collision_factors.n_coll_performed > 0
        @test collision_factors.n_eq_w_coll_performed == collision_factors.n_coll_performed
        
        # Check that no new particle was created (even though weights are unequal)
        @test pia.n_total[1] == 2
        
        # Check that weights remain unchanged
        @test particles[1].w == w1_initial
        @test particles[2].w == w2_initial
        
        # Check energy conservation (weighted)
        ke_final = 0.5 * m * (particles[1].w * sum(particles[1].v.^2) + 
                            particles[2].w * sum(particles[2].v.^2))
        @test isapprox(ke_initial, ke_final, atol=1e-14)
    end

    # Test case 3: Unequal weight particles (no splitting should occur)
    @testset "Unequal weight particles (w1<w2)" begin
        # Reset particles with unequal weights
        particles[1].w = 1.0
        particles[2].w = 2.0
        
        # Store initial state
        v1_initial = particles[1].v
        v2_initial = particles[2].v
        w1_initial = particles[1].w
        w2_initial = particles[2].w
        
        # Compute initial kinetic energy (weighted)
        m = species_data[1].mass
        ke_initial = 0.5 * m * (w1_initial * sum(v1_initial.^2) + w2_initial * sum(v2_initial.^2))

        # compute relative velocity and center-of-mass velocity
        Merzbild.compute_g!(collision_data, particles[1], particles[2])
        Merzbild.compute_com!(collision_data, interaction_data[1,1], particles[1], particles[2])
        
        # Perform collision
        Merzbild.collide_2particles_vhs_equal_weight!(rng, collision_data, collision_factors, 
                                                    interaction_data[1,1], particles[1], particles[2],
                                                    particles, particles, pia, 1, 1, 1)
        
        # Check that collision was performed
        @test collision_factors.n_coll_performed > 0
        @test collision_factors.n_eq_w_coll_performed == collision_factors.n_coll_performed
        
        # Check that no new particle was created (even though weights are unequal)
        @test pia.n_total[1] == 2
        
        # Check that weights remain unchanged
        @test particles[1].w == w1_initial
        @test particles[2].w == w2_initial
        
        # Check energy conservation (weighted)
        ke_final = 0.5 * m * (particles[1].w * sum(particles[1].v.^2) + 
                            particles[2].w * sum(particles[2].v.^2))
        @test isapprox(ke_initial, ke_final, atol=1e-14)
    end

    # Test case 4: Different species collision
    @testset "Different species collision" begin
        # Load two species
        species_data_2 = load_species_data(particles_data_path, ["Ar", "He"])
        interaction_data_2 = load_interaction_data(interaction_data_path, species_data_2)
        
        # Create particle indexer for 2 species
        pia_2 = ParticleIndexerArray(1, 2)
        pia_2.n_total[1] = 1  # Ar particles
        pia_2.n_total[2] = 1  # He particles
        
        # Initialize particles for both species
        particles_ar = ParticleVector(1)
        Merzbild.update_particle_buffer_new_particle!(particles_ar, 1)
        particles_ar[1] = Particle(1.0, SVector{3}(1.0, 0.0, 0.0), SVector{3}(0.0, 0.0, 0.0))
        
        particles_he = ParticleVector(1)
        Merzbild.update_particle_buffer_new_particle!(particles_he, 1)
        particles_he[1] = Particle(1.0, SVector{3}(-1.0, 0.0, 0.0), SVector{3}(0.0, 0.0, 0.0))
        
        # Store initial velocities and weights
        v_ar_initial = particles_ar[1].v
        v_he_initial = particles_he[1].v
        w_ar_initial = particles_ar[1].w
        w_he_initial = particles_he[1].w
        
        # Compute initial kinetic energy
        m_ar = species_data_2[1].mass
        m_he = species_data_2[2].mass
        ke_initial = 0.5 * (m_ar * sum(v_ar_initial.^2) + m_he * sum(v_he_initial.^2))

        # compute relative velocity and center-of-mass velocity
        Merzbild.compute_g!(collision_data, particles_ar[1], particles_he[1])
        Merzbild.compute_com!(collision_data, interaction_data_2[1,2], particles[1], particles[2])
        
        # Perform collision
        Merzbild.collide_2particles_vhs_equal_weight!(rng, collision_data, collision_factors, 
                                                    interaction_data_2[1,2], particles_ar[1], particles_he[1],
                                                    particles_ar, particles_he, pia_2, 1, 1, 2)
        
        # Check that collision was performed
        @test collision_factors.n_coll_performed > 0
        @test collision_factors.n_eq_w_coll_performed == collision_factors.n_coll_performed
        
        # Check that no new particles were created
        @test pia_2.n_total[1] == 1
        @test pia_2.n_total[2] == 1
        
        # Check that weights remain unchanged
        @test particles_ar[1].w == w_ar_initial
        @test particles_he[1].w == w_he_initial
        
        # Check energy conservation
        ke_final = 0.5 * (m_ar * sum(particles_ar[1].v.^2) + m_he * sum(particles_he[1].v.^2))
        @test isapprox(ke_initial, ke_final, atol=1e-14)
    end

    # Test case 5: Unequal weight particles with different species (no splitting should occur)
    @testset "Unequal weight particles with different species" begin
        # Load two species
        species_data_2 = load_species_data(particles_data_path, ["Ar", "He"])
        interaction_data_2 = load_interaction_data(interaction_data_path, species_data_2)
        
        # Create particle indexer for 2 species
        pia_2 = ParticleIndexerArray(1, 2)
        pia_2.n_total[1] = 1  # Ar particles
        pia_2.n_total[2] = 1  # He particles
        
        # Initialize particles for both species
        particles_ar = ParticleVector(1)
        Merzbild.update_particle_buffer_new_particle!(particles_ar, 1)
        particles_ar[1] = Particle(1.0, SVector{3}(1.0, 0.0, 0.0), SVector{3}(0.0, 0.0, 0.0))
        
        particles_he = ParticleVector(1)
        Merzbild.update_particle_buffer_new_particle!(particles_he, 1)
        particles_he[1] = Particle(1.0, SVector{3}(-1.0, 0.0, 0.0), SVector{3}(0.0, 0.0, 0.0))
        
        # Reset particles with unequal weights
        particles_ar[1].w = 2.0
        particles_he[1].w = 1.0
        
        # Store initial state
        v_ar_initial = particles_ar[1].v
        v_he_initial = particles_he[1].v
        w_ar_initial = particles_ar[1].w
        w_he_initial = particles_he[1].w
        
        # Compute initial kinetic energy (weighted)
        m_ar = species_data_2[1].mass
        m_he = species_data_2[2].mass
        ke_initial = 0.5 * (m_ar * w_ar_initial * sum(v_ar_initial.^2) + m_he * w_he_initial * sum(v_he_initial.^2))

        # compute relative velocity and center-of-mass velocity
        Merzbild.compute_g!(collision_data, particles_ar[1], particles_he[1])
        Merzbild.compute_com!(collision_data, interaction_data_2[1,2], particles_ar[1], particles_he[1])
        
        # Perform collision
        Merzbild.collide_2particles_vhs_equal_weight!(rng, collision_data, collision_factors, 
                                                    interaction_data_2[1,2], particles_ar[1], particles_he[1],
                                                    particles_ar, particles_he, pia_2, 1, 1, 2)
        
        # Check that collision was performed
        @test collision_factors.n_coll_performed > 0
        @test collision_factors.n_eq_w_coll_performed == collision_factors.n_coll_performed
        
        # Check that no new particles were created
        @test pia_2.n_total[1] == 1
        @test pia_2.n_total[2] == 1
        
        # Check that weights remain unchanged
        @test particles_ar[1].w == w_ar_initial
        @test particles_he[1].w == w_he_initial
        
        # Check energy conservation (weighted)
        ke_final = 0.5 * (m_ar * particles_ar[1].w * sum(particles_ar[1].v.^2) + 
                        m_he * particles_he[1].w * sum(particles_he[1].v.^2))
        @test isapprox(ke_initial, ke_final, atol=1e-14)
    end
end