@testset "ionizing collision test" begin

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, ["Ar", "e-"])

    # fill dummy data for Ar+e- interaction
    interaction_data_path = joinpath(@__DIR__, "..", "data", "vhs.toml")
    interaction_data = load_interaction_data(interaction_data_path, species_data, 1e-10, 1.0, 273.0)

    collision_data = CollisionData()

    particles_neutral = ParticleVector(1)
    Merzbild.update_particle_buffer_new_particle!(particles_neutral, 1)
    particles_neutral[1] = Particle(1.0, [0.0, 0.0, 0.0], [0.0, 1.0, -2.0])

    particles_electron = ParticleVector(2)
    Merzbild.update_particle_buffer_new_particle!(particles_electron, 1)
    particles_electron[1] = Particle(1.0, [5e6, 0.0, 0.0], [-3.0, 2.0, 4.0])

    # this will be the emitted electron
    Merzbild.update_particle_buffer_new_particle!(particles_electron, 2)
    particles_electron[2] = Particle(1.0, [0.0, 0.0, 0.0], [0.0, 1.0, -2.0])

    Merzbild.compute_g!(collision_data, particles_neutral[1], particles_electron[1])

    E_coll_electron_eV = 0.5 * collision_data.g^2 * Merzbild.e_mass_div_electron_volt  # convert to eV
    E_coll = 0.5 * collision_data.g^2 * interaction_data[1,2].m_r * Merzbild.eV_J_inv

    @test abs(E_coll_electron_eV - E_coll) < 1e-3
    @test E_coll < E_coll_electron_eV  # m_n * m_e / (m_n + m_e) < m_e

    # E_total = 0.5 * ()
    # println("e- only: $E_coll_electron_eV exact: $E_coll")
    # compute_g_new_ionization!(collision_data, interaction, E_i, energy_splitting)
    # scatter_ionization_electrons!(rng, collision_data, particles_electron, i1, i2)

end