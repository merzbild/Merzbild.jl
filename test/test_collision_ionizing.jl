@testset "ionizing collision test" begin
    seed = 1234
    rng = StableRNG(seed)

    E_ion_eV = 15.76
    vy_neutral = 0.0

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, ["Ar", "e-"])

    # fill dummy data for Ar+e- interaction
    interaction_data_path = joinpath(@__DIR__, "..", "data", "vhs.toml")
    interaction_data = load_interaction_data(interaction_data_path, species_data, 1e-10, 1.0, 273.0)

    collision_data = CollisionData()

    particles_neutral = ParticleVector(1)
    Merzbild.update_particle_buffer_new_particle!(particles_neutral, 1)
    particles_neutral[1] = Particle(1.0, [0.0, vy_neutral, 0.0], [0.0, 1.0, -2.0])

    particles_electron = ParticleVector(2)
    Merzbild.update_particle_buffer_new_particle!(particles_electron, 1)
    particles_electron[1] = Particle(1.0, [5e6, 0.0, 0.0], [-3.0, 2.0, 4.0])

    # this will be the emitted electron
    Merzbild.update_particle_buffer_new_particle!(particles_electron, 2)
    particles_electron[2] = Particle(1.0, [0.0, 0.0, 0.0], [0.0, 1.0, -2.0])

    Merzbild.compute_g!(collision_data, particles_neutral[1], particles_electron[1])

    E_coll_electron_eV = 0.5 * collision_data.g^2 * Merzbild.e_mass_div_electron_volt  # convert to eV
    E_coll_eV = 0.5 * collision_data.g^2 * interaction_data[1,2].m_r * Merzbild.eV_J_inv

    collision_data.E_coll_eV = E_coll_eV
    collision_data.E_coll_electron_eV = E_coll_electron_eV

    @test abs(E_coll_electron_eV - E_coll_eV) < 1e-3
    @test E_coll_eV < E_coll_electron_eV  # m_n * m_e / (m_n + m_e) < m_e

    E_neutral_eV = Merzbild.eV_J_inv * 0.5 * species_data[1].mass * norm(particles_neutral[1].v)^2
    E_electron_pre_eV = Merzbild.eV_J_inv * 0.5 * species_data[2].mass * norm(particles_electron[1].v)^2

    @test abs(E_electron_pre_eV - E_coll_electron_eV) < 4e-12  # some minor differences due to constants possible

    E_total_pre_eV = E_neutral_eV + E_electron_pre_eV
    # ionization energy of argon is 15.76 eV
    Merzbild.compute_g_new_ionization!(collision_data, interaction_data[1,2], E_ion_eV, ElectronEnergySplitEqual)
    Merzbild.scatter_ionization_electrons!(rng, collision_data, particles_electron, 1, 2)

    # neutral particles don't change (we delete them later anyway but need the data to set the ions properties)
    @test maximum(abs.(particles_neutral[1].x - [0.0, 1.0, -2.0])) < 2*eps()
    @test maximum(abs.(particles_neutral[1].v - [0.0, vy_neutral, 0.0])) < 2*eps()
    @test abs.(particles_neutral[1].w - 1.0) < 2*eps()

    # positions of electrons don't change
    @test maximum(abs.(particles_electron[1].x - [-3.0, 2.0, 4.0])) < 2*eps()
    @test abs.(particles_electron[1].w - 1.0) < 2*eps()
    @test maximum(abs.(particles_electron[2].x - [0.0, 1.0, -2.0])) < 2*eps()
    @test abs.(particles_electron[2].w - 1.0) < 2*eps()

    E_electron1_post_eV = Merzbild.eV_J_inv * 0.5 * species_data[2].mass * norm(particles_electron[1].v)^2
    E_electron2_post_eV = Merzbild.eV_J_inv * 0.5 * species_data[2].mass * norm(particles_electron[2].v)^2

    @test abs(E_electron1_post_eV - E_electron2_post_eV) / E_electron1_post_eV < 2 * eps()

    E_total_post_eV = E_neutral_eV + E_ion_eV + E_electron1_post_eV + E_electron2_post_eV
    # TODO
    # test energies
    # test other split
    # test vy ! = 0.0
end