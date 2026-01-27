@testset "ionizing collision test" begin

    function create_n(vx, vy)
        particles_n = ParticleVector(1)
        Merzbild.update_particle_buffer_new_particle!(particles_n, 1)
        particles_n[1] = Particle(1.0, [vx, vy, 0.0], [0.0, 1.0, -2.0])

        return particles_n
    end

    function create_e(vx, vy)
        particles_e = ParticleVector(2)
        Merzbild.update_particle_buffer_new_particle!(particles_e, 1)
        particles_e[1] = Particle(1.0, [vx, vy, 0.0], [-3.0, 2.0, 4.0])

        # this will be the emitted electron
        Merzbild.update_particle_buffer_new_particle!(particles_e, 2)
        particles_e[2] = Particle(1.0, [0.0, 0.0, 0.0], [0.0, 1.0, -2.0])

        return particles_e
    end

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

    particles_neutral = create_n(0.0, vy_neutral)
    particles_electron = create_e(5e6, 0.0)

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

    # test energy conservation; this is the simplest case with a stationary neutral
    @test abs(E_total_post_eV - E_total_pre_eV) / E_total_pre_eV < 1e-14

    # reset data and try out the other energy splitting (one electron takes all energy)
    particles_neutral = create_n(0.0, vy_neutral)
    particles_electron = create_e(5e6, 0.0)

    Merzbild.compute_g!(collision_data, particles_neutral[1], particles_electron[1])

    E_coll_electron_eV = 0.5 * collision_data.g^2 * Merzbild.e_mass_div_electron_volt  # convert to eV
    E_coll_eV = 0.5 * collision_data.g^2 * interaction_data[1,2].m_r * Merzbild.eV_J_inv

    collision_data.E_coll_eV = E_coll_eV
    collision_data.E_coll_electron_eV = E_coll_electron_eV

    E_neutral_eV = Merzbild.eV_J_inv * 0.5 * species_data[1].mass * norm(particles_neutral[1].v)^2
    E_electron_pre_eV = Merzbild.eV_J_inv * 0.5 * species_data[2].mass * norm(particles_electron[1].v)^2
    E_total_pre_eV = E_neutral_eV + E_electron_pre_eV

    Merzbild.compute_g_new_ionization!(collision_data, interaction_data[1,2], E_ion_eV, ElectronEnergySplitZeroE)
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

    # some loss of precision due to conversions, etc.
    @test abs(E_electron1_post_eV - (E_electron_pre_eV - E_ion_eV)) / E_electron1_post_eV < 1.25e-14
    @test abs(E_electron2_post_eV) < 2 * eps()

    E_total_post_eV = E_neutral_eV + E_ion_eV + E_electron1_post_eV + E_electron2_post_eV

    # test energy conservation; this is the simplest case with a stationary neutral
    @test abs(E_total_post_eV - E_total_pre_eV) / E_total_pre_eV < 1e-14


    # now we switch to stationary electron and moving neutral
    vy_neutral = 5e6

    particles_neutral = create_n(0.0, vy_neutral)
    particles_electron = create_e(0.0, 0.0)

    Merzbild.compute_g!(collision_data, particles_neutral[1], particles_electron[1])

    E_coll_electron_eV = 0.5 * collision_data.g^2 * Merzbild.e_mass_div_electron_volt  # convert to eV
    E_coll_eV = 0.5 * collision_data.g^2 * interaction_data[1,2].m_r * Merzbild.eV_J_inv

    collision_data.E_coll_eV = E_coll_eV
    collision_data.E_coll_electron_eV = E_coll_electron_eV

    @test abs(E_coll_electron_eV - E_coll_eV) < 1e-3
    @test E_coll_eV < E_coll_electron_eV  # m_n * m_e / (m_n + m_e) < m_e

    E_neutral_eV = Merzbild.eV_J_inv * 0.5 * species_data[1].mass * norm(particles_neutral[1].v)^2
    E_electron_pre_eV = Merzbild.eV_J_inv * 0.5 * species_data[2].mass * norm(particles_electron[1].v)^2

    # energy is very different from the energy of the electron itself (since that assumes stationary neutral)
    @test abs(E_electron_pre_eV - E_coll_electron_eV) > 50

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

    # test energy conservation
    @test abs(E_total_post_eV - E_total_pre_eV) / E_total_pre_eV < 1e-14

    # final test: both neutral and electron moving
    vx_neutral = -4e3
    vy_neutral = 5e6
    vx_electron = 5e5
    vy_electron = -4e5

    particles_neutral = create_n(vx_neutral, vy_neutral)
    particles_electron = create_e(vx_electron, vy_electron)

    Merzbild.compute_g!(collision_data, particles_neutral[1], particles_electron[1])

    E_coll_electron_eV = 0.5 * collision_data.g^2 * Merzbild.e_mass_div_electron_volt  # convert to eV
    E_coll_eV = 0.5 * collision_data.g^2 * interaction_data[1,2].m_r * Merzbild.eV_J_inv

    collision_data.E_coll_eV = E_coll_eV
    collision_data.E_coll_electron_eV = E_coll_electron_eV

    @test E_coll_eV < E_coll_electron_eV  # m_n * m_e / (m_n + m_e) < m_e

    E_neutral_eV = Merzbild.eV_J_inv * 0.5 * species_data[1].mass * norm(particles_neutral[1].v)^2
    E_electron_pre_eV = Merzbild.eV_J_inv * 0.5 * species_data[2].mass * norm(particles_electron[1].v)^2

    # energy is very different from the energy of the electron itself (since that assumes stationary neutral)
    @test abs(E_electron_pre_eV - E_coll_electron_eV) > 1.5

    E_total_pre_eV = E_neutral_eV + E_electron_pre_eV
    # ionization energy of argon is 15.76 eV
    Merzbild.compute_g_new_ionization!(collision_data, interaction_data[1,2], E_ion_eV, ElectronEnergySplitEqual)
    Merzbild.scatter_ionization_electrons!(rng, collision_data, particles_electron, 1, 2)

    # neutral particles don't change (we delete them later anyway but need the data to set the ions properties)
    @test maximum(abs.(particles_neutral[1].x - [0.0, 1.0, -2.0])) < 2*eps()
    @test maximum(abs.(particles_neutral[1].v - [vx_neutral, vy_neutral, 0.0])) < 2*eps()
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

    # test energy conservation
    @test abs(E_total_post_eV - E_total_pre_eV) / E_total_pre_eV < 1e-14
end