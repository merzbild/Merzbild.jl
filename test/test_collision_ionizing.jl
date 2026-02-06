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

    function create_ion()
        particles_i = ParticleVector(1)
        Merzbild.update_particle_buffer_new_particle!(particles_i, 1)
        particles_i[1] = Particle(1.0, [0.0, 0.0, 0.0], [-3.0, 2.0, 4.0])
        return particles_i
    end

    seed = 1234
    rng = StableRNG(seed)

    E_ion_eV = 15.76
    vy_neutral = 0.0

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, ["Ar", "e-"])
    m_ratio = species_data[2].mass / species_data[1].mass

    # fill dummy data for Ar+e- interaction
    interaction_data_path = joinpath(@__DIR__, "..", "data", "vhs.toml")
    interaction_data = load_interaction_data(interaction_data_path, species_data, 1e-10, 1.0, 273.0)

    collision_data = CollisionData()

    particles_neutral = create_n(0.0, vy_neutral)
    particles_electron = create_e(3e6, 0.0)
    particles_ion = create_ion()

    Merzbild.compute_g!(collision_data, particles_neutral[1], particles_electron[1])
    Merzbild.compute_com!(collision_data, interaction_data[1,2], particles_neutral[1], particles_electron[1])

    E_coll_eV = 0.5 * collision_data.g^2 * interaction_data[1,2].m_r * Merzbild.eV_J_inv

    collision_data.E_coll_eV = E_coll_eV

    E_neutral_eV = Merzbild.eV_J_inv * 0.5 * species_data[1].mass * norm(particles_neutral[1].v)^2
    E_electron_pre_eV = Merzbild.eV_J_inv * 0.5 * species_data[2].mass * norm(particles_electron[1].v)^2

    E_total_pre_eV = E_neutral_eV + E_electron_pre_eV
    # ionization energy of argon is 15.76 eV
    Merzbild.compute_g_new_ionization!(collision_data, interaction_data[1,2], E_ion_eV, ElectronEnergySplitEqual)

    Merzbild.scatter_ionization_electrons_and_ion!(rng, collision_data, particles_electron, particles_ion, 1, 2, 1, m_ratio)

    # neutral particles don't change (we delete them later anyway but need the data to set the ions properties)
    @test maximum(abs.(particles_neutral[1].x - [0.0, 1.0, -2.0])) < 2*eps()
    @test maximum(abs.(particles_neutral[1].v - [0.0, vy_neutral, 0.0])) < 2*eps()
    @test abs.(particles_neutral[1].w - 1.0) < 2*eps()

    # positions of electrons don't change
    @test maximum(abs.(particles_electron[1].x - [-3.0, 2.0, 4.0])) < 2*eps()
    @test abs.(particles_electron[1].w - 1.0) < 2*eps()
    @test maximum(abs.(particles_electron[2].x - [0.0, 1.0, -2.0])) < 2*eps()
    @test abs.(particles_electron[2].w - 1.0) < 2*eps()

    # energy wrt absolute frame of reference
    E_electron1_post_eV = Merzbild.eV_J_inv * 0.5 * species_data[2].mass * norm(particles_electron[1].v)^2
    E_electron2_post_eV = Merzbild.eV_J_inv * 0.5 * species_data[2].mass * norm(particles_electron[2].v)^2

    @test norm(particles_ion[1].v) > 0.0

    # assume ion mass ~ neutral mass
    E_i_eV = Merzbild.eV_J_inv * 0.5 * species_data[1].mass * norm(particles_ion[1].v)^2

    E_total_post_eV = E_i_eV + E_ion_eV + E_electron1_post_eV + E_electron2_post_eV
    # test energy conservation; this is the simplest case with a stationary neutral
    # error on the order of 1e-5 since we assume mass(ion) = mass(electron) to simplify collision mechanics
    # for argon this leads to an error on the order of 1.25e-5
    @test abs(E_total_post_eV - E_total_pre_eV) / E_total_pre_eV < 1e-5

    # now we shift to v_com frame of reference and neet to check equal energies
    E_electron1_post_eV = Merzbild.eV_J_inv * 0.5 * species_data[2].mass * norm(particles_electron[1].v - collision_data.v_com)^2
    E_electron2_post_eV = Merzbild.eV_J_inv * 0.5 * species_data[2].mass * norm(particles_electron[2].v - collision_data.v_com)^2
    @test abs(E_electron1_post_eV - E_electron2_post_eV) / E_electron1_post_eV < 2 * eps()

    # Testcase 2
    # reset data and try out the other energy splitting (one electron takes all energy)
    particles_neutral = create_n(0.0, vy_neutral)
    particles_electron = create_e(-6e6, 0.0)

    Merzbild.compute_g!(collision_data, particles_neutral[1], particles_electron[1])
    Merzbild.compute_com!(collision_data, interaction_data[1,2], particles_neutral[1], particles_electron[1])

    E_coll_eV = 0.5 * collision_data.g^2 * interaction_data[1,2].m_r * Merzbild.eV_J_inv

    collision_data.E_coll_eV = E_coll_eV

    E_neutral_eV = Merzbild.eV_J_inv * 0.5 * species_data[1].mass * norm(particles_neutral[1].v)^2
    E_electron_pre_eV = Merzbild.eV_J_inv * 0.5 * species_data[2].mass * norm(particles_electron[1].v)^2
    E_total_pre_eV = E_neutral_eV + E_electron_pre_eV

    Merzbild.compute_g_new_ionization!(collision_data, interaction_data[1,2], E_ion_eV, ElectronEnergySplitZeroE)
    Merzbild.scatter_ionization_electrons_and_ion!(rng, collision_data, particles_electron, particles_ion, 1, 2, 1, m_ratio)

    # neutral particles don't change (we delete them later anyway but need the data to set the ions properties)
    @test maximum(abs.(particles_neutral[1].x - [0.0, 1.0, -2.0])) < 2*eps()
    @test maximum(abs.(particles_neutral[1].v - [0.0, vy_neutral, 0.0])) < 2*eps()
    @test abs.(particles_neutral[1].w - 1.0) < 2*eps()

    # positions of electrons don't change
    @test maximum(abs.(particles_electron[1].x - [-3.0, 2.0, 4.0])) < 2*eps()
    @test abs.(particles_electron[1].w - 1.0) < 2*eps()
    @test maximum(abs.(particles_electron[2].x - [0.0, 1.0, -2.0])) < 2*eps()
    @test abs.(particles_electron[2].w - 1.0) < 2*eps()

    @test norm(particles_ion[1].v) > 0.0

    E_electron1_post_eV = Merzbild.eV_J_inv * 0.5 * species_data[2].mass * norm(particles_electron[1].v)^2
    E_electron2_post_eV = Merzbild.eV_J_inv * 0.5 * species_data[2].mass * norm(particles_electron[2].v)^2
    @test E_electron1_post_eV > 0
    @test E_electron2_post_eV > 0

    # assume ion mass ~ neutral mass
    E_i_eV = Merzbild.eV_J_inv * 0.5 * species_data[1].mass * norm(particles_ion[1].v)^2

    E_total_post_eV = E_i_eV + E_ion_eV + E_electron1_post_eV + E_electron2_post_eV
    @test abs(E_total_post_eV - E_total_pre_eV) / E_total_pre_eV < 1.7e-5

    # switch to c-o-m frame of reference
    E_electron1_post_eV = Merzbild.eV_J_inv * 0.5 * species_data[2].mass * norm(particles_electron[1].v - collision_data.v_com)^2
    E_electron2_post_eV = Merzbild.eV_J_inv * 0.5 * species_data[2].mass * norm(particles_electron[2].v - collision_data.v_com)^2

    @test E_electron1_post_eV > 1
    @test abs(E_electron2_post_eV) < 2 * eps()

    # Testcase 3
    # now we switch to stationary electron and moving neutral
    vy_neutral = 4e6

    particles_neutral = create_n(0.0, vy_neutral)
    particles_electron = create_e(0.0, 0.0)

    Merzbild.compute_g!(collision_data, particles_neutral[1], particles_electron[1])
    Merzbild.compute_com!(collision_data, interaction_data[1,2], particles_neutral[1], particles_electron[1])

    E_coll_eV = 0.5 * collision_data.g^2 * interaction_data[1,2].m_r * Merzbild.eV_J_inv

    collision_data.E_coll_eV = E_coll_eV

    E_neutral_eV = Merzbild.eV_J_inv * 0.5 * species_data[1].mass * norm(particles_neutral[1].v)^2
    E_electron_pre_eV = Merzbild.eV_J_inv * 0.5 * species_data[2].mass * norm(particles_electron[1].v)^2

    E_total_pre_eV = E_neutral_eV + E_electron_pre_eV
    # ionization energy of argon is 15.76 eV
    Merzbild.compute_g_new_ionization!(collision_data, interaction_data[1,2], E_ion_eV, ElectronEnergySplitEqual)
    Merzbild.scatter_ionization_electrons_and_ion!(rng, collision_data, particles_electron, particles_ion, 1, 2, 1, m_ratio)

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
    E_i_eV = Merzbild.eV_J_inv * 0.5 * species_data[1].mass * norm(particles_ion[1].v)^2
    E_total_post_eV = E_i_eV + E_ion_eV + E_electron1_post_eV + E_electron2_post_eV

    # test energy conservation
    @test abs(E_total_post_eV - E_total_pre_eV) / E_total_pre_eV < 1.5e-5

    # switch to c-o-m to test energy split
    E_electron1_post_eV = Merzbild.eV_J_inv * 0.5 * species_data[2].mass * norm(particles_electron[1].v - collision_data.v_com)^2
    E_electron2_post_eV = Merzbild.eV_J_inv * 0.5 * species_data[2].mass * norm(particles_electron[2].v - collision_data.v_com)^2
    @test abs(E_electron1_post_eV - E_electron2_post_eV) / E_electron1_post_eV < 2 * eps()


    # final test: both neutral and electron moving
    vx_neutral = -4e3
    vy_neutral = 5e6
    vx_electron = 5e5
    vy_electron = -4e5

    particles_neutral = create_n(vx_neutral, vy_neutral)
    particles_electron = create_e(vx_electron, vy_electron)

    Merzbild.compute_g!(collision_data, particles_neutral[1], particles_electron[1])
    Merzbild.compute_com!(collision_data, interaction_data[1,2], particles_neutral[1], particles_electron[1])

    E_coll_eV = 0.5 * collision_data.g^2 * interaction_data[1,2].m_r * Merzbild.eV_J_inv

    collision_data.E_coll_eV = E_coll_eV

    E_neutral_eV = Merzbild.eV_J_inv * 0.5 * species_data[1].mass * norm(particles_neutral[1].v)^2
    E_electron_pre_eV = Merzbild.eV_J_inv * 0.5 * species_data[2].mass * norm(particles_electron[1].v)^2

    E_total_pre_eV = E_neutral_eV + E_electron_pre_eV
    # ionization energy of argon is 15.76 eV
    Merzbild.compute_g_new_ionization!(collision_data, interaction_data[1,2], E_ion_eV, ElectronEnergySplitEqual)
    Merzbild.scatter_ionization_electrons_and_ion!(rng, collision_data, particles_electron, particles_ion, 1, 2, 1, m_ratio)

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

    E_i_eV = Merzbild.eV_J_inv * 0.5 * species_data[1].mass * norm(particles_ion[1].v)^2
    E_total_post_eV = E_i_eV + E_ion_eV + E_electron1_post_eV + E_electron2_post_eV

    # test energy conservation
    @test abs(E_total_post_eV - E_total_pre_eV) / E_total_pre_eV < 1.5e-5

    # test energy split
    E_electron1_post_eV = Merzbild.eV_J_inv * 0.5 * species_data[2].mass * norm(particles_electron[1].v - collision_data.v_com)^2
    E_electron2_post_eV = Merzbild.eV_J_inv * 0.5 * species_data[2].mass * norm(particles_electron[2].v - collision_data.v_com)^2

    @test abs(E_electron1_post_eV - E_electron2_post_eV) / E_electron1_post_eV < 2 * eps()


    # final test
    # to verify that error in energy conservation is due to mass approximation
    # we create an artificial argon species with a mass 10000 times that of argon
    # this should improve conservation by 4 orders of magnitude at least

    species_data[1] = Species("Ar", 66.3e-23, 0.0, 0.0)

    m_ratio = species_data[2].mass / species_data[1].mass

    # fill dummy data for Ar+e- interaction
    interaction_data_path = joinpath(@__DIR__, "..", "data", "vhs.toml")
    interaction_data = load_interaction_data(interaction_data_path, species_data, 1e-10, 1.0, 273.0)

    collision_data = CollisionData()

    particles_neutral = create_n(0.0, vy_neutral)
    particles_electron = create_e(3e6, 0.0)
    particles_ion = create_ion()

    Merzbild.compute_g!(collision_data, particles_neutral[1], particles_electron[1])
    Merzbild.compute_com!(collision_data, interaction_data[1,2], particles_neutral[1], particles_electron[1])

    E_coll_eV = 0.5 * collision_data.g^2 * interaction_data[1,2].m_r * Merzbild.eV_J_inv

    collision_data.E_coll_eV = E_coll_eV

    E_neutral_eV = Merzbild.eV_J_inv * 0.5 * species_data[1].mass * norm(particles_neutral[1].v)^2
    E_electron_pre_eV = Merzbild.eV_J_inv * 0.5 * species_data[2].mass * norm(particles_electron[1].v)^2

    E_total_pre_eV = E_neutral_eV + E_electron_pre_eV
    # ionization energy of argon is 15.76 eV
    Merzbild.compute_g_new_ionization!(collision_data, interaction_data[1,2], E_ion_eV, ElectronEnergySplitEqual)

    Merzbild.scatter_ionization_electrons_and_ion!(rng, collision_data, particles_electron, particles_ion, 1, 2, 1, m_ratio)

    # neutral particles don't change (we delete them later anyway but need the data to set the ions properties)
    @test maximum(abs.(particles_neutral[1].x - [0.0, 1.0, -2.0])) < 2*eps()
    @test maximum(abs.(particles_neutral[1].v - [0.0, vy_neutral, 0.0])) < 2*eps()
    @test abs.(particles_neutral[1].w - 1.0) < 2*eps()

    # positions of electrons don't change
    @test maximum(abs.(particles_electron[1].x - [-3.0, 2.0, 4.0])) < 2*eps()
    @test abs.(particles_electron[1].w - 1.0) < 2*eps()
    @test maximum(abs.(particles_electron[2].x - [0.0, 1.0, -2.0])) < 2*eps()
    @test abs.(particles_electron[2].w - 1.0) < 2*eps()

    # energy wrt absolute frame of reference
    E_electron1_post_eV = Merzbild.eV_J_inv * 0.5 * species_data[2].mass * norm(particles_electron[1].v)^2
    E_electron2_post_eV = Merzbild.eV_J_inv * 0.5 * species_data[2].mass * norm(particles_electron[2].v)^2

    @test norm(particles_ion[1].v) > 0.0

    # assume ion mass ~ neutral mass
    E_i_eV = Merzbild.eV_J_inv * 0.5 * species_data[1].mass * norm(particles_ion[1].v)^2

    E_total_post_eV = E_i_eV + E_ion_eV + E_electron1_post_eV + E_electron2_post_eV
    # test energy conservation with the heavier neutral, should work better
    @test abs(E_total_post_eV - E_total_pre_eV) / E_total_pre_eV < 1.4e-9
end