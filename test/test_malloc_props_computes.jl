@testset "malloc: props computes and averaging" begin

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

    particles = [ParticleVector(n_particles_Ar), ParticleVector(n_particles_He)]

    pia = ParticleIndexerArray([0, 0])

    sample_particles_equal_weight!(rng, particles[1], pia, 1, 1, 
                                   n_particles_Ar, species_data[1].mass, T0_Ar, Fnum, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
    sample_particles_equal_weight!(rng, particles[2], pia, 1, 2, 
                                   n_particles_He, species_data[2].mass, T0_He, Fnum, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)


    phys_props = PhysProps(1, 2)
    phys_props_avg = PhysProps(1, 2)
    compute_props_sorted!(particles, pia, species_data, phys_props)
    avg_props!(phys_props_avg, phys_props, 8)

    collision_factors::Array{CollisionFactors, 3} = create_collision_factors_array(n_species)
    collision_data::CollisionData = CollisionData()

    estimate_sigma_g_w_max!(collision_factors, interaction_data, species_data, T0_list, Fnum)

    Δt = 2.5e-3
    V = 1.0

    for ts in 1:2
        compute_props_sorted!(particles, pia, species_data, phys_props)
    end

    for ts in 1:5
        bytes_props = @allocated compute_props_sorted!(particles, pia, species_data, phys_props)
        @test bytes_props == 0
        bytes_props_avg = @allocated avg_props!(phys_props_avg, phys_props, 8)
        @test bytes_props_avg == 0
    end
end