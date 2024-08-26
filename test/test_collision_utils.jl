@testset "collision utils" begin
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, "Ar")

    interaction_data_path = joinpath(@__DIR__, "..", "data", "vhs.toml")
    interaction_data = load_interaction_data(interaction_data_path, species_data)

    @test length(species_data) == 1
    @test size(interaction_data) == (1,1)
    @test abs(interaction_data[1,1].μ1 - 0.5) < eps()
    @test abs(interaction_data[1,1].μ2 - 0.5) < eps()

    p1 = Particle(1e10, [2.0, 1.0, 0.0], [0.0, 0.0, 0.0])
    p2 = Particle(1e10, [0.0, -1.0, -1.0], [0.0, 0.0, 0.0])

    collision_data = CollisionData()

    Merzbild.compute_com!(collision_data, interaction_data[1,1], p1, p2)
    @test abs(collision_data.v_com[1] - 1.0) < eps()
    @test abs(collision_data.v_com[2] - 0.0) < eps()
    @test abs(collision_data.v_com[3] - (-0.5)) < eps()

    Merzbild.compute_g!(collision_data, p1, p2)
    @test abs(collision_data.g - 3.0) < eps()

    seed = 1234
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    # test that scattering doesn't change v_com and g magnitude
    Merzbild.scatter_vhs!(rng, collision_data, interaction_data[1,1], p1, p2)
    Merzbild.compute_com!(collision_data, interaction_data[1,1], p1, p2)

    @test abs(collision_data.v_com[1] - 1.0) < eps()
    @test abs(collision_data.v_com[2] - 0.0) < eps()
    @test abs(collision_data.v_com[3] - (-0.5)) < eps()

    Merzbild.compute_g!(collision_data, p1, p2)
    @test abs(collision_data.g - 3.0) < 2.1 * eps()

    species_data_2 = load_species_data(particles_data_path, ["Ar", "He"])
    interaction_data_2 = load_interaction_data(interaction_data_path, species_data_2)

    @test length(species_data_2) == 2
    @test size(interaction_data_2) == (2,2)
    @test abs(interaction_data_2[1,1].μ1 - 0.5) < eps()
    @test abs(interaction_data_2[1,1].μ2 - 0.5) < eps()
    @test abs(interaction_data_2[2,2].μ1 - 0.5) < eps()
    @test abs(interaction_data_2[2,2].μ2 - 0.5) < eps()
    @test abs(interaction_data_2[1,2].μ1 - interaction_data_2[2,1].μ2) < eps()
    @test abs(interaction_data_2[1,2].μ2 - interaction_data_2[2,1].μ1) < eps()
    @test abs(interaction_data_2[1,2].μ1 + interaction_data_2[1,2].μ2 - 1.0) < eps()



end