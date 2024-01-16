@testset "collision utils" begin
    species_list = load_species_list("data/particles.toml", "Ar")
    interaction_data = load_interaction_data("data/vhs.toml", species_list)

    @test size(interaction_data) == (1,1)
    @test abs(interaction_data[1,1].μ1 - 0.5) < eps()
    @test abs(interaction_data[1,1].μ2 - 0.5) < eps()

    p1 = Particle(1e10, [2.0, 1.0, 0.0], [0.0, 0.0, 0.0])
    p2 = Particle(1e10, [0.0, -1.0, -1.0], [0.0, 0.0, 0.0])

    collision_data = create_collision_data()

    Merzbild.compute_com!(collision_data, interaction_data[1,1], p1, p2)
    @test abs(collision_data.v_com[1] - 1.0) < eps()
    @test abs(collision_data.v_com[2] - 0.0) < eps()
    @test abs(collision_data.v_com[3] - (-0.5)) < eps()

    Merzbild.compute_g!(collision_data, p1, p2)
    @test abs(collision_data.g - 3.0) < eps()
# function compute_com(collision_data::CollisionData, interaction::Interaction, p1, p2)
#     collision_data.v_com = interaction.μ1 * p1.v + interaction.μ2 * p2.v
# end

# function compute_g(collision_data::CollisionData, p1, p2)
#     collision_data.g_vec = p1.v - p2.v
#     collision_data.g = norm(collision_data.g_vec)
# end
end