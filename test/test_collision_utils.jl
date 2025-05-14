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


    sp_data, int_data = load_species_and_interaction_data(particles_data_path, interaction_data_path,
                                                          ["Ar", "He"]; fill_dummy=false)
    
    @test length(sp_data) == 2
    for i in 1:2
        for j in 1:2
            @test abs(interaction_data_2[i,j].μ1 - int_data[i,j].μ1) < eps()
            @test abs(interaction_data_2[i,j].μ2 - int_data[i,j].μ2) < eps()
            @test abs(interaction_data_2[i,j].m_r - int_data[i,j].m_r) < eps()
            @test abs(interaction_data_2[i,j].vhs_d - int_data[i,j].vhs_d) < eps()
            @test abs(interaction_data_2[i,j].vhs_o - int_data[i,j].vhs_o) < eps()
            @test abs(interaction_data_2[i,j].vhs_Tref - int_data[i,j].vhs_Tref) < eps()
            @test abs(interaction_data_2[i,j].vhs_muref - int_data[i,j].vhs_muref) < eps()
            @test abs(interaction_data_2[i,j].vhs_factor - int_data[i,j].vhs_factor) < eps()
        end
    end

    e_sp_data = load_species_data(particles_data_path, "e-")
    @test_throws KeyError e_int_data = load_interaction_data(interaction_data_path, e_sp_data)

    e_int_data = load_interaction_data_with_dummy(interaction_data_path, e_sp_data)
    @test abs(e_int_data[1,1].vhs_d - 1e-10) < eps()
    @test abs(e_int_data[1,1].vhs_o - 1.0) < eps()
    @test abs(e_int_data[1,1].vhs_Tref - 273.0) < eps()

    e_sp_data2, e_int_data2 = load_species_and_interaction_data(particles_data_path, interaction_data_path,
                                                                ["e-"]; fill_dummy=true)

    @test abs(e_int_data2[1,1].μ1 - e_int_data[1,1].μ1) < eps()
    @test abs(e_int_data2[1,1].μ2 - e_int_data[1,1].μ2) < eps()
    @test abs(e_int_data2[1,1].m_r - e_int_data[1,1].m_r) < eps()
    @test abs(e_int_data2[1,1].vhs_d - e_int_data[1,1].vhs_d) < eps()
    @test abs(e_int_data2[1,1].vhs_o - e_int_data[1,1].vhs_o) < eps()
    @test abs(e_int_data2[1,1].vhs_Tref - e_int_data[1,1].vhs_Tref) < eps()
    @test abs(e_int_data2[1,1].vhs_muref - e_int_data[1,1].vhs_muref) < eps()
    @test abs(e_int_data2[1,1].vhs_factor - e_int_data[1,1].vhs_factor) < eps()

    species_data2 = load_species_data(particles_data_path, ["Ar", "He"])
    interaction_data2 = load_interaction_data(interaction_data_path, species_data2)
    pia2 = ParticleIndexerArray(40, 2) # 40 cells, 2 species

    collision_factors1 = create_collision_factors_array(2, 40)
    collision_factors2 = create_collision_factors_array(pia2)
    @test size(collision_factors1) == size(collision_factors2)

    # we will use as the baseline
    estimate_sigma_g_w_max!(collision_factors1, interaction_data2, species_data2, [300.0, 300.0], 1e10; mult_factor=2.0)
    collision_factors3 = create_collision_factors_array(pia2, interaction_data2, species_data2, 300.0, 1e10; mult_factor=2.0)
    @test size(collision_factors3) == size(collision_factors1)

    flag = false
    for i1 in 1:size(collision_factors3)[1]
        for i2 in 1:size(collision_factors3)[2]
            for i3 in 1:size(collision_factors3)[3]
                if abs(collision_factors3[i1, i2, i3].sigma_g_w_max - collision_factors1[i1, i2, i3].sigma_g_w_max) > 2 * eps()
                    flag = true
                end
            end
        end
    end
    @test flag == false

    estimate_sigma_g_w_max!(collision_factors1, interaction_data2, species_data2, [300.0, 600.0], 1e10; mult_factor=2.0)
    collision_factors4 = create_collision_factors_array(pia2, interaction_data2, species_data2, [300.0, 600.0],
                                                        1e10; mult_factor=2.0)
    @test size(collision_factors4) == size(collision_factors1)
    

    flag = false
    for i1 in 1:size(collision_factors4)[1]
        for i2 in 1:size(collision_factors4)[2]
            for i3 in 1:size(collision_factors4)[3]
                if abs(collision_factors4[i1, i2, i3].sigma_g_w_max - collision_factors1[i1, i2, i3].sigma_g_w_max) > 2 * eps()
                    flag = true
                end
            end
        end
    end
    @test flag == false
    # collision_data = CollisionData()
end