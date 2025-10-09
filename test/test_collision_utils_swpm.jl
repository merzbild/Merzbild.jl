@testset "collision utils for SWPM" begin
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    interaction_data_path = joinpath(@__DIR__, "..", "data", "vhs.toml")
    
    species_data2 = load_species_data(particles_data_path, ["Ar", "He"])
    interaction_data2 = load_interaction_data(interaction_data_path, species_data2)

    pia2 = ParticleIndexerArray(40, 2) # 40 cells, 2 species

    collision_factors_swpm = create_collision_factors_swpm_array(2)
    @test size(collision_factors_swpm) == (2,2,1)
    flag = false
    for i1 in 1:size(collision_factors_swpm)[1]
        for i2 in 1:size(collision_factors_swpm)[2]
            for i3 in 1:size(collision_factors_swpm)[3]
                if abs(collision_factors_swpm[i1, i2, i3].sigma_g_max) > 2 * eps()
                    flag = true
                end
            end
        end
    end
    @test flag == false

    collision_factors_swpm = create_collision_factors_swpm_array(2,40)
    @test size(collision_factors_swpm) == (2,2,40)
    flag = false
    for i1 in 1:size(collision_factors_swpm)[1]
        for i2 in 1:size(collision_factors_swpm)[2]
            for i3 in 1:size(collision_factors_swpm)[3]
                if abs(collision_factors_swpm[i1, i2, i3].sigma_g_max) > 2 * eps()
                    flag = true
                end
            end
        end
    end
    @test flag == false

    collision_factors_swpm = create_collision_factors_swpm_array(pia2)
    @test size(collision_factors_swpm) == (2,2,40)
    flag = false
    for i1 in 1:size(collision_factors_swpm)[1]
        for i2 in 1:size(collision_factors_swpm)[2]
            for i3 in 1:size(collision_factors_swpm)[3]
                if abs(collision_factors_swpm[i1, i2, i3].sigma_g_max) > 2 * eps()
                    flag = true
                end
            end
        end
    end
    @test flag == false

    # we use NTC version with Fnum = 1.0 as baseline
    collision_factors_ntc = create_collision_factors_array(2, 40)
    estimate_sigma_g_w_max!(collision_factors_ntc, interaction_data2, species_data2, [300.0, 300.0], 1.0; mult_factor=2.0)

    collision_factors_swpm = create_collision_factors_swpm_array(pia2, interaction_data2, species_data2, 300.0; mult_factor=2.0)
    @test size(collision_factors_swpm) == size(collision_factors_ntc)

    flag = false
    for i1 in 1:size(collision_factors_swpm)[1]
        for i2 in 1:size(collision_factors_swpm)[2]
            for i3 in 1:size(collision_factors_swpm)[3]
                if abs(collision_factors_swpm[i1, i2, i3].sigma_g_max - collision_factors_ntc[i1, i2, i3].sigma_g_w_max) > 2 * eps()
                    flag = true
                end
            end
        end
    end
    @test flag == false

    estimate_sigma_g_w_max!(collision_factors_ntc, interaction_data2, species_data2, [300.0, 600.0], 1; mult_factor=2.0)
    collision_factors_swpm = create_collision_factors_swpm_array(pia2, interaction_data2, species_data2, [300.0, 600.0]; mult_factor=2.0)
    @test size(collision_factors_swpm) == size(collision_factors_ntc)
    

    flag = false
    for i1 in 1:size(collision_factors_swpm)[1]
        for i2 in 1:size(collision_factors_swpm)[2]
            for i3 in 1:size(collision_factors_swpm)[3]
                if abs(collision_factors_swpm[i1, i2, i3].sigma_g_max - collision_factors_ntc[i1, i2, i3].sigma_g_w_max) > 2 * eps()
                    flag = true
                end
            end
        end
    end
    @test flag == false
end