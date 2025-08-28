@testset "collision_fp" begin
    #scale_norm_rands!
    using Statistics

    seed = 1234
    rng = StableRNG(seed)
    n_particles = 100
    mean_v = zeros(3)
    stddev = zeros(3)

    collision_data_fp = CollisionDataFP()

    resize!(collision_data_fp.xvel_rand, n_particles)
    resize!(collision_data_fp.yvel_rand, n_particles)
    resize!(collision_data_fp.zvel_rand, n_particles)

    # norm_rands = randn(rng, Float64, (n_particles, 3))
    Merzbild.sample_normal_rands!(rng, collision_data_fp, n_particles)
    Merzbild.scale_norm_rands!(collision_data_fp, n_particles)

    mean_v = [mean(collision_data_fp.xvel_rand),
              mean(collision_data_fp.yvel_rand),
              mean(collision_data_fp.zvel_rand)]
    
    stddev = [std(collision_data_fp.xvel_rand; corrected=false),
              std(collision_data_fp.yvel_rand; corrected=false),
              std(collision_data_fp.zvel_rand; corrected=false)]

    @test isapprox(mean_v, zeros(3); atol=1e-15)
    @test isapprox(stddev, ones(3); atol=1e-14)

    # test that scaling accounts only for the first n_local particles!
    collision_data_fp.xvel_rand[:] .= 10000.0
    collision_data_fp.yvel_rand[:] .= 10000.0
    collision_data_fp.zvel_rand[:] .= 10000.0
    
    collision_data_fp.xvel_rand[1:4] .= [1.0, -3.0, 1.0, 2.0]
    collision_data_fp.yvel_rand[1:4] .= [2.5, -2.5, 0.5, 1.0]
    collision_data_fp.zvel_rand[1:4] .= [4.0, 3.0, 1.0, 2.0]

    Merzbild.scale_norm_rands!(collision_data_fp, 4)
    @test isapprox(sum(collision_data_fp.xvel_rand[1:4])/4.0, 0.0; atol=1e-15)
    @test isapprox(sum(collision_data_fp.yvel_rand[1:4])/4.0, 0.0; atol=1e-15)
    @test isapprox(sum(collision_data_fp.zvel_rand[1:4])/4.0, 0.0; atol=1e-15)

    @test isapprox(std(collision_data_fp.xvel_rand[1:4]; corrected=false), 1.0; atol=1e-15)
    @test isapprox(std(collision_data_fp.yvel_rand[1:4]; corrected=false), 1.0; atol=1e-15)
    @test isapprox(std(collision_data_fp.zvel_rand[1:4]; corrected=false), 1.0; atol=1e-15)
end