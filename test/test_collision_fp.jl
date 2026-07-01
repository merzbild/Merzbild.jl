@testset "collision_fp" begin
    #scale_norm_rands!
    using Statistics

    seed = 1234
    rng = StableRNG(seed)
    n_particles = 100

    collision_data_fp = CollisionDataFP()

    resize!(collision_data_fp.xvel_rand, n_particles)
    resize!(collision_data_fp.yvel_rand, n_particles)
    resize!(collision_data_fp.zvel_rand, n_particles)

    # equal-weight particles: the weighted scaling reduces to the plain (unweighted) case
    particles = ParticleVector(n_particles)
    for i in 1:n_particles
        Merzbild.add_particle!(particles, i, 1.0,
                               SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0))
    end
    indexer = ParticleIndexer(n_particles)

    # norm_rands = randn(rng, Float64, (n_particles, 3))
    Merzbild.sample_normal_rands!(rng, collision_data_fp, n_particles)
    Merzbild.scale_norm_rands!(collision_data_fp, particles, indexer, 1.0 * n_particles)

    mean_v = [mean(collision_data_fp.xvel_rand),
              mean(collision_data_fp.yvel_rand),
              mean(collision_data_fp.zvel_rand)]

    stddev = [std(collision_data_fp.xvel_rand; corrected=false),
              std(collision_data_fp.yvel_rand; corrected=false),
              std(collision_data_fp.zvel_rand; corrected=false)]

    @test isapprox(mean_v, zeros(3); atol=1e-15)
    @test isapprox(stddev, ones(3); atol=1e-14)

    # test that scaling accounts only for the particles in the indexer!
    collision_data_fp.xvel_rand[:] .= 10000.0
    collision_data_fp.yvel_rand[:] .= 10000.0
    collision_data_fp.zvel_rand[:] .= 10000.0

    collision_data_fp.xvel_rand[1:4] .= [1.0, -3.0, 1.0, 2.0]
    collision_data_fp.yvel_rand[1:4] .= [2.5, -2.5, 0.5, 1.0]
    collision_data_fp.zvel_rand[1:4] .= [4.0, 3.0, 1.0, 2.0]

    Merzbild.scale_norm_rands!(collision_data_fp, particles, ParticleIndexer(4), 4.0)
    @test isapprox(sum(collision_data_fp.xvel_rand[1:4])/4.0, 0.0; atol=1e-15)
    @test isapprox(sum(collision_data_fp.yvel_rand[1:4])/4.0, 0.0; atol=1e-15)
    @test isapprox(sum(collision_data_fp.zvel_rand[1:4])/4.0, 0.0; atol=1e-15)

    @test isapprox(std(collision_data_fp.xvel_rand[1:4]; corrected=false), 1.0; atol=1e-15)
    @test isapprox(std(collision_data_fp.yvel_rand[1:4]; corrected=false), 1.0; atol=1e-15)
    @test isapprox(std(collision_data_fp.zvel_rand[1:4]; corrected=false), 1.0; atol=1e-15)

    # variable-weight particles: the weight-averaged mean must be 0 and the weight-averaged variance 1
    weights = [1.0, 2.0, 3.0, 4.0]
    local_w = sum(weights)
    particles_w = ParticleVector(4)
    for i in 1:4
        Merzbild.add_particle!(particles_w, i, weights[i],
                               SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0))
    end

    collision_data_fp.xvel_rand[1:4] .= [1.0, -3.0, 1.0, 2.0]
    collision_data_fp.yvel_rand[1:4] .= [2.5, -2.5, 0.5, 1.0]
    collision_data_fp.zvel_rand[1:4] .= [4.0, 3.0, 1.0, 2.0]

    Merzbild.scale_norm_rands!(collision_data_fp, particles_w, ParticleIndexer(4), local_w)

    for rand_comp in (collision_data_fp.xvel_rand, collision_data_fp.yvel_rand, collision_data_fp.zvel_rand)
        wmean = sum(weights[i] * rand_comp[i] for i in 1:4) / local_w
        wvar = sum(weights[i] * rand_comp[i]^2 for i in 1:4) / local_w
        @test isapprox(wmean, 0.0; atol=1e-15)
        @test isapprox(wvar, 1.0; atol=1e-14)
    end
end
