@testset "collision_fp" begin
    #scale_norm_rands!
    seed = 1234
    rng = StableRNG(seed)
    n_particles = 100
    mean = zeros(3)
    stddev = zeros(3)

    collision_data_fp = CollisionDataFP()
    norm_rands = randn(rng, Float64, (n_particles, 3))
    Merzbild.scale_norm_rands!(collision_data_fp, norm_rands)

    for i in 1:n_particles
        rand = norm_rands[i, :]
        mean += rand / n_particles
    end

    for i in 1:n_particles
        stddev[1] += norm_rands[i, 1] * norm_rands[i, 1] / n_particles
        stddev[2] += norm_rands[i, 2] * norm_rands[i, 2] / n_particles
        stddev[3] += norm_rands[i, 3] * norm_rands[i, 3] / n_particles
    end

    @test isapprox(mean, zeros(3); atol=1e-15)
    @test isapprox(stddev, ones(3); atol=1e-14)
end