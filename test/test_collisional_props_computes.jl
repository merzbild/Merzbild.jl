@testset "collisional prop computes" begin
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, "Ar")

    interaction_data_path = joinpath(@__DIR__, "..", "data", "vhs.toml")
    interaction_data = load_interaction_data(interaction_data_path, species_data)

    mfp1 = mean_free_path(interaction_data, 1, 1e23, 300.0)
    cf1 = mean_collision_frequency(interaction_data, species_data, 1, 1e23, 300.0)
    @test abs(mfp1 - 1.2940597144669493e-5) <= 2*eps()
    @test abs(cf1 - 2.907146431686186e7) <= 2*eps()

    # higher temperatures lead to smaller mfp and high coll frequency
    mfp2T = mean_free_path(interaction_data, 1, 1e23, 500.0)
    cf2T = mean_collision_frequency(interaction_data, species_data, 1, 1e23, 500.0)
    @test mfp2T < mfp1
    @test cf2T > cf1

    # same with density, and density scaling is linear
    mfp2n = mean_free_path(interaction_data, 1, 2e23, 300.0)
    cf2n = mean_collision_frequency(interaction_data, species_data, 1, 2e23, 300.0)
    @test mfp2n < mfp1
    @test cf2n > cf1
    @test abs(mfp2n/mfp1 - 0.5) <= 2*eps()
    @test abs(cf2n/cf1 - 2.0) <= 2*eps()

    # create pseudo-Maxwell molecule gas by hand, then MFP is independent of T
    interaction_data[1,1] = Interaction(species_data[1].mass/2, 0.5, 0.5, 4.11e-10, 0.5, 300.0, 0.0, 0.0)
    @test abs(mean_free_path(interaction_data, 1, 1e23, 300.0) - mean_free_path(interaction_data, 1, 1e23, 5000.0)) <= 2*eps()
end