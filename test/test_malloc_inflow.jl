
@testset "malloc tests: inflow velocity sampling" begin

    T = 300.0
    ndens = 1e20
    v_wall = 500.0

    seed = 1234
    rng = StableRNG(seed)

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")

    normal_left = SVector{3,Float64}(1.0, 0.0, 0.0)
    normal_right = SVector{3,Float64}(-1.0, 0.0, 0.0)

    bc_no_drift = InflowBC1D(1, species_data, T, SVector{3,Float64}(0.0, v_wall, 0.0), normal_left, ndens)
    bc_envelope4 = InflowBC1D(1, species_data, T, SVector{3,Float64}(v_wall, 0.0, 0.0), normal_left, ndens)
    bc_envelope4_right = InflowBC1D(1, species_data, T, SVector{3,Float64}(-v_wall, 0.0, 0.0), normal_right, ndens)

    @test bc_no_drift.sampler === InflowNoDrift()
    @test bc_envelope4.sampler === InflowEnvelope4()

    # sampling many particles in a loop, as an emission routine would
    function sample_many(rng, inflowbc, n)
        v = SVector{3,Float64}(0.0, 0.0, 0.0)
        for _ in 1:n
            v += sample_inflow_velocity(rng, inflowbc)
        end
        return v
    end

    for bc in (bc_no_drift, bc_envelope4, bc_envelope4_right)
        sample_many(rng, bc, 100)

        # the acceptance-rejection loop must not allocate over many samples either
        bytes = @allocated sample_many(rng, bc, 10000)
        @test bytes == 0
    end

    compute_inflow_number_flux(1, species_data, T, SVector{3,Float64}(v_wall, 0.0, 0.0), normal_left, ndens)
    bytes = @allocated compute_inflow_number_flux(1, species_data, T, SVector{3,Float64}(v_wall, 0.0, 0.0),
                                                  normal_left, ndens)
    @test bytes == 0
end
