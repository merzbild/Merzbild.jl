@testset "malloc: NNLS merging internal computes" begin

    # particle setup with electrons (species 2) and neutrals (species 1),
    # split across both particle index groups so that the group2 branches are exercised
    function create_particles(weight_total; vel_scale=1.0, v_neutral_scale=1.0)
        pv = ParticleVector(24)

        weight_per_particle = weight_total / 24

        ii = 1
        for octant in 1:8
            vx_s = Merzbild.vx_sign(octant)
            vy_s = Merzbild.vy_sign(octant)
            vz_s = Merzbild.vz_sign(octant)

            pv[ii] = Particle(weight_per_particle, [vx_s * 1e3, vy_s * 1.2e3, vz_s * 0.9e3] .* vel_scale, [0.0, 0.0, 0.0])
            Merzbild.update_particle_buffer_new_particle!(pv, ii)
            ii += 1

            pv[ii] = Particle(weight_per_particle, [vx_s * 3.2e5, vy_s * 0.9e4, vz_s * 4.5e5] .* vel_scale, [0.0, 0.0, 0.0])
            Merzbild.update_particle_buffer_new_particle!(pv, ii)
            ii += 1

            pv[ii] = Particle(weight_per_particle, [vx_s * 2.5e6, vy_s * 1.2e3, vz_s * 1.5e4] .* vel_scale, [0.0, 0.0, 0.0])
            Merzbild.update_particle_buffer_new_particle!(pv, ii)
            ii += 1
        end

        pia_ = ParticleIndexerArray([24, 24])

        pia_.n_total[2] = 24
        pia_.indexer[1,2].n_local = 24
        pia_.indexer[1,2].n_group1 = 4
        pia_.indexer[1,2].start1 = 1
        pia_.indexer[1,2].end1 = 4
        pia_.indexer[1,2].n_group2 = 20
        pia_.indexer[1,2].start2 = 5
        pia_.indexer[1,2].end2 = 24

        pv_neutral = ParticleVector(24)
        weight_per_particle = weight_total * 3 / 24

        pia_.n_total[1] = 24
        pia_.indexer[1,1].n_local = 24
        pia_.indexer[1,1].n_group1 = 8
        pia_.indexer[1,1].start1 = 1
        pia_.indexer[1,1].end1 = 8
        pia_.indexer[1,1].n_group2 = 16
        pia_.indexer[1,1].start2 = 9
        pia_.indexer[1,1].end2 = 24

        ii = 1
        for octant in 1:8
            vx_s = Merzbild.vx_sign(octant)
            vy_s = Merzbild.vy_sign(octant)
            vz_s = Merzbild.vz_sign(octant)

            pv_neutral[ii] = Particle(weight_per_particle, [-vx_s * 1e2, -vy_s * 1.2e3, -vz_s * 0.9e2] .* v_neutral_scale, [0.0, 0.0, 0.0])
            Merzbild.update_particle_buffer_new_particle!(pv_neutral, ii)
            ii += 1

            pv_neutral[ii] = Particle(weight_per_particle, [-vx_s * 3.2e3, -vy_s * 0.9e3, -vz_s * 4.5e2] .* v_neutral_scale, [0.0, 0.0, 0.0])
            Merzbild.update_particle_buffer_new_particle!(pv_neutral, ii)
            ii += 1

            pv_neutral[ii] = Particle(weight_per_particle, [-vx_s * 2.5e2, -vy_s * 1.2e2, -vz_s * 1.5e2] .* v_neutral_scale, [0.0, 0.0, 0.0])
            Merzbild.update_particle_buffer_new_particle!(pv_neutral, ii)
            ii += 1
        end

        return [pv_neutral, pv], pia_
    end

    seed = 1234
    Random.seed!(seed)
    rng = StableRNG(seed)

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, ["Ar", "e-"])

    interaction_data_path = joinpath(@__DIR__, "..", "data", "vhs.toml")
    interaction_data = load_interaction_data(interaction_data_path, species_data, 1e-10, 1.0, 273.0)

    e_n_data_path = joinpath(@__DIR__, "..", "data", "test_neutral_electron_data.xml")
    n_e_interactions = ElectronNeutralInteractions(species_data, e_n_data_path,
                                                   Dict("Ar" => "ConstantDB"),
                                                   Dict("Ar" => ScatteringIsotropic),
                                                   Dict("Ar" => ElectronEnergySplitEqual))

    computed_cs = create_computed_crosssections(n_e_interactions)

    ndens = 1e15
    vref = 5e5
    cs_ref = 1e-19

    mim = [[0,0,0], [1,0,0], [0,1,0], [0,0,1], [2,0,0], [0,2,0], [0,0,2]]

    lhs_ncols = 24

    @testset "velocity-moment computes (non rate-preserving)" begin
        nnls = NNLSMerge(mim, 30)
        nnls.vref = vref
        nnls.inv_vref = 1.0 / vref

        particles, pia = create_particles(ndens)
        pv = particles[2]

        lhs_matrix = zeros(nnls.n_total_conserved, lhs_ncols)
        vel_pos_matrix = zeros(6, lhs_ncols)

        # warmup (compile) + verify the computes actually run
        Merzbild.compute_w_total_v0!(nnls, pv, pia, 1, 2)
        Merzbild.compute_lhs_and_rhs!(nnls, lhs_matrix, vel_pos_matrix, pv, pia, 1, 2)
        Merzbild.scale_lhs_rhs!(nnls, lhs_matrix, :vref, lhs_ncols)
        Merzbild.compute_lhs_and_rhs!(nnls, lhs_matrix, vel_pos_matrix, pv, pia, 1, 2)
        Merzbild.scale_lhs_rhs!(nnls, lhs_matrix, :variance, lhs_ncols)

        bytes_wv0 = @allocated Merzbild.compute_w_total_v0!(nnls, pv, pia, 1, 2)
        @test bytes_wv0 == 0

        bytes_lhs = @allocated Merzbild.compute_lhs_and_rhs!(nnls, lhs_matrix, vel_pos_matrix, pv, pia, 1, 2)
        @test bytes_lhs == 0

        bytes_scale_vref = @allocated Merzbild.scale_lhs_rhs!(nnls, lhs_matrix, :vref, lhs_ncols)
        @test bytes_scale_vref == 0

        Merzbild.compute_lhs_and_rhs!(nnls, lhs_matrix, vel_pos_matrix, pv, pia, 1, 2)
        bytes_scale_variance = @allocated Merzbild.scale_lhs_rhs!(nnls, lhs_matrix, :variance, lhs_ncols)
        @test bytes_scale_variance == 0
    end

    @testset "approximate rate-preserving computes" begin
        nnls_rp = NNLSMerge(mim, 30; rate_preserving=true)
        nnls_rp.vref = vref
        nnls_rp.inv_vref = 1.0 / vref

        particles, pia = create_particles(ndens)
        pv = particles[2]

        lhs_matrix = zeros(nnls_rp.n_total_conserved, lhs_ncols)
        vel_pos_matrix = zeros(6, lhs_ncols)

        # warmup (compile)
        Merzbild.compute_lhs_and_rhs_rate_preserving!(nnls_rp, lhs_matrix, vel_pos_matrix,
                                                      interaction_data[1,2], n_e_interactions, computed_cs,
                                                      pv, pia, 1, 2, 1, CSExtendConstant)
        Merzbild.scale_lhs_rhs_rate_preserving!(nnls_rp, lhs_matrix, cs_ref, cs_ref, :variance, lhs_ncols)

        bytes_lhs = @allocated Merzbild.compute_lhs_and_rhs_rate_preserving!(nnls_rp, lhs_matrix, vel_pos_matrix,
                                                                            interaction_data[1,2], n_e_interactions, computed_cs,
                                                                            pv, pia, 1, 2, 1, CSExtendConstant)
        @test bytes_lhs == 0

        bytes_scale = @allocated Merzbild.scale_lhs_rhs_rate_preserving!(nnls_rp, lhs_matrix, cs_ref, cs_ref, :variance, lhs_ncols)
        @test bytes_scale == 0
    end

    @testset "exact rate-preserving computes" begin
        nnls_rp = NNLSMerge(mim, 30; rate_preserving=true)
        nnls_rp.vref = vref
        nnls_rp.inv_vref = 1.0 / vref

        particles, pia = create_particles(ndens; v_neutral_scale=2.0)
        pv = particles[2]
        pv_neutral = particles[1]

        lhs_matrix = zeros(nnls_rp.n_total_conserved, lhs_ncols)
        vel_pos_matrix = zeros(6, lhs_ncols)

        # warmup (compile)
        Merzbild.compute_lhs_and_rhs_rate_preserving!(nnls_rp, lhs_matrix, vel_pos_matrix,
                                                      interaction_data[1,2], n_e_interactions, computed_cs,
                                                      pv, pv_neutral, pia, 1, 2, 1, CSExtendConstant)

        bytes_lhs = @allocated Merzbild.compute_lhs_and_rhs_rate_preserving!(nnls_rp, lhs_matrix, vel_pos_matrix,
                                                                            interaction_data[1,2], n_e_interactions, computed_cs,
                                                                            pv, pv_neutral, pia, 1, 2, 1, CSExtendConstant)
        @test bytes_lhs == 0
    end
end
