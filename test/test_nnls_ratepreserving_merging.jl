@testset "nnls rate-preserving merging" begin

    function create_particles(weight_total)
        pv = ParticleVector(24)

        weight_per_particle = weight_total / 24

        # create some electron particles
        # E = 1/2 mv^2
        # 15.76 eV = 2.525E-18 J
        # v_ion = sqrt(2E/m)
        # m = 9.11e-31
        # v = 2.3545e6 m/s - any particle with v>this value will ionize
        
        ii = 1
        for octant in 1:8
            vx_s = Merzbild.vx_sign(octant)
            vy_s = Merzbild.vy_sign(octant)
            vz_s = Merzbild.vz_sign(octant)

            pv[ii] = Particle(weight_per_particle, [vx_s * 1e3, vy_s * 1.2e3, vz_s * 0.9e3], [0.0, 0.0, 0.0])
            Merzbild.update_particle_buffer_new_particle!(pv, ii)
            ii += 1

            pv[ii] = Particle(weight_per_particle, [vx_s * 3.2e5, vy_s * 0.9e4, vz_s * 4.5e5], [0.0, 0.0, 0.0])
            Merzbild.update_particle_buffer_new_particle!(pv, ii)
            ii += 1

            pv[ii] = Particle(weight_per_particle, [vx_s * 2.5e6, vy_s * 1.2e3, vz_s * 1.5e4], [0.0, 0.0, 0.0])
            Merzbild.update_particle_buffer_new_particle!(pv, ii)
            ii += 1
        end

        pia_ = ParticleIndexerArray(24)

        pia_.n_total[1] = 24
        pia_.indexer[1,1].n_local = 24
        pia_.indexer[1,1].n_group1 = 12
        pia_.indexer[1,1].start1 = 1
        pia_.indexer[1,1].end1 = 12
        pia_.indexer[1,1].n_group2 = 12
        pia_.indexer[1,1].start2 = 13
        pia_.indexer[1,1].end2 = 24

        return pv, pia_
    end

    seed = 1234
    Random.seed!(seed)
    rng = StableRNG(seed)

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, ["Ar", "e-"])

    interaction_data_path = joinpath(@__DIR__, "..", "data", "vhs.toml")
    interaction_data = load_interaction_data(interaction_data_path, species_data, 1e-10, 1.0, 273.0)

    e_n_data_path = joinpath(@__DIR__, "..", "data", "test_neutral_electron_data.xml")
    n_e_interactions = load_electron_neutral_interactions(species_data, e_n_data_path,
                                                          Dict("Ar" => "ConstantDB"),
                                                          Dict("Ar" => ScatteringIsotropic),
                                                          Dict("Ar" => ElectronEnergySplitEqual))

    computed_cs = create_computed_crosssections(n_e_interactions)

    # conserve only the lowest-order moments
    mim = [[0,0,0], [1,0,0], [0,1,0], [0,0,1], [2,0,0], [0,2,0], [0,0,2]]
    nnls = NNLSMerge(mim, 20; rate_preserving=false)
    nnls_rp = NNLSMerge(mim, 20; rate_preserving=true)

    # 7 moments only
    @test length(nnls.rhs_vector) == 7
    @test size(nnls.work[1].QA) == (7,20)
    @test length(nnls.work[1].Qb) == 7

    # 7 moments + 2 rates
    @test length(nnls_rp.rhs_vector) == 9
    @test size(nnls_rp.work[1].QA) == (9,20)
    @test length(nnls_rp.work[1].Qb) == 9

    ndens = 1e15
    particles, pia = create_particles(ndens)
    # set some reference values
    vref = 5e5
    cs_ref = 1e-19

    result = merge_nnls_based!(rng, nnls, particles, pia, 1, 1; vref=vref, scaling=:variance)

    @test result == 1
    @test pia.n_total[1] < 24

    # reset particles
    particles, pia = create_particles(ndens)

    w0 = sum([particles[i].w for i in 1:24])
    @test abs(w0 - ndens)/ndens < 2*eps()

    v0 = sum([particles[i].w * particles[i].v for i in 1:24])/w0
    @test maximum(abs.(v0)) < 2*eps()

    # mean velocity is 0, so we can just sum up 0.5 v^2 and get energy (per unit mass)
    E0 = 0.5 * sum([particles[i].w * sum(particles[i].v.^2) for i in 1:24])/w0
    
    rate_elastic = 0.0
    rate_ionization = 0.0

    for i in 1:24
        v_magnitude = norm(particles[i].v)
        
        Merzbild.compute_cross_sections!(computed_cs, interaction_data[1,2], v_magnitude,
                                         n_e_interactions, 1; extend=CSExtendConstant)

        cs_elastic = Merzbild.get_cs_elastic(n_e_interactions, computed_cs, 1)
        @test abs(cs_elastic - 1e-19) / 1e-19 < 2 * eps()

        rate_elastic += particles[i].w * cs_elastic * v_magnitude

        rate_ionization += particles[i].w * Merzbild.get_cs_ionization(n_e_interactions, computed_cs, 1) * v_magnitude
    end

    @test rate_ionization > 0.0

    result = merge_nnls_based_rate_preserving!(rng, nnls_rp,
                                               interaction_data, n_e_interactions, computed_cs,
                                               particles, pia, 1, 1, 1,
                                               cs_ref, cs_ref; scaling=:variance,
                                               vref=vref, n_rand_pairs=0, max_err=1e-11,
                                               centered_at_mean=false, v_multipliers=[], iteration_mult=2,
                                               extend=CSExtendConstant)

    @test result == 1
    np_new = pia.n_total[1]

    @test np_new < 24

    k_rate_elastic = rate_elastic / w0
    k_rate_ionization = rate_ionization / w0

    # test RHS rate coefficients
    @test abs(nnls_rp.rhs_vector[8] - k_rate_elastic)/k_rate_elastic < 4*eps()
    @test abs(nnls_rp.rhs_vector[9] - k_rate_ionization)/k_rate_ionization < 4*eps()

    @test maximum(abs.(nnls.rhs_vector - nnls_rp.rhs_vector[1:7])) < 4*eps()

    wnew = sum([particles[i].w for i in 1:np_new])
    @test abs(wnew - w0)/w0 < 1e-15

    v_new = sum([particles[i].w * particles[i].v for i in 1:np_new])/wnew
    @test maximum(abs.(v0)) < 1e-14

    # mean velocity is 0, so we can just sum up 0.5 v^2 and get energy (per unit mass)
    E_new = 0.5 * sum([particles[i].w * sum(particles[i].v.^2) for i in 1:np_new])/wnew
    @test abs(E_new - E0)/E0 < 1e-15

    # test conservation of rates
    rate_elastic_new = 0.0
    rate_ionization_new = 0.0

    for i in 1:np_new
        v_magnitude = norm(particles[i].v)
        
        Merzbild.compute_cross_sections!(computed_cs, interaction_data[1,2], v_magnitude,
                                         n_e_interactions, 1; extend=CSExtendConstant)

        cs_elastic = Merzbild.get_cs_elastic(n_e_interactions, computed_cs, 1)
        @test abs(cs_elastic - 1e-19) / 1e-19 < 2 * eps()

        rate_elastic_new += particles[i].w * cs_elastic * v_magnitude

        rate_ionization_new += particles[i].w * Merzbild.get_cs_ionization(n_e_interactions, computed_cs, 1) * v_magnitude
    end

    @test abs(rate_elastic_new - rate_elastic)/rate_elastic < 1e-15
    @test abs(rate_ionization_new - rate_ionization)/rate_ionization < 1e-15

    # finally test the edge case of init_np = total_np and no pre-allocated matrices, v_multipliers=[]
    # reset particles
    particles, pia = create_particles(ndens)

    nnls_rp = NNLSMerge(mim, 24; rate_preserving=true)

    result = merge_nnls_based_rate_preserving!(rng, nnls_rp,
                                               interaction_data, n_e_interactions, computed_cs,
                                               particles, pia, 1, 1, 1,
                                               cs_ref, cs_ref; scaling=:variance,
                                               vref=vref, n_rand_pairs=0, max_err=1e-11,
                                               centered_at_mean=false, v_multipliers=[], iteration_mult=2,
                                               extend=CSExtendConstant)

    @test result == 1
    @test pia.n_total[1] < 24
end