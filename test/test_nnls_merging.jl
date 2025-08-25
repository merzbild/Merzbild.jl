@testset "nnls_merging" begin

    function create_particles()
        pv = ParticleVector(50)

        for i in 1:50
            Merzbild.update_particle_buffer_new_particle!(pv, i)
            pv.particles[i] = Particle(1.0, [i, 5 - 2.0 * i, -100.0 + 0.25 * i^2], [-0.5 * i, 27 + i, 400.0 - 0.5 * i^2])
        end

        return pv
    end

    function create_particles2()
        pv = ParticleVector(16)

        pv.particles = [Particle(1.0, [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]),
                        Particle(1.0, [1.0, 1.0, -1.0], [0.0, 0.0, 0.0]),
                        Particle(1.0, [1.0, -1.0, 1.0], [0.0, 0.0, 0.0]),
                        Particle(1.0, [1.0, -1.0, -1.0], [0.0, 0.0, 0.0]),
                        Particle(1.0, [-1.0, 1.0, 1.0], [0.0, 0.0, 0.0]),
                        Particle(1.0, [-1.0, 1.0, -1.0], [0.0, 0.0, 0.0]),
                        Particle(1.0, [-1.0, -1.0, 1.0], [0.0, 0.0, 0.0]),
                        Particle(1.0, [-1.0, -1.0, -1.0], [0.0, 0.0, 0.0]),
                        Particle(1.0, [2.0, 2.0, 2.0], [0.0, 0.0, 0.0]),
                        Particle(1.0, [2.0, 2.0, -2.0], [0.0, 0.0, 0.0]),
                        Particle(1.0, [2.0, -2.0, 2.0], [0.0, 0.0, 0.0]),
                        Particle(1.0, [2.0, -2.0, -2.0], [0.0, 0.0, 0.0]),
                        Particle(1.0, [-2.0, 2.0, 2.0], [0.0, 0.0, 0.0]),
                        Particle(1.0, [-2.0, 2.0, -2.0], [0.0, 0.0, 0.0]),
                        Particle(1.0, [-2.0, -2.0, 2.0], [0.0, 0.0, 0.0]),
                        Particle(1.0, [-2.0, -2.0, -2.0], [0.0, 0.0, 0.0])]

        for i in 1:16
            Merzbild.update_particle_buffer_new_particle!(pv, i)
        end
        
        return pv
    end



    function create_particles3()
        pv = ParticleVector(6)

        pv.particles = [Particle(1.0, [3.0, 1.0, 1.0], [0.0, 0.0, 0.0]),
                        Particle(1.0, [1.0, 4.0, -1.0], [0.0, 0.0, 0.0]),
                        Particle(1.0, [1.0, -8.0, 1.0], [0.0, 0.0, 0.0]),
                        Particle(1.0, [1.0, -1.0, -9.0], [0.0, 0.0, 0.0]),
                        Particle(1.0, [-7.0, -11.0, 1.0], [0.0, 0.0, 0.0]),
                        Particle(1.0, [-9.0, 6.0, -1.0], [0.0, 0.0, 0.0])]

        for i in 1:6
            Merzbild.update_particle_buffer_new_particle!(pv, i)
        end
        
        return pv
    end


    function create_particles4()
        pv = ParticleVector(2)

        pv.particles = [Particle(1.0, [3.0, 1.0, -4.0], [0.0, 0.0, 0.0]),
                        Particle(1.0, [-3.0, 1.0, 4.0], [0.0, 0.0, 0.0]),]

        for i in 1:2
            Merzbild.update_particle_buffer_new_particle!(pv, i)
        end
        
        return pv
    end

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, "Ar")

    seed = 1234
    Random.seed!(seed)
    rng = StableRNG(seed)

    Nx = 2
    Ny = 2
    Nz = 2

    phys_props::PhysProps = PhysProps(1, 1, [4], Tref=1)

    Δabs = 2.5
    Δrel_xsmall = 5e-13
    
    particles = [create_particles()]
    pia = ParticleIndexerArray(length(particles[1]))

    compute_props!(particles, pia, species_data, phys_props)

    mim = []
    n_moms = 4
    for i in 1:n_moms
        append!(mim, compute_multi_index_moments(i))
    end

    mnnls = NNLSMerge(mim, 30)
    vref = sqrt(2 * k_B * 300.0 / species_data[1].mass)
    result = merge_nnls_based!(rng, mnnls, particles[1], pia, 1, 1; vref=vref, scaling=:vref)
    @test length(mnnls.rhs_vector) == length(mim) + 1  # we didn't construct 0-th order moment

    n0 = phys_props.n[1, 1]
    np0 = phys_props.np[1, 1]
    v0 = phys_props.v[:, 1, 1]
    T0 = phys_props.T[1, 1]
    M40 = phys_props.moments[1, 1, 1]

    mixed_order_3 = compute_multi_index_moments(3)
    moms3 = zeros(length(mixed_order_3))

    for (i, mim) in enumerate(mixed_order_3)
        moms3[i] = Merzbild.compute_mixed_moment(particles, pia, 1, 1, mim)
    end

    # test some internal computes used in the merging
    @test sum(abs.(mnnls.v0 .- v0)) < eps()
    @test sum(abs.(mnnls.w_total .- n0)) < eps()

    @test abs(mnnls.minvx + v0[1] - 1.0) < eps()
    @test abs(mnnls.maxvx + v0[1] - 50.0) < eps()
    @test abs(mnnls.minvy + v0[2] - (5.0 - 2.0 * 50)) < eps()
    @test abs(mnnls.maxvy + v0[2] - (5.0 - 2.0 * 1)) < eps()
    @test abs(mnnls.minvz + v0[3] - (-100 + 0.25)) < eps()
    @test abs(mnnls.maxvz + v0[3] - (-100 + 0.25 * 2500)) < eps()

    compute_props!(particles, pia, species_data, phys_props)
    # test that merging conserves mass / momentum / energy
    @test phys_props.np[1, 1] < np0
    @test abs(n0 - phys_props.n[1, 1]) < 1.5e-14
    @test sum(abs.(v0 - phys_props.v[:, 1, 1])) < 2e-13
    @test abs(T0 - phys_props.T[1, 1]) < 3.6e-14

    # test higher-order moment in g
    @test abs(M40 - phys_props.moments[1, 1, 1]) < 7.4e-12

    # test moments of order 3
    moms3_post = zeros(length(mixed_order_3))
    
    for (i, mim) in enumerate(mixed_order_3)
        moms3_post[i] = Merzbild.compute_mixed_moment(particles, pia, 1, 1, mim)
    end

    @test sum(abs.(moms3_post - moms3))/length(moms3_post) < 1e-14


    mim = []
    n_moms = 2
    for i in 1:n_moms
        append!(mim, compute_multi_index_moments(i))
    end

    mnnls2 = NNLSMerge(mim, 30)
    particles2 = [create_particles2()]
    pia2 = ParticleIndexerArray(length(particles2[1]))
    vref = 1.0

    result = merge_nnls_based!(rng, mnnls2, particles2[1], pia2, 1, 1; vref=vref, scaling=:vref)
    @test abs(mnnls2.rhs_vector[1] - 1.0) < eps()

    # mean velocity is 0.0
    @test abs(mnnls2.rhs_vector[2] - 0.0) < eps()
    @test abs(mnnls2.rhs_vector[3] - 0.0) < eps()
    @test abs(mnnls2.rhs_vector[4] - 0.0) < eps()

    # mixed second-order moments are also 0.0
    @test abs(mnnls2.rhs_vector[8] - 0.0) < eps()
    @test abs(mnnls2.rhs_vector[9] - 0.0) < eps()
    @test abs(mnnls2.rhs_vector[10] - 0.0) < eps()


    # check we always preserve the lowest-order moments
    mnnls = NNLSMerge([], 30)
    @test [0,0,0] in mnnls.mim
    @test [1,0,0] in mnnls.mim
    @test [0,1,0] in mnnls.mim
    @test [0,0,1] in mnnls.mim
    @test [2,0,0] in mnnls.mim
    @test [0,2,0] in mnnls.mim
    @test [0,0,2] in mnnls.mim
    @test mnnls.n_moments == 7

    # check computation of LHS/RHS of the system

    particles = [create_particles3()]
    pia = ParticleIndexerArray(length(particles[1]))

    centered_at_mean = true
    v_multipliers = [0.25, 0.5, 1.0]
    n_add = centered_at_mean ? 1 : 0
    n_rand_pairs = 0
    n_add += 8 * length(v_multipliers)
    lhs_ncols = pia.indexer[1, 1].n_local + n_add + n_rand_pairs
    lhs_matrix = zeros(mnnls.n_moments, lhs_ncols)

    Merzbild.compute_lhs_and_rhs!(mnnls, lhs_matrix,
                                  particles[1], pia, 1, 1)
    # we didn't compute anything using the additional particles
    @test sum(abs.(lhs_matrix[:, pia.indexer[1, 1].n_local+1:end])) == 0

    Merzbild.compute_lhs_particles_additional!(rng, pia.indexer[1, 1].n_local+1, mnnls, lhs_matrix,
                                               particles[1], pia, 1, 1,
                                               n_rand_pairs, centered_at_mean, v_multipliers)

    # centered particle adds zero values
    @test sum(abs.(lhs_matrix[:, pia.indexer[1, 1].n_local+1])) == 0

    # the other particles add non-zero values
    @test sum(abs.(lhs_matrix[:, pia.indexer[1, 1].n_local+2:end])) != 0

    centered_at_mean = false
    v_multipliers = [0.25, 0.5, 1.0]
    n_add = centered_at_mean ? 1 : 0
    n_rand_pairs = 0
    n_add += 8 * length(v_multipliers)
    lhs_ncols = pia.indexer[1, 1].n_local + n_add + n_rand_pairs
    lhs_matrix = zeros(mnnls.n_moments, lhs_ncols)

    Merzbild.compute_lhs_and_rhs!(mnnls, lhs_matrix,
                                  particles[1], pia, 1, 1)
    # we didn't compute anything using the additional particles
    @test sum(abs.(lhs_matrix[:, pia.indexer[1, 1].n_local+1:end])) == 0

    Merzbild.compute_lhs_particles_additional!(rng, pia.indexer[1, 1].n_local+1, mnnls, lhs_matrix,
                                               particles[1], pia, 1, 1,
                                               n_rand_pairs, centered_at_mean, v_multipliers)

    # no centered particle, so no zero values
    @test sum(abs.(lhs_matrix[:, pia.indexer[1, 1].n_local+1:end])) != 0

    mim = []
    #  [[0, 0, 0],
    # [1, 0, 0], [0, 1, 0], [0, 0, 1],
    # [2, 0, 0], [0, 2, 0], [0, 0, 2]]
    mnnls3 = NNLSMerge(mim, 30)

    lhs_ncols = 2
    # 7, 
    lhs_matrix = zeros(mnnls3.n_moments, lhs_ncols)

    particles = [create_particles4()]
    pia = ParticleIndexerArray(length(particles[1]))

    Merzbild.compute_lhs_and_rhs!(mnnls3, lhs_matrix,
                                  particles[1], pia, 1, 1)

    lhs_matrix_2 = copy(lhs_matrix)
    rhs_vec_2 = copy(mnnls3.rhs_vector)
    mnnls3.vref = 2.0
    mnnls3.inv_vref = 1.0 / mnnls3.vref
    Merzbild.scale_lhs_rhs!(mnnls3, lhs_matrix, :vref)
    
    @test abs(rhs_vec_2[1] - mnnls3.rhs_vector[1]) < 2*eps()
    @test abs(rhs_vec_2[2] / mnnls3.vref - mnnls3.rhs_vector[2]) < 2*eps()
    @test abs(rhs_vec_2[3] / mnnls3.vref - mnnls3.rhs_vector[3]) < 2*eps()
    @test abs(rhs_vec_2[4] / mnnls3.vref - mnnls3.rhs_vector[4]) < 2*eps()
    @test abs(rhs_vec_2[5] / mnnls3.vref^2 - mnnls3.rhs_vector[5]) < 2*eps()
    @test abs(rhs_vec_2[6] / mnnls3.vref^2 - mnnls3.rhs_vector[6]) < 2*eps()
    @test abs(rhs_vec_2[7] / mnnls3.vref^2 - mnnls3.rhs_vector[7]) < 2*eps()

    for i in 1:2
        @test maximum(abs.(lhs_matrix_2[1,i] - lhs_matrix[1,i])) < 2*eps()
        @test maximum(abs.(lhs_matrix_2[2,i] / mnnls3.vref - lhs_matrix[2,i])) < 2*eps()
        @test maximum(abs.(lhs_matrix_2[3,i] / mnnls3.vref - lhs_matrix[3,i])) < 2*eps()
        @test maximum(abs.(lhs_matrix_2[4,i] / mnnls3.vref - lhs_matrix[4,i])) < 2*eps()
        @test maximum(abs.(lhs_matrix_2[5,i] / mnnls3.vref^2 - lhs_matrix[5,i])) < 2*eps()
        @test maximum(abs.(lhs_matrix_2[6,i] / mnnls3.vref^2 - lhs_matrix[6,i])) < 2*eps()
        @test maximum(abs.(lhs_matrix_2[7,i] / mnnls3.vref^2 - lhs_matrix[7,i])) < 2*eps()
    end

    mnnls3.rhs_vector = copy(rhs_vec_2)
    lhs_matrix = copy(lhs_matrix_2)

    mnnls3.Ex = 2.0
    mnnls3.Ey = 3.0
    mnnls3.Ez = 4.0

    Merzbild.scale_lhs_rhs!(mnnls3, lhs_matrix, :variance)

    @test abs(rhs_vec_2[1] - mnnls3.rhs_vector[1]) < 2*eps()
    @test abs(rhs_vec_2[2] / 2.0 - mnnls3.rhs_vector[2]) < 2*eps()
    @test abs(rhs_vec_2[3] / 3.0 - mnnls3.rhs_vector[3]) < 2*eps()
    @test abs(rhs_vec_2[4] / 4.0 - mnnls3.rhs_vector[4]) < 2*eps()
    @test abs(rhs_vec_2[5] / 4.0 - mnnls3.rhs_vector[5]) < 2*eps()
    @test abs(rhs_vec_2[6] / 9.0 - mnnls3.rhs_vector[6]) < 2*eps()
    @test abs(rhs_vec_2[7] / 16.0 - mnnls3.rhs_vector[7]) < 2*eps()

    for i in 1:2
        @test maximum(abs.(lhs_matrix_2[1,i] - lhs_matrix[1,i])) < 2*eps()
        @test maximum(abs.(lhs_matrix_2[2,i] / 2.0 - lhs_matrix[2,i])) < 2*eps()
        @test maximum(abs.(lhs_matrix_2[3,i] / 3.0 - lhs_matrix[3,i])) < 2*eps()
        @test maximum(abs.(lhs_matrix_2[4,i] / 4.0 - lhs_matrix[4,i])) < 2*eps()
        @test maximum(abs.(lhs_matrix_2[5,i] / 4.0 - lhs_matrix[5,i])) < 2*eps()
        @test maximum(abs.(lhs_matrix_2[6,i] / 9.0 - lhs_matrix[6,i])) < 2*eps()
        @test maximum(abs.(lhs_matrix_2[7,i] / 16.0 - lhs_matrix[7,i])) < 2*eps()
    end
end