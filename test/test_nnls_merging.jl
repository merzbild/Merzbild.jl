@testset "nnls_merging" begin

    function create_particles()
        vp::Vector{Particle} = []

        for i in 1:50
            push!(vp, Particle(1.0, [i, 5 - 2.0 * i, -100.0 + 0.25 * i^2], [-0.5 * i, 27 + i, 400.0 - 0.5 * i^2]))
        end

        return vp
    end

    function create_particles2()
        vp::Vector{Particle} = [Particle(1.0, [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]),
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

        
        return vp
    end

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, "Ar")

    seed = 1234
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    Nx = 2
    Ny = 2
    Nz = 2

    phys_props::PhysProps = PhysProps(1, 1, [4], Tref=1)

    Δabs = 2.5
    Δrel_xsmall = 5e-13
    
    particles::Vector{Vector{Particle}} = [create_particles()]
    pia = ParticleIndexerArray(length(particles[1]))

    compute_props!(particles, pia, species_data, phys_props)

    mim = []
    n_moms = 4
    for i in 1:n_moms
        append!(mim, compute_multi_index_moments(i))
    end

    mnnls = NNLSMerge(mim, 30)
    vref = sqrt(2 * k_B * 300.0 / species_data[1].mass)
    result = merge_nnls_based!(rng, mnnls, particles[1], pia, 1, 1, vref)
    
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
    particles2::Vector{Vector{Particle}} = [create_particles2()]
    pia2 = ParticleIndexerArray(length(particles2[1]))
    vref = 1.0

    result = merge_nnls_based!(rng, mnnls2, particles2[1], pia2, 1, 1, vref)
    @test abs(mnnls2.rhs_vector[1] - 1.0) < eps()

    # mean velocity is 0.0
    @test abs(mnnls2.rhs_vector[2] - 0.0) < eps()
    @test abs(mnnls2.rhs_vector[3] - 0.0) < eps()
    @test abs(mnnls2.rhs_vector[4] - 0.0) < eps()

    # mixed second-order moments are also 0.0
    @test abs(mnnls2.rhs_vector[8] - 0.0) < eps()
    @test abs(mnnls2.rhs_vector[9] - 0.0) < eps()
    @test abs(mnnls2.rhs_vector[10] - 0.0) < eps()
end