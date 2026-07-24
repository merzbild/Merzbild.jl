@testset "nnls_merging, 2D particles" begin

    function create_particles()
        pv = ParticleVector{2}(50)

        for i in 1:50
            Merzbild.update_particle_buffer_new_particle!(pv, i)
            pv.particles[i] = Particle{2}(1.0, SVector{3,Float64}(i, 5 - 2.0 * i, -100.0 + 0.25 * i^2), SVector{2,Float64}(-0.5 * i, 27 + i))
        end

        return pv
    end

    function create_particles2(; mult=0.0, mult2=1.0)
        pv = ParticleVector{2}(16)

        pv.particles = [Particle{2}(1.0, SVector{3,Float64}(1.0, 1.0, 1.0), SVector{2,Float64}(1.0 * mult, 0.0)),
                        Particle{2}(1.0, SVector{3,Float64}(1.0, 1.0, -1.0), SVector{2,Float64}(0.0, 1.0 * mult)),
                        Particle{2}(1.0, SVector{3,Float64}(1.0, -1.0, 1.0), SVector{2,Float64}(0.0, -1.0 * mult)),
                        Particle{2}(1.0, SVector{3,Float64}(1.0, -1.0, -1.0), SVector{2,Float64}(0.0, 0.0)),
                        Particle{2}(1.0, SVector{3,Float64}(-1.0, 1.0, 1.0), SVector{2,Float64}(0.0, 0.5 * mult)),
                        Particle{2}(1.0, SVector{3,Float64}(-1.0, 1.0, -1.0), SVector{2,Float64}(0.5 * mult, 0.0)),
                        Particle{2}(1.0, SVector{3,Float64}(-1.0, -1.0, 1.0), SVector{2,Float64}(0.0, -0.5 * mult)),
                        Particle{2}(1.0, SVector{3,Float64}(-1.0, -1.0, -1.0), SVector{2,Float64}(0.0, 0.0)),
                        Particle{2}(1.0, SVector{3,Float64}(2.0, 2.0, 2.0), SVector{2,Float64}(0.0, 0.0)),
                        Particle{2}(1.0, SVector{3,Float64}(2.0, 2.0, -2.0), SVector{2,Float64}(0.0, 0.0)),
                        Particle{2}(1.0, SVector{3,Float64}(2.0, -2.0, 2.0), SVector{2,Float64}(0.0, 0.0)),
                        Particle{2}(1.0, SVector{3,Float64}(2.0, -2.0, -2.0), SVector{2,Float64}(0.0, 0.0)),
                        Particle{2}(1.0, SVector{3,Float64}(-2.0, 2.0, 2.0), SVector{2,Float64}(0.0, 0.0)),
                        Particle{2}(1.0*mult2, SVector{3,Float64}(-2.0, 2.0, -2.0), SVector{2,Float64}(0.0, 1.5 * mult)),
                        Particle{2}(1.0*mult2, SVector{3,Float64}(-2.0, -2.0, 2.0), SVector{2,Float64}(0.0, 0.0)),
                        Particle{2}(1.0*mult2, SVector{3,Float64}(-2.0, -2.0, -2.0), SVector{2,Float64}(-2.0 * mult, 0.0))]

        for i in 1:16
            Merzbild.update_particle_buffer_new_particle!(pv, i)
        end
        
        return pv
    end

    function create_particles3()
        pv = ParticleVector{2}(6)

        pv.particles = [Particle{2}(1.0, SVector{3,Float64}(3.0, 1.0, 1.0), SVector{2,Float64}(0.0, 0.0)),
                        Particle{2}(1.0, SVector{3,Float64}(1.0, 4.0, -1.0), SVector{2,Float64}(0.0, 0.0)),
                        Particle{2}(1.0, SVector{3,Float64}(1.0, -8.0, 1.0), SVector{2,Float64}(0.0, 0.0)),
                        Particle{2}(1.0, SVector{3,Float64}(1.0, -1.0, -9.0), SVector{2,Float64}(0.0, 0.0)),
                        Particle{2}(1.0, SVector{3,Float64}(-7.0, -11.0, 1.0), SVector{2,Float64}(0.0, 0.0)),
                        Particle{2}(1.0, SVector{3,Float64}(-9.0, 6.0, -1.0), SVector{2,Float64}(0.0, 0.0))]

        for i in 1:6
            Merzbild.update_particle_buffer_new_particle!(pv, i)
        end
        
        return pv
    end

    function create_particles4()
        pv = ParticleVector{2}(2)

        pv.particles = [Particle{2}(1.0, SVector{3,Float64}(3.0, 1.0, -4.0), SVector{2,Float64}(0.0, 0.0)),
                        Particle{2}(1.0, SVector{3,Float64}(-3.0, 1.0, 4.0), SVector{2,Float64}(0.0, 0.0)),]

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

    T_mom_ref = 1
    phys_props::PhysProps = PhysProps(1, 1)

    mlist = [4]
    mscaling = zeros(1)
    mvals_list = zeros(1)
    compute_moment_scaling!(mscaling, mlist, 1, species_data, T_mom_ref)

    Δabs = 2.5
    Δrel_xsmall = 5e-13
    
    particles = [create_particles()]
    pia = ParticleIndexerArray(length(particles[1]))

    compute_props!(particles, pia, species_data, phys_props)
    compute_moments!(mvals_list, mscaling, mlist, particles, pia, 1, 1, species_data, phys_props)

    mim = []
    n_moms = 4
    for i in 1:n_moms
        append!(mim, compute_multi_index_moments(i))
    end

    mnnls = NNLSMerge{2}(mim, 30)

    vref = sqrt(2 * k_B * 300.0 / species_data[1].mass)
    result = merge_nnls_based!(rng, mnnls, particles[1], pia, 1, 1; vref=vref, scaling=:vref)
    @test length(mnnls.rhs_vector) == length(mim) + 1  # we didn't construct 0-th order moment

    n0 = phys_props.n[1, 1]
    np0 = phys_props.np[1, 1]
    v0 = phys_props.v[:, 1, 1]
    T0 = phys_props.T[1, 1]
    M40 = mvals_list[1]

    mixed_order_3 = compute_multi_index_moments(3)
    moms3 = zeros(length(mixed_order_3))

    for (i, mim) in enumerate(mixed_order_3)
        moms3[i] = Merzbild.compute_mixed_moment(particles, pia, 1, 1, mim)
    end

    # test some internal computes used in the merging
    @test sum(abs.(mnnls.v0 .- v0)) < eps()
    @test sum(abs.(mnnls.w_total .- n0)) < 2.5e-14

    compute_props!(particles, pia, species_data, phys_props)
    compute_moments!(mvals_list, mscaling, mlist, particles, pia, 1, 1, species_data, phys_props)
    # test that merging conserves mass / momentum / energy
    @test phys_props.np[1, 1] < np0
    @test abs(n0 - phys_props.n[1, 1]) < 2.5e-14
    @test sum(abs.(v0 - phys_props.v[:, 1, 1])) < 2e-13
    @test abs(T0 - phys_props.T[1, 1]) < 3.6e-14

    # test higher-order moment in g
    @test abs(M40 - mvals_list[1]) < 7.4e-12

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

    mnnls2 = NNLSMerge{2}(mim, 30)
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
    mnnls = NNLSMerge{2}([], 30)
    @test [0,0,0] in mnnls.mim
    @test [1,0,0] in mnnls.mim
    @test [0,1,0] in mnnls.mim
    @test [0,0,1] in mnnls.mim
    @test [2,0,0] in mnnls.mim
    @test [0,2,0] in mnnls.mim
    @test [0,0,2] in mnnls.mim
    @test mnnls.n_moments_vel == 7
    @test mnnls.n_moments_pos == 0
    @test mnnls.pos_i_x == -1
    @test mnnls.pos_i_y == -1
    @test mnnls.pos_i_z == -1

    # check computation of LHS/RHS of the system

    particles = [create_particles3()]
    pia = ParticleIndexerArray(length(particles[1]))

    lhs_ncols = pia.indexer[1, 1].n_local
    lhs_matrix = zeros(mnnls.n_moments_vel, lhs_ncols)
    vel_pos_matrix = zeros(6, lhs_ncols)

    Merzbild.compute_lhs_and_rhs!(mnnls, lhs_matrix, vel_pos_matrix,
                                  particles[1], pia, 1, 1)

    max_v_diff = 0.0
    max_x_diff = 0.0
    for i in 1:length(particles[1])
        max_v_diff = max(abs.(vel_pos_matrix[1,i] - particles[1][i].v[1]),
                         abs.(vel_pos_matrix[2,i] - particles[1][i].v[2]),
                         abs.(vel_pos_matrix[3,i] - particles[1][i].v[3]),
                         max_v_diff)
        max_x_diff = max(abs.(vel_pos_matrix[4,i] - particles[1][i].x[1]),
                         abs.(vel_pos_matrix[5,i] - particles[1][i].x[2]),
                         max_x_diff)
    end
    @test max_v_diff <= 2*eps()
    @test max_x_diff <= 2*eps()

    mim = []
    #  [[0, 0, 0],
    # [1, 0, 0], [0, 1, 0], [0, 0, 1],
    # [2, 0, 0], [0, 2, 0], [0, 0, 2]]
    mnnls3 = NNLSMerge{2}(mim, 30)

    lhs_ncols = 2
    # 7, 
    lhs_matrix = zeros(mnnls3.n_moments_vel, lhs_ncols)
    vel_pos_matrix = zeros(6, lhs_ncols)

    particles = [create_particles4()]
    pia = ParticleIndexerArray(length(particles[1]))

    Merzbild.compute_lhs_and_rhs!(mnnls3, lhs_matrix, vel_pos_matrix,
                                  particles[1], pia, 1, 1)

    lhs_matrix_2 = copy(lhs_matrix)
    rhs_vec_2 = copy(mnnls3.rhs_vector)
    mnnls3.vref = 2.0
    mnnls3.inv_vref = 1.0 / mnnls3.vref
    Merzbild.scale_lhs_rhs!(mnnls3, lhs_matrix, :vref, size(lhs_matrix, 2))

    @test maximum(abs.(mnnls3.scalev - [2.0, 2.0, 2.0])) < 2*eps()
    
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

    mnnls3.Ev = SVector{3,Float64}(2.0, 3.0, 4.0)

    Merzbild.scale_lhs_rhs!(mnnls3, lhs_matrix, :variance, size(lhs_matrix, 2))

    @test abs(rhs_vec_2[1] - mnnls3.rhs_vector[1]) < 2*eps()
    @test abs(rhs_vec_2[2] / 2.0 - mnnls3.rhs_vector[2]) < 2*eps()
    @test abs(rhs_vec_2[3] / 3.0 - mnnls3.rhs_vector[3]) < 2*eps()
    @test abs(rhs_vec_2[4] / 4.0 - mnnls3.rhs_vector[4]) < 2*eps()
    @test abs(rhs_vec_2[5] / 4.0 - mnnls3.rhs_vector[5]) < 2*eps()
    @test abs(rhs_vec_2[6] / 9.0 - mnnls3.rhs_vector[6]) < 2*eps()
    @test abs(rhs_vec_2[7] / 16.0 - mnnls3.rhs_vector[7]) < 2*eps()

    @test maximum(abs.(mnnls3.scalev - [2.0, 3.0, 4.0])) < 2*eps()
    for i in 1:2
        @test maximum(abs.(lhs_matrix_2[1,i] - lhs_matrix[1,i])) < 2*eps()
        @test maximum(abs.(lhs_matrix_2[2,i] / 2.0 - lhs_matrix[2,i])) < 2*eps()
        @test maximum(abs.(lhs_matrix_2[3,i] / 3.0 - lhs_matrix[3,i])) < 2*eps()
        @test maximum(abs.(lhs_matrix_2[4,i] / 4.0 - lhs_matrix[4,i])) < 2*eps()
        @test maximum(abs.(lhs_matrix_2[5,i] / 4.0 - lhs_matrix[5,i])) < 2*eps()
        @test maximum(abs.(lhs_matrix_2[6,i] / 9.0 - lhs_matrix[6,i])) < 2*eps()
        @test maximum(abs.(lhs_matrix_2[7,i] / 16.0 - lhs_matrix[7,i])) < 2*eps()
    end


    # test conservation of 1st spatial moment
    particles = [create_particles2(mult=2.0,mult2=2.0)]
   
    # first one is discarded
    pos_moments = [[0,0,1],[1,0,0],[0,1,0]]

    mnnls_pos = NNLSMerge{2}([], 30; multi_index_moments_pos=pos_moments)
    @test mnnls_pos.pos_i_x == 1
    @test mnnls_pos.pos_i_y == 2
    @test mnnls_pos.pos_i_z == -1

    x_mean = [0.0, 0.0]
    v_mean = [0.0, 0.0, 0.0]
    tw = 0.0
    for i in 1:16
        x_mean += particles[1][i].x * particles[1][i].w
        v_mean += particles[1][i].v * particles[1][i].w
        tw += particles[1][i].w
    end

    x_mean = x_mean / tw
    v_mean = v_mean / tw

    v_std = [0.0, 0.0, 0.0]
    x_std = [0.0, 0.0]
    for i in 1:16
        v_std += particles[1][i].w * (particles[1][i].v - v_mean).^2
        x_std += particles[1][i].w * (particles[1][i].x - x_mean).^2
    end

    v_std = sqrt.(v_std / tw)
    x_std = sqrt.(x_std / tw)

    pia = ParticleIndexerArray(length(particles[1]))
    result = merge_nnls_based!(rng, mnnls_pos, particles[1], pia, 1, 1)
    
    @test length(mnnls_pos.rhs_vector) == 9
    @test mnnls_pos.n_total_conserved == mnnls_pos.n_moments_vel + mnnls_pos.n_moments_pos
    
    @test result == 1
    @test pia.indexer[1,1].n_local[1] < 16
    @test maximum(abs.(mnnls_pos.x0 - x_mean)) < 4*eps()

    n_new = pia.indexer[1,1].n_local[1]

    @test maximum(abs.(mnnls_pos.Ex - x_std)) < 4*eps()
    @test maximum(abs.(mnnls_pos.scalex - x_std)) < 4*eps()

    # we only conserve mean position
    @test maximum(abs.(mnnls_pos.rhs_vector[8:9])) < 4*eps() 
    
    new_w = 0.0
    new_v = [0.0, 0.0, 0.0]
    new_x = [0.0, 0.0]
    s1 = pia.indexer[1,1].start1
    e1 = pia.indexer[1,1].end1
    for i in s1:e1
        new_w += particles[1][i].w
        new_v += particles[1][i].v * particles[1][i].w
        new_x += particles[1][i].x * particles[1][i].w
    end

    if pia.indexer[1,1].n_group2 > 0
        s2 = pia.indexer[1,1].start2
        e2 = pia.indexer[1,1].end2
        for i in s2:e2
            new_w += particles[1][i].w
            new_v += particles[1][i].v * particles[1][i].w
            new_x += particles[1][i].x * particles[1][i].w
        end
    end

    new_v /= new_w
    new_x /= new_w

    @test abs(new_w - tw) < 3.75e-15
    @test maximum(abs.(new_v - v_mean)) < 1e-15
    @test maximum(abs.(new_x - x_mean)) < 1e-15

    # now we try a higher threshold so that we get fewer particles
    # this will mess up everything except mass conservation though
    particles = [create_particles2(mult=2.0,mult2=2.0)]
    pia = ParticleIndexerArray(length(particles[1]))
    result = merge_nnls_based!(rng, mnnls_pos, particles[1], pia, 1, 1; w_threshold=0.1)

    @test result == 1
    @test pia.indexer[1,1].n_local[1] < n_new

    new_w = 0.0
    new_v = [0.0, 0.0, 0.0]
    new_x = [0.0, 0.0]
    s1 = pia.indexer[1,1].start1
    e1 = pia.indexer[1,1].end1
    for i in s1:e1
        new_w += particles[1][i].w
        # new_v += particles[1][i].v * particles[1][i].w
        # new_x += particles[1][i].x * particles[1][i].w
    end

    if pia.indexer[1,1].n_group2 > 0
        s2 = pia.indexer[1,1].start2
        e2 = pia.indexer[1,1].end2
        for i in s2:e2
            new_w += particles[1][i].w
            # new_v += particles[1][i].v * particles[1][i].w
            # new_x += particles[1][i].x * particles[1][i].w
        end
    end

    # new_v /= new_w
    # new_x /= new_w

    @test abs(new_w - tw) < 4.5e-15
    # @test maximum(abs.(new_v - v_mean)) < 4e-15
    # @test maximum(abs.(new_x - x_mean)) < 4e-15

    # test that computation of spatial moments does not affect the LHS and RHS entries corresponding to the velocity moments
    particles = [create_particles2(mult=2.0,mult2=2.0)]
    pia = ParticleIndexerArray(length(particles[1]))
    lhs_ncols = 16
   
    add_vel_moments = [[3,1,0],[0,2,1]]
    pos_moments = [[1,0,0],[0,1,0],[2,0,2]]

    mnnls_pos = NNLSMerge{2}(add_vel_moments, 30; multi_index_moments_pos=pos_moments)
    mnnls_no_pos = NNLSMerge{2}(add_vel_moments, 30)
    @test mnnls_pos.n_total_conserved == mnnls_pos.n_moments_vel + mnnls_pos.n_moments_pos
    @test mnnls_no_pos.n_total_conserved == mnnls_no_pos.n_moments_vel + mnnls_no_pos.n_moments_pos

    @test mnnls_pos.n_moments_vel == mnnls_no_pos.n_moments_vel
    @test mnnls_no_pos.n_moments_pos == 0
    @test mnnls_pos.n_moments_pos == 3
    @test mnnls_pos.pos_i_x == 1
    @test mnnls_pos.pos_i_y == 2
    @test mnnls_pos.pos_i_z == -1
    @test length(mnnls_pos.rhs_vector) == length(mnnls_no_pos.rhs_vector) + 3 == 12

    lhs_matrix_pos = zeros(mnnls_pos.n_moments_vel + mnnls_pos.n_moments_pos, lhs_ncols)
    lhs_matrix_no_pos = zeros(mnnls_no_pos.n_moments_vel, lhs_ncols)
    vel_pos_matrix = zeros(6, lhs_ncols)
    Merzbild.compute_lhs_and_rhs!(mnnls_pos, lhs_matrix_pos, vel_pos_matrix,
                                  particles[1], pia, 1, 1)
    Merzbild.compute_lhs_and_rhs!(mnnls_no_pos, lhs_matrix_no_pos, vel_pos_matrix,
                                  particles[1], pia, 1, 1)

    @test maximum(abs.(mnnls_pos.rhs_vector[1:mnnls_pos.n_moments_vel] - mnnls_no_pos.rhs_vector)) < 2*eps()
    @test maximum(abs.(lhs_matrix_pos[1:mnnls_pos.n_moments_vel,:] - lhs_matrix_no_pos)) < 2*eps()

    # check computation of central spatial moments
    for i in 1:16
        @test abs(lhs_matrix_pos[mnnls_pos.n_moments_vel+1,i] - (particles[1][i].x[1] - mnnls_pos.x0[1])) < 2*eps()
        @test abs(lhs_matrix_pos[mnnls_pos.n_moments_vel+2,i] - (particles[1][i].x[2] - mnnls_pos.x0[2])) < 2*eps()
        @test abs(lhs_matrix_pos[mnnls_pos.n_moments_vel+3,i] - (particles[1][i].x[1] - mnnls_pos.x0[1])^2) < 2*eps()
    end

    old_val = lhs_matrix_pos[mnnls_pos.n_moments_vel+2,3]
    old_val2 = lhs_matrix_pos[mnnls_pos.n_moments_vel+3,3]

    # scaling of spatial moments also doesn't interfere with scaling of velocity moments
    Merzbild.scale_lhs_rhs!(mnnls_pos, lhs_matrix_pos, :vref, size(lhs_matrix_pos, 2))
    Merzbild.scale_lhs_rhs!(mnnls_no_pos, lhs_matrix_no_pos, :vref, size(lhs_matrix_no_pos, 2))

    @test maximum(abs.(mnnls_pos.rhs_vector[1:mnnls_pos.n_moments_vel] - mnnls_no_pos.rhs_vector)) < 2*eps()
    @test maximum(abs.(lhs_matrix_pos[1:mnnls_pos.n_moments_vel,:] - lhs_matrix_no_pos)) < 2*eps()

    # test correct scaling of spatial moments
    @test abs(lhs_matrix_pos[mnnls_pos.n_moments_vel+2,3] * x_std[2] - old_val) < 2*eps()
    @test abs(lhs_matrix_pos[mnnls_pos.n_moments_vel+3,3] * x_std[1]^2 - old_val2) < 2*eps()

    # finally check that we write 0.0 + mean_position to the position data in case 1st moment not specified
    result = merge_nnls_based!(rng, mnnls_pos, particles[1], pia, 1, 1; w_threshold=1e-12)
    @test result == 1

    mnnls_pos = NNLSMerge{2}(add_vel_moments, 25; multi_index_moments_pos=pos_moments, matrix_ncol_nprealloc=10)
    @test length(mnnls_pos.vel_pos_matrices) == 11
    @test length(mnnls_pos.column_norms) == 11
    @test length(mnnls_pos.work) == 12
    @test mnnls_pos.lhs_matrix_ncols_start == 25
    @test mnnls_pos.lhs_matrix_ncols_end == 35

    # scratch buffers for column counts outside the pre-allocated range
    @test size(mnnls_pos.vel_pos_matrix_scratch) == (3+2, 36)
    @test length(mnnls_pos.column_norms_scratch) == 36

    for i in 25:35
        @test size(mnnls_pos.vel_pos_matrices[i-25+1]) == (3+2, i)
        @test length(mnnls_pos.column_norms[i-25+1]) == i
        @test size(mnnls_pos.work[i-25+1].QA,2) == i

        @test size(mnnls_pos.work[i-25+1].x) == (i,)
        @test size(mnnls_pos.work[i-25+1].w) == (i,)
        @test size(mnnls_pos.work[i-25+1].idx) == (i,)
    end

    # test merging with pre-allocated NNLS work arrays
    particles = [create_particles2(mult=2.0,mult2=2.0)]
    pia = ParticleIndexerArray(length(particles[1]))
    result = merge_nnls_based!(rng, mnnls_pos, particles[1], pia, 1, 1; w_threshold=1e-12)
    @test result == 1

    particles = [create_particles2(mult=2.0,mult2=2.0)]
    pia = ParticleIndexerArray(length(particles[1]))
    result = merge_nnls_based!(rng, mnnls_pos, particles[1], pia, 1, 1; w_threshold=1e-12)
    @test result == 1
    @test pia.index_last[1] == pia.n_total[1]

    wtot = 0.0
    s1 = pia.indexer[1,1].start1
    e1 = pia.indexer[1,1].end1
    for i in s1:e1
        wtot += particles[1][i].w
    end

    if pia.indexer[1,1].n_group2 > 0
        s2 = pia.indexer[1,1].start2
        e2 = pia.indexer[1,1].end2
        for i in s2:e2
            wtot += particles[1][i].w
        end
    end

    @test abs(wtot - tw) < 4.5e-15

    # finally test the edge case of init_np = total_np and no pre-allocated matrices, v_multipliers=[]
    # reset particles
    particles = [create_particles2(mult=2.0,mult2=2.0)]
    pia = ParticleIndexerArray(length(particles[1]))

    mnnls = NNLSMerge{2}(add_vel_moments, length(particles[1]); multi_index_moments_pos=[])
    result = merge_nnls_based!(rng, mnnls, particles[1], pia, 1, 1; w_threshold=1e-12)

    @test result == 1
    @test pia.n_total[1] < 16
end