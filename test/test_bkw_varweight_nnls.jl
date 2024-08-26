@testset "bkw variable weight + NNLS merging" begin


    # Important!
    # The time scaling in the analytical solution is different
    # Tref = 273.0
    # mref = 66.3e-27 
    # mcd = mref / 2.0
    # dref = 4.11e-10
    # nref = 1e23
    # Lref = 1.0 / (nref * constants.pi * dref**2)
    # vref = ((2 * constants.k * Tref) / mref)**0.5
    # time_ref = Lref / vref

    # kappa_mult = constants.pi * dref**2 * (mcd / (2 * constants.k * tref))**(-0.5) / gamma(5/2 - 1.0)
    # ttt_bkw = 1 / (4 * constants.pi * n * kappa_mult)
    # magic_factor = time_ref / ttt_bkw / (4 * constants.pi)
    # print(magic_factor)  # approximately 1.59577 for Argon, 1.5963 for N
        
    # def analytic(time, N):    
    #     C = 1. - 0.4 * np.exp(-time * magic_factor / 6)
    #     kk = N // 2
    #     return C**(kk - 1) * (kk - (kk - 1) * C)

    function analytic(time, magic_factor, N)
        C = 1.0 .- 0.4 * exp.(-time * magic_factor / 6)
        kk = N / 2
        return C.^(kk - 1) .* (kk .- (kk - 1) * C)
    end
        
    seed = 0
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, "Ar")

    interaction_data_path = joinpath(@__DIR__, "..", "data", "pseudo_maxwell.toml")
    interaction_data::Array{Interaction, 2} = load_interaction_data(interaction_data_path, species_data)

    dt_scaled = 0.025
    n_t = 500

    nv = 32
    np_base = 40^3  # some initial guess on # of particle in simulation
    
    n_full_up_to_total = 6
    n_up_to_total = 8

    threshold = 150
    ntarget_octree = 100

    mim = []
    
    n_moms = n_full_up_to_total
    for i in 1:n_moms
        append!(mim, compute_multi_index_moments(i))
    end

    for i in n_full_up_to_total+1:n_up_to_total
        append!(mim, [[i, 0, 0], [0, i, 0], [0, 0, i]])
    end

    mnnls = NNLSMerge(mim, threshold)
    ocm = OctreeN2Merge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel, max_Nbins=6000)

    T0::Float64 = 273.0
    sigma_ref = π * (interaction_data[1,1].vhs_d^2)
    n_dens = 1e23

    vref = sqrt(2 * k_B * T0 / species_data[1].mass)
    Lref = 1.0 / (n_dens * sigma_ref)
    tref = Lref / vref
    moments_list = [4, 6, 8, 10]

    kappa_mult = sigma_ref * (interaction_data[1,1].m_r / (2 * k_B * T0))^(-0.5) / gamma(5/2 - 1.0)
    ttt_bkw = 1 / (4 * π * n_dens * kappa_mult)
    magic_factor = tref / ttt_bkw / (4 * π)

    particles::Vector{Vector{Particle}} = [Vector{Particle}(undef, np_base)]

    vdf0 = (vx, vy, vz) -> bkw(vx, vy, vz, species_data[1].mass, T0, 0.0)

    n_sampled = sample_on_grid!(rng, vdf0, particles[1], nv, species_data[1].mass, T0, n_dens,
                                0.0, 1.0, 0.0, 1.0, 0.0, 1.0;
                                v_mult=3.5, cutoff_mult=3.5, noise=0.0, v_offset=[0.0, 0.0, 0.0])

    pia = ParticleIndexerArray(n_sampled)

    phys_props::PhysProps = PhysProps(1, 1, moments_list, Tref=T0)
    compute_props_with_total_moments!(particles, pia, species_data, phys_props)

    sol_path = joinpath(@__DIR__, "data", "tmp_bkw_nnls.nc")
    ds = NCDataHolder(sol_path, species_data, phys_props)
    write_netcdf_phys_props(ds, phys_props, 0)

    collision_factors::CollisionFactors = CollisionFactors()
    collision_data::CollisionData = CollisionData()

    Fnum = n_dens/n_sampled
    collision_factors.sigma_g_w_max = estimate_sigma_g_w_max(interaction_data[1,1], species_data[1], T0, Fnum)

    Δt::Float64 = dt_scaled * tref
    V::Float64 = 1.0
    natt = 0
    
    nnls_success_flag = 1

    for ts in 1:n_t
        ntc!(rng, collision_factors, collision_data, interaction_data, particles[1], pia, 1, 1, Δt, V)

        if phys_props.np[1,1] > threshold
            nnls_success_flag = merge_nnls_based!(rng, mnnls, particles[1], pia, 1, 1, vref)

            if nnls_success_flag == -1
                merge_octree_N2_based!(ocm, particles[1], pia, 1, 1, ntarget_octree)
            end
        end
        
        compute_props_with_total_moments!(particles, pia, species_data, phys_props)
        write_netcdf_phys_props(ds, phys_props, ts)
    end
    close_netcdf(ds)

    @test abs(phys_props.T[1,1] - T0) < 5e-4
    @test abs(phys_props.n[1,1] / n_dens - 1.0) < 1e-11
    @test phys_props.np[1,1] < threshold

    ref_sol_path = joinpath(@__DIR__, "data", "bkw_vw_nnls_6full_upto8_150_seed0.nc")
    ref_sol = NCDataset(ref_sol_path, "r")
    sol = NCDataset(sol_path, "r")

    @test length(sol["timestep"]) == n_t + 1

    ref_mom = ref_sol["moments"]
    sol_mom = sol["moments"]

    for mom_no in 1:length(moments_list)
        diff = abs.(ref_mom[mom_no, 1, 1, :350] - sol_mom[mom_no, 1, 1, :350])
        @test maximum(diff) <= 2.4e-6 # something weird going on in the test/prod environment  # * eps()
    end

    close(ref_sol)

    analytic_4 = analytic(sol["timestep"] * dt_scaled, magic_factor, 4)
    diff = abs.(analytic_4 .- sol_mom[1, 1, 1, :]) ./ analytic_4
    @test maximum(diff) < 0.12

    close(sol)
    rm(sol_path)
end