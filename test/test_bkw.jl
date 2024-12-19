@testset "bkw" begin


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
        
    seed = 1234
    rng = StableRNG(seed)

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, "Ar")

    interaction_data_path = joinpath(@__DIR__, "..", "data", "pseudo_maxwell.toml")
    interaction_data::Array{Interaction, 2} = load_interaction_data(interaction_data_path, species_data)

    dt_scaled = 0.025
    n_t = 500
    n_particles = 20000

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

    Fnum::Float64 = n_dens / n_particles

    particles::Vector{ParticleVector} = [ParticleVector(n_particles)]

    pia = ParticleIndexerArray(0)

    sample_particles_equal_weight!(rng, particles[1], pia, 1, 1, n_particles,
                                   species_data[1].mass, T0, Fnum, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0;
                                   distribution=:BKW)


    phys_props::PhysProps = PhysProps(1, 1, moments_list, Tref=T0)
    compute_props_with_total_moments!(particles, pia, species_data, phys_props)

    sol_path = joinpath(@__DIR__, "data", "tmp_bkw.nc")
    ds = NCDataHolder(sol_path, species_data, phys_props)
    write_netcdf_phys_props(ds, phys_props, 0)

    collision_factors::CollisionFactors = CollisionFactors()
    collision_data::CollisionData = CollisionData()

    collision_factors.sigma_g_w_max = estimate_sigma_g_w_max(interaction_data[1,1], species_data[1], T0, Fnum)

    Δt::Float64 = dt_scaled * tref
    V::Float64 = 1.0

    for ts in 1:n_t
        ntc!(rng, collision_factors, collision_data, interaction_data, particles[1], pia, 1, 1, Δt, V)
        
        compute_props_with_total_moments!(particles, pia, species_data, phys_props)
        write_netcdf_phys_props(ds, phys_props, ts)
    end
    close_netcdf(ds)

    ref_sol_path = joinpath(@__DIR__, "data", "bkw_20k_seed1234.nc")
    ref_sol = NCDataset(ref_sol_path, "r")
    sol = NCDataset(sol_path, "r")

    @test length(sol["timestep"]) == n_t + 1

    @test maximum(sol["np"]) == n_particles
    @test minimum(sol["np"]) == n_particles

    ref_mom = ref_sol["moments"]
    sol_mom = sol["moments"]

    for mom_no in 1:length(moments_list)
        diff = abs.(ref_mom[mom_no, 1, 1, :] - sol_mom[mom_no, 1, 1, :])
        @test maximum(diff) <= 2 * eps()
    end

    close(ref_sol)

    analytic_4 = analytic(sol["timestep"] * dt_scaled, magic_factor, 4)
    diff = abs.(analytic_4 .- sol_mom[1, 1, 1, :]) ./ analytic_4
    @test maximum(diff) < 0.05

    analytic_6 = analytic(sol["timestep"] * dt_scaled, magic_factor, 6)
    diff = abs.(analytic_6 .- sol_mom[2, 1, 1, :]) ./ analytic_6
    @test maximum(diff) < 0.055

    analytic_8 = analytic(sol["timestep"] * dt_scaled, magic_factor, 8)
    diff = abs.(analytic_8 .- sol_mom[3, 1, 1, :]) ./ analytic_8
    @test maximum(diff) < 0.15

    close(sol)
    rm(sol_path)
end