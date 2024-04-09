@testset "bkw" begin

    function analytic(time, magic_factor, N)
        C = 1.0 .- 0.4 * exp.(-time * magic_factor / 6)
        kk = N / 2
        return C.^(kk - 1) .* (kk .- (kk - 1) * C)
    end
        
    seed = 1234
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    species_list::Vector{Species} = load_species_list("data/particles.toml", "Ar")
    interaction_data::Array{Interaction, 2} = load_interaction_data("data/pseudo_maxwell.toml", species_list)

    dt_scaled = 0.025
    n_t = 500
    n_particles = 20000

    T0::Float64 = 273.0
    sigma_ref = π * (interaction_data[1,1].vhs_d^2)
    n_dens = 1e23

    vref = sqrt(2 * k_B * T0 / species_list[1].mass)
    Lref = 1.0 / (n_dens * sigma_ref)
    tref = Lref / vref
    moments_list = [4, 6, 8, 10]

    kappa_mult = sigma_ref * (interaction_data[1,1].m_r / (2 * k_B * T0))^(-0.5) / gamma(5/2 - 1.0)
    ttt_bkw = 1 / (4 * π * n_dens * kappa_mult)
    magic_factor = tref / ttt_bkw / (4 * π)

    Fnum::Float64 = n_dens / n_particles

    particles::Vector{Vector{Particle}} = [Vector{Particle}(undef, n_particles)]
    sample_particles_equal_weight!(rng, particles[1], n_particles, T0, species_list[1].mass, Fnum,
                                0.0, 1.0, 0.0, 1.0, 0.0, 1.0; distribution=:BKW)

    particle_indexer::Array{ParticleIndexer, 2} = Array{ParticleIndexer, 2}(undef, 1, 1)
    particle_indexer[1,1] = create_particle_indexer(n_particles)

    phys_props::PhysProps = create_props(1, 1, moments_list, Tref=T0)
    compute_props!(phys_props, particle_indexer, particles, species_list)

    ds = create_netcdf_phys_props("test/data/tmp_bkw.nc", phys_props, species_list)
    write_netcdf_phys_props(ds, phys_props, 0)

    collision_factors::CollisionFactors = create_collision_factors()
    collision_data::CollisionData = create_collision_data()

    collision_factors.sigma_g_w_max = estimate_sigma_g_w_max(interaction_data[1,1], species_list[1], T0, Fnum)

    Δt::Float64 = dt_scaled * tref
    V::Float64 = 1.0

    for ts in 1:n_t
        ntc!(rng, collision_factors, particle_indexer[1,1], collision_data, interaction_data[1,1], particles[1],
            Δt, V)
        
        compute_props!(phys_props, particle_indexer, particles, species_list)
        write_netcdf_phys_props(ds, phys_props, ts)
    end
    close(ds)

    ref_sol = NCDataset("test/data/bkw_20k_seed1234.nc", "r")
    sol = NCDataset("test/data/tmp_bkw.nc", "r")

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
    @test maximum(diff) < 0.05

    analytic_8 = analytic(sol["timestep"] * dt_scaled, magic_factor, 8)
    diff = abs.(analytic_8 .- sol_mom[3, 1, 1, :]) ./ analytic_8
    @test maximum(diff) < 0.15

    close(sol)
    rm("test/data/tmp_bkw.nc")
end