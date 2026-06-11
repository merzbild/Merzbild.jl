@testset "malloc: bkw variable weight + octree N:2 merging with particles with dim(x)=0,1,2,3" begin


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
        
    seed = 1234
    rng = StableRNG(seed)

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, "Ar")

    interaction_data_path = joinpath(@__DIR__, "..", "data", "pseudo_maxwell.toml")
    interaction_data::Array{Interaction, 2} = load_interaction_data(interaction_data_path, species_data)

    dt_scaled = 0.025

    nv = 30
    np_base = 30^3  # some initial guess on # of particle in simulation

    threshold = 6000
    Ntarget = 5000

    @testset "0D particles, 0D merging" begin
        oc = OctreeN2Merge{0}(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel, max_Nbins=6000)

        T0::Float64 = 273.0
        moments_list = [4, 6, 8, 10]
        
        sigma_ref = π * (interaction_data[1,1].vhs_d^2)
        n_dens = 1e23

        vref = sqrt(2 * k_B * T0 / species_data[1].mass)
        Lref = 1.0 / (n_dens * sigma_ref)
        tref = Lref / vref

        Δt = dt_scaled * tref

        particles = [ParticleVector{0}(np_base)]

        vdf0 = (vx, vy, vz) -> bkw(vx, vy, vz, species_data[1].mass, T0, 0.0)

        n_sampled = sample_on_grid!(rng, vdf0, particles[1], nv, species_data[1].mass, T0, n_dens,
                                    0.0, 1.0, 0.0, 1.0, 0.0, 1.0;
                                    v_mult=3.5, cutoff_mult=3.5, noise=0.0, v_offset=[0.0, 0.0, 0.0])

        pia = ParticleIndexerArray(n_sampled)

        phys_props::PhysProps = PhysProps(1, 1, moments_list, Tref=T0)
        compute_props_with_total_moments!(particles, pia, species_data, phys_props)

        collision_factors::CollisionFactors = CollisionFactors()
        collision_data::CollisionData = CollisionData()

        Fnum = n_dens/n_sampled
        collision_factors.sigma_g_w_max = estimate_sigma_g_w_max(interaction_data[1,1], species_data[1], T0, Fnum)

        Δt::Float64 = dt_scaled * tref
        V::Float64 = 1.0

        merges = false

        for ts in 1:10
            ntc!(rng, collision_factors, collision_data, interaction_data, particles[1], pia, 1, 1, Δt, V)

            if phys_props.np[1,1] > threshold
                merge_octree_N2_based!(rng, oc, particles[1], pia, 1, 1, Ntarget)
                merges = true
            end
            
            compute_props_with_total_moments!(particles, pia, species_data, phys_props)
        end
        @test merges == true

        merges = false

        for ts in 1:30
            bytes_coll = @allocated ntc!(rng, collision_factors, collision_data, interaction_data, particles[1], pia, 1, 1, Δt, V)
            @test bytes_coll == 0

            if phys_props.np[1,1] > threshold
                bytes_merge = @allocated merge_octree_N2_based!(rng, oc, particles[1], pia, 1, 1, Ntarget)
                merges = true
                @test bytes_merge == 0
            end
            
            bytes_props = @allocated compute_props_with_total_moments!(particles, pia, species_data, phys_props)
            @test bytes_props == 0
        end

        @test merges == true
    end

    @testset "1D particles, 1D merging" begin
        # 1d particle x vector
        oc1 = OctreeN2Merge{1}(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel, max_Nbins=6000)
        particles3 = [ParticleVector{1}(np_base)]

        T0 = 273.0
        n_dens = 1e23
        moments_list = [4, 6, 8, 10]
        
        sigma_ref = π * (interaction_data[1,1].vhs_d^2)

        vref = sqrt(2 * k_B * T0 / species_data[1].mass)
        Lref = 1.0 / (n_dens * sigma_ref)
        tref = Lref / vref

        Δt = dt_scaled * tref
        V = 1.0

        phys_props::PhysProps = PhysProps(1, 1, moments_list, Tref=T0)

        vdf0 = (vx, vy, vz) -> bkw(vx, vy, vz, species_data[1].mass, T0, 0.0)

        n_sampled = sample_on_grid!(rng, vdf0, particles3[1], nv, species_data[1].mass, T0, n_dens,
                                    0.0, 1.0, 0.0, 1.0, 0.0, 1.0;
                                    v_mult=3.5, cutoff_mult=3.5, noise=0.0, v_offset=[0.0, 0.0, 0.0])

        pia = ParticleIndexerArray(n_sampled)

        collision_factors::CollisionFactors = CollisionFactors()
        collision_data::CollisionData = CollisionData()

        Fnum = n_dens/n_sampled
        collision_factors.sigma_g_w_max = estimate_sigma_g_w_max(interaction_data[1,1], species_data[1], T0, Fnum)

        merges = false

        for ts in 1:3
            ntc!(rng, collision_factors, collision_data, interaction_data, particles3[1], pia, 1, 1, Δt, V)

            if phys_props.np[1,1] > threshold
                merge_octree_N2_based!(rng, oc1, particles3[1], pia, 1, 1, Ntarget)
                merges = true
            end
            
            compute_props_with_total_moments!(particles3, pia, species_data, phys_props)
        end
        @test merges == true

        merges = false

        for ts in 1:20
            bytes_coll = @allocated ntc!(rng, collision_factors, collision_data, interaction_data, particles3[1], pia, 1, 1, Δt, V)
            @test bytes_coll == 0

            if phys_props.np[1,1] > threshold
                bytes_merge = @allocated merge_octree_N2_based!(rng, oc1, particles3[1], pia, 1, 1, Ntarget)
                merges = true
                @test bytes_merge == 0
            end
            
            bytes_props = @allocated compute_props_with_total_moments!(particles3, pia, species_data, phys_props)
            @test bytes_props == 0
        end

        @test merges == true
    end

    @testset "2D particles, 2D merging" begin
        # 2d particle x vector
        oc2 = OctreeN2Merge{2}(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel, max_Nbins=6000)
        particles3 = [ParticleVector{2}(np_base)]

        T0 = 273.0
        n_dens = 1e23
        moments_list = [4, 6, 8, 10]
        
        sigma_ref = π * (interaction_data[1,1].vhs_d^2)

        vref = sqrt(2 * k_B * T0 / species_data[1].mass)
        Lref = 1.0 / (n_dens * sigma_ref)
        tref = Lref / vref

        Δt = dt_scaled * tref
        V = 1.0

        phys_props::PhysProps = PhysProps(1, 1, moments_list, Tref=T0)

        vdf0 = (vx, vy, vz) -> bkw(vx, vy, vz, species_data[1].mass, T0, 0.0)

        n_sampled = sample_on_grid!(rng, vdf0, particles3[1], nv, species_data[1].mass, T0, n_dens,
                                    0.0, 1.0, 0.0, 1.0, 0.0, 1.0;
                                    v_mult=3.5, cutoff_mult=3.5, noise=0.0, v_offset=[0.0, 0.0, 0.0])

        pia = ParticleIndexerArray(n_sampled)

        collision_factors::CollisionFactors = CollisionFactors()
        collision_data::CollisionData = CollisionData()

        Fnum = n_dens/n_sampled
        collision_factors.sigma_g_w_max = estimate_sigma_g_w_max(interaction_data[1,1], species_data[1], T0, Fnum)

        merges = false

        for ts in 1:3
            ntc!(rng, collision_factors, collision_data, interaction_data, particles3[1], pia, 1, 1, Δt, V)

            if phys_props.np[1,1] > threshold
                merge_octree_N2_based!(rng, oc2, particles3[1], pia, 1, 1, Ntarget)
                merges = true
            end
            
            compute_props_with_total_moments!(particles3, pia, species_data, phys_props)
        end
        @test merges == true

        merges = false

        for ts in 1:20
            bytes_coll = @allocated ntc!(rng, collision_factors, collision_data, interaction_data, particles3[1], pia, 1, 1, Δt, V)
            @test bytes_coll == 0

            if phys_props.np[1,1] > threshold
                bytes_merge = @allocated merge_octree_N2_based!(rng, oc2, particles3[1], pia, 1, 1, Ntarget)
                merges = true
                @test bytes_merge == 0
            end
            
            bytes_props = @allocated compute_props_with_total_moments!(particles3, pia, species_data, phys_props)
            @test bytes_props == 0
        end

        @test merges == true
    end

    @testset "3D particles, 3D merging" begin
        # 3d particle x vector
        oc3 = OctreeN2Merge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel, max_Nbins=6000)
        threshold = 6000
        Ntarget = 5000
        particles3 = [ParticleVector{3}(np_base)]

        T0 = 273.0
        n_dens = 1e23
        moments_list = [4, 6, 8, 10]
        
        sigma_ref = π * (interaction_data[1,1].vhs_d^2)

        vref = sqrt(2 * k_B * T0 / species_data[1].mass)
        Lref = 1.0 / (n_dens * sigma_ref)
        tref = Lref / vref

        Δt = dt_scaled * tref
        V = 1.0

        phys_props::PhysProps = PhysProps(1, 1, moments_list, Tref=T0)

        vdf0 = (vx, vy, vz) -> bkw(vx, vy, vz, species_data[1].mass, T0, 0.0)

        n_sampled = sample_on_grid!(rng, vdf0, particles3[1], nv, species_data[1].mass, T0, n_dens,
                                    0.0, 1.0, 0.0, 1.0, 0.0, 1.0;
                                    v_mult=3.5, cutoff_mult=3.5, noise=0.0, v_offset=[0.0, 0.0, 0.0])

        pia = ParticleIndexerArray(n_sampled)

        collision_factors::CollisionFactors = CollisionFactors()
        collision_data::CollisionData = CollisionData()

        Fnum = n_dens/n_sampled
        collision_factors.sigma_g_w_max = estimate_sigma_g_w_max(interaction_data[1,1], species_data[1], T0, Fnum)

        merges = false

        for ts in 1:3
            ntc!(rng, collision_factors, collision_data, interaction_data, particles3[1], pia, 1, 1, Δt, V)

            if phys_props.np[1,1] > threshold
                merge_octree_N2_based!(rng, oc3, particles3[1], pia, 1, 1, Ntarget)
                merges = true
            end
            
            compute_props_with_total_moments!(particles3, pia, species_data, phys_props)
        end
        @test merges == true

        merges = false

        for ts in 1:20
            bytes_coll = @allocated ntc!(rng, collision_factors, collision_data, interaction_data, particles3[1], pia, 1, 1, Δt, V)
            @test bytes_coll == 0

            if phys_props.np[1,1] > threshold
                bytes_merge = @allocated merge_octree_N2_based!(rng, oc3, particles3[1], pia, 1, 1, Ntarget)
                merges = true
                @test bytes_merge == 0
            end
            
            bytes_props = @allocated compute_props_with_total_moments!(particles3, pia, species_data, phys_props)
            @test bytes_props == 0
        end

        @test merges == true
    end
end