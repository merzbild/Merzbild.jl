@testset "VSS scattering" begin
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    vss_data_path = joinpath(@__DIR__, "..", "data", "vss.toml")

    species_data = load_species_data(particles_data_path, ["Ar", "He"])
    interaction_data = load_interaction_data(vss_data_path, species_data)
    interaction = interaction_data[1, 2]
    α = interaction.vss_alpha

    rng = StableRNG(1234)

    p1 = Particle(1e10, [700.0, -200.0, 350.0], [0.0, 0.0, 0.0])
    p2 = Particle(1e10, [-100.0, 400.0, 50.0], [0.0, 0.0, 0.0])

    collision_data = CollisionData()
    Merzbild.compute_com!(collision_data, interaction, p1, p2)
    Merzbild.compute_g!(collision_data, p1, p2)

    v_com = collision_data.v_com
    g_vec = collision_data.g_vec
    g = collision_data.g

    # orthonormal frame with n_par along the pre-collisional relative velocity
    n_par = g_vec / g
    gyz = sqrt(g_vec[2]^2 + g_vec[3]^2)
    n_perp1 = SVector{3,Float64}(0.0, g_vec[3], -g_vec[2]) / gyz
    n_perp2 = SVector{3,Float64}(n_par[2] * n_perp1[3] - n_par[3] * n_perp1[2],
                                 n_par[3] * n_perp1[1] - n_par[1] * n_perp1[3],
                                 n_par[1] * n_perp1[2] - n_par[2] * n_perp1[1])

    collision_data_check = CollisionData()

    n_samples = 1000000

    n_cchi_bins = 40
    cchi_counts = zeros(Int64, n_cchi_bins)

    n_eps_bins = 40
    eps_counts = zeros(Int64, n_eps_bins)

    mean_cchi = 0.0
    mean_cchi2 = 0.0

    max_g_err = 0.0
    max_v_com_err = 0.0

    for _ in 1:n_samples
        Merzbild.scatter_vss!(rng, collision_data, interaction, p1, p2)

        n_vec = collision_data.g_vec_new / g

        cchi = n_vec[1] * n_par[1] + n_vec[2] * n_par[2] + n_vec[3] * n_par[3]
        mean_cchi += cchi
        mean_cchi2 += cchi * cchi

        cchi_index = min(n_cchi_bins, 1 + floor(Int64, 0.5 * (cchi + 1.0) * n_cchi_bins))
        cchi_counts[cchi_index] += 1

        # the azimuthal angle around the pre-collisional relative velocity is uniform
        proj1 = n_vec[1] * n_perp1[1] + n_vec[2] * n_perp1[2] + n_vec[3] * n_perp1[3]
        proj2 = n_vec[1] * n_perp2[1] + n_vec[2] * n_perp2[2] + n_vec[3] * n_perp2[3]
        eps_index = min(n_eps_bins, 1 + floor(Int64, (atan(proj2, proj1) + pi) * n_eps_bins / (2 * pi)))
        eps_counts[eps_index] += 1

        # scattering is elastic: |g| and v_com are unchanged
        Merzbild.compute_g!(collision_data_check, p1, p2)
        Merzbild.compute_com!(collision_data_check, interaction, p1, p2)
        max_g_err = max(max_g_err, abs(collision_data_check.g - g) / g)
        max_v_com_err = max(max_v_com_err, maximum(abs.(collision_data_check.v_com - v_com)) / g)
    end

    mean_cchi /= n_samples
    mean_cchi2 /= n_samples

    @test max_g_err < 1e-13
    @test max_v_com_err < 1e-13

    # <cos(chi)> = 2 α / (1 + α) - 1, <cos^2(chi)> = 4 α / (2 + α) - 4 α / (1 + α) + 1
    # the latter defines the viscosity cross-section σ_μ = σ_T (1 - <cos^2(chi)>)
    @test abs(mean_cchi - (2 * α / (1 + α) - 1)) < 5e-3
    @test abs(mean_cchi2 - (4 * α / (2 + α) - 4 * α / (1 + α) + 1)) < 5e-3

    # cos(chi) = 2 R^(1/α) - 1, so P(cos(chi) < c) = ((1 + c) / 2)^α
    chi2 = 0.0
    for index in 1:n_cchi_bins
        cdf_low = (index - 1.0)^α / n_cchi_bins^α
        cdf_high = index^α / n_cchi_bins^α
        expected_count = n_samples * (cdf_high - cdf_low)
        chi2 += (cchi_counts[index] - expected_count)^2 / expected_count
    end
    dof = n_cchi_bins - 1
    @test chi2 < dof + 5 * sqrt(2 * dof)

    expected_eps_count = n_samples / n_eps_bins
    chi2_eps = 0.0
    for index in eachindex(eps_counts)
        chi2_eps += (eps_counts[index] - expected_eps_count)^2 / expected_eps_count
    end
    dof_eps = n_eps_bins - 1
    @test chi2_eps < dof_eps + 5 * sqrt(2 * dof_eps)
end

@testset "VSS scattering isotropy at alpha=1" begin
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_data(particles_data_path, ["Ar", "He"])

    # alpha = 1 recovers isotropic scattering
    interaction = Interaction(VSS(), species_data[1].mass, species_data[2].mass, 3.25e-10, 0.735, 273.0, 1.0)
    @test interaction.vss_inv_alpha == 1.0

    rng = StableRNG(1234)

    p1 = Particle(1e10, [700.0, -200.0, 350.0], [0.0, 0.0, 0.0])
    p2 = Particle(1e10, [-100.0, 400.0, 50.0], [0.0, 0.0, 0.0])

    collision_data = CollisionData()
    Merzbild.compute_com!(collision_data, interaction, p1, p2)
    Merzbild.compute_g!(collision_data, p1, p2)
    g = collision_data.g

    n_samples = 1000000

    mean_dir = zeros(Float64, 3)
    second_moments = zeros(Float64, 3, 3)

    # for an isotropic distribution the projection onto any fixed axis is uniform in [-1,1]
    tilted_axis = SVector{3,Float64}(1.0, 2.0, -2.0) / 3.0
    n_tilted_bins = 50
    tilted_counts = zeros(Int64, n_tilted_bins)

    for _ in 1:n_samples
        Merzbild.scatter_vss!(rng, collision_data, interaction, p1, p2)

        n_vec = collision_data.g_vec_new / g

        tilted_proj = n_vec[1] * tilted_axis[1] + n_vec[2] * tilted_axis[2] + n_vec[3] * tilted_axis[3]
        tilted_index = min(n_tilted_bins, 1 + floor(Int64, 0.5 * (tilted_proj + 1.0) * n_tilted_bins))
        tilted_counts[tilted_index] += 1

        for i in 1:3
            mean_dir[i] += n_vec[i]
            for j in 1:3
                second_moments[i, j] += n_vec[i] * n_vec[j]
            end
        end
    end

    mean_dir /= n_samples
    second_moments /= n_samples

    @test maximum(abs.(mean_dir)) < 5 * sqrt(1.0 / (3 * n_samples))

    for i in 1:3
        for j in 1:3
            @test abs(second_moments[i, j] - (i == j ? 1.0 / 3.0 : 0.0)) < 5e-3
        end
    end

    expected_tilted_count = n_samples / n_tilted_bins
    chi2_tilted = 0.0
    for index in eachindex(tilted_counts)
        chi2_tilted += (tilted_counts[index] - expected_tilted_count)^2 / expected_tilted_count
    end
    dof_tilted = n_tilted_bins - 1
    @test chi2_tilted < dof_tilted + 5 * sqrt(2 * dof_tilted)
end
