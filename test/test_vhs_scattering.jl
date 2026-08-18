@testset "VHS scattering isotropy" begin
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    interaction_data_path = joinpath(@__DIR__, "..", "data", "vhs.toml")

    species_data = load_species_data(particles_data_path, ["Ar", "He"])
    interaction_data = load_interaction_data(interaction_data_path, species_data)
    interaction = interaction_data[1, 2]

    rng = StableRNG(1234)

    p1 = Particle(1e10, [700.0, -200.0, 350.0], [0.0, 0.0, 0.0])
    p2 = Particle(1e10, [-100.0, 400.0, 50.0], [0.0, 0.0, 0.0])

    collision_data = CollisionData()
    Merzbild.compute_com!(collision_data, interaction, p1, p2)
    Merzbild.compute_g!(collision_data, p1, p2)

    v_com = collision_data.v_com
    g = collision_data.g

    collision_data_check = CollisionData()

    n_samples = 1000000

    # bins uniform in cos(theta) and phi are of equal area on the unit sphere
    n_ctheta_bins = 40
    n_phi_bins = 40
    counts = zeros(Int64, n_ctheta_bins, n_phi_bins)

    # for an isotropic distribution the projection onto any fixed axis is uniform in [-1,1],
    # not just onto the coordinate axes
    # pick an axis and normalize its length
    tilted_axis = SVector{3,Float64}(1.0, 2.0, -2.0) / 3.0
    n_tilted_bins = 50
    tilted_counts = zeros(Int64, n_tilted_bins)

    mean_dir = zeros(Float64, 3)
    second_moments = zeros(Float64, 3, 3)

    max_g_err = 0.0
    max_v_com_err = 0.0

    for _ in 1:n_samples
        Merzbild.scatter_vhs!(rng, collision_data, interaction, p1, p2)

        n_vec = collision_data.g_vec_new / g

        ctheta_index = min(n_ctheta_bins, 1 + floor(Int64, 0.5 * (n_vec[3] + 1.0) * n_ctheta_bins))
        phi_index = min(n_phi_bins, 1 + floor(Int64, (atan(n_vec[2], n_vec[1]) + pi) * n_phi_bins / (2 * pi)))
        counts[ctheta_index, phi_index] += 1

        tilted_proj = n_vec[1] * tilted_axis[1] + n_vec[2] * tilted_axis[2] + n_vec[3] * tilted_axis[3]
        tilted_index = min(n_tilted_bins, 1 + floor(Int64, 0.5 * (tilted_proj + 1.0) * n_tilted_bins))
        tilted_counts[tilted_index] += 1

        for i in 1:3
            mean_dir[i] += n_vec[i]
            for j in 1:3
                second_moments[i, j] += n_vec[i] * n_vec[j]
            end
        end

        # scattering is elastic: |g| and v_com are unchanged
        Merzbild.compute_g!(collision_data_check, p1, p2)
        Merzbild.compute_com!(collision_data_check, interaction, p1, p2)
        max_g_err = max(max_g_err, abs(collision_data_check.g - g) / g)
        max_v_com_err = max(max_v_com_err, maximum(abs.(collision_data_check.v_com - v_com)) / g)
    end

    mean_dir /= n_samples
    second_moments /= n_samples

    @test max_g_err < 1e-13
    @test max_v_com_err < 1e-13

    # <n> = 0, standard deviation of a single component is sqrt(1/(3 n_samples))
    @test maximum(abs.(mean_dir)) < 5 * sqrt(1.0 / (3 * n_samples))

    # <n_i n_j> = delta_ij / 3
    for i in 1:3
        for j in 1:3
            @test abs(second_moments[i, j] - (i == j ? 1.0 / 3.0 : 0.0)) < 5e-3
        end
    end

    expected_count = n_samples / (n_ctheta_bins * n_phi_bins)
    chi2 = 0.0
    for index in eachindex(counts)
        chi2 += (counts[index] - expected_count)^2 / expected_count
    end
    dof = n_ctheta_bins * n_phi_bins - 1
    @test chi2 < dof + 5 * sqrt(2 * dof)

    expected_tilted_count = n_samples / n_tilted_bins
    chi2_tilted = 0.0
    for index in eachindex(tilted_counts)
        chi2_tilted += (tilted_counts[index] - expected_tilted_count)^2 / expected_tilted_count
    end
    dof_tilted = n_tilted_bins - 1
    @test chi2_tilted < dof_tilted + 5 * sqrt(2 * dof_tilted)
end
