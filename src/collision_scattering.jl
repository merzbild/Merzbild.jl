function scatter_vhs(rng, collision_data, interaction, p1, p2)
    ϕ = 2.0 * π * rand(rng, Float64)
    cphi = cos(ϕ)
    sphi = sin(ϕ)

    ctheta = 2.0 * rand(rng, Float64) - 1.0
    stheta = sqrt(1.0 - ctheta^2)

    collision_data.g_vec_new[1] = collision_data.g * stheta * cphi
    collision_data.g_vec_new[2] = collision_data.g * stheta * sphi
    collision_data.g_vec_new[3] = collision_data.g * ctheta

    p1.v = collision_data.v_com .+ interaction.μ2 * collision_data.g_vec_new
    p2.v = collision_data.v_com .- interaction.μ1 * collision_data.g_vec_new
end