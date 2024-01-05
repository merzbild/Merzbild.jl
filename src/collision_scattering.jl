function scatter_vhs(rng, collision_data, interaction, p1, p2)
    ϕ = twopi * rand(rng, Float64)
    cphi = cos(ϕ)
    sphi = sin(ϕ)

    ctheta = 2.0 * rand(rng, Float64) - 1.0
    stheta = sqrt(1.0 - ctheta^2)

    collision_data.g_vec_new = collision_data.g * SVector{3,Float64}(stheta * cphi, stheta * sphi, ctheta)

    p1.v = collision_data.v_com + interaction.μ2 * collision_data.g_vec_new
    p2.v = collision_data.v_com - interaction.μ1 * collision_data.g_vec_new
end