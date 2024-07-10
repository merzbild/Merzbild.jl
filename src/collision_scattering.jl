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

function scatter_electron_vhs!(rng, collision_data, particles_electron, g_new, i)
    ϕ = twopi * rand(rng, Float64)
    cphi = cos(ϕ)
    sphi = sin(ϕ)

    ctheta = 2.0 * rand(rng, Float64) - 1.0
    stheta = sqrt(1.0 - ctheta^2)

    ratio = g_new / collision_data.g
    ur = particles_electron[i].v[1] * ratio
    vr = particles_electron[i].v[2] * ratio
    wr = particles_electron[i].v[3] * ratio

    denominator = sqrt(vr*vr + wr*wr)    

    if (denominator < eps()*g_new)
        particles_electron[i].v = SVector{3, Float64}(ctheta * ur + stheta * sphi * denominator,
                                                      ctheta * vr + stheta * (g_new * wr * cphi - ur * vr * sphi) / denominator,
                                                      ctheta * wr - stheta * (g_new * vr * cphi + ur * wr * sphi) / denominator) + collision_data.v_com
    else
        particles_electron[i].v = SVector{3, Float64}(ctheta * ur, stheta * cphi * ur, stheta * sphi * ur) + collision_data.v_com
    end
end

function scatter_ionization_electrons!(rng, collision_data, particles_electron, i1, i2)
    scatter_electron_vhs!(rng, collision_data, particles_electron, collision_data.g_new_1, i1)
    scatter_electron_vhs!(rng, collision_data, particles_electron, collision_data.g_new_2, i2)
end