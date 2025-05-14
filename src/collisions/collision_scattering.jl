"""
    scatter_vhs!(rng, collision_data, interaction, p1, p2)

Scatter two particles using VHS (isotropic) scattering.

# Positional arguments
* `rng`: the random number generator
* `collision_data`: the `CollisionData` instance to which stores the
    center-of-mass velocity, the magnitude of the pre-collisional relative velocity of the particles and to which
    the new post-collisional velocity will be written
* `interaction`: the `Interaction` instance for the colliding particles
* `p1`: the first colliding particle
* `p2`: the second colliding particle
"""
@inline function scatter_vhs!(rng, collision_data, interaction, p1, p2)
    ϕ = twopi * rand(rng, Float64)
    cphi = cos(ϕ)
    sphi = sin(ϕ)

    ctheta = 2.0 * rand(rng, Float64) - 1.0
    stheta = sqrt(1.0 - ctheta^2)

    collision_data.g_vec_new = collision_data.g * SVector{3,Float64}(stheta * cphi, stheta * sphi, ctheta)

    p1.v = collision_data.v_com + interaction.μ2 * collision_data.g_vec_new
    p2.v = collision_data.v_com - interaction.μ1 * collision_data.g_vec_new
end

"""
    scatter_vhs!(rng, collision_data, interaction, p1, p2)

Scatter an electron off of a neutral particle using VHS (isotropic) scattering and
re-scale its relative velocity to `g_new`.

# Positional arguments
* `rng`: the random number generator
* `collision_data`: the `CollisionData` instance to which stores the
    center-of-mass velocity and the magnitude of the pre-collisional relative velocity
    of the electron and the neutral
* `particles_electron`: the electron particle to scatter off of the neutral particle
* `g_new`: the magnitude of the post-collisional relative velocity
"""
@inline function scatter_electron_vhs!(rng, collision_data, particles_electron, g_new)
    ϕ = twopi * rand(rng, Float64)
    cphi = cos(ϕ)
    sphi = sin(ϕ)

    ctheta = 2.0 * rand(rng, Float64) - 1.0
    stheta = sqrt(1.0 - ctheta^2)

    ratio = g_new / collision_data.g
    ur = particles_electron.v[1] * ratio
    vr = particles_electron.v[2] * ratio
    wr = particles_electron.v[3] * ratio

    denominator = sqrt(vr*vr + wr*wr)    

    if (denominator < eps()*g_new)
        particles_electron.v = SVector{3, Float64}(ctheta * ur + stheta * sphi * denominator,
                                                      ctheta * vr + stheta * (g_new * wr * cphi - ur * vr * sphi) / denominator,
                                                      ctheta * wr - stheta * (g_new * vr * cphi + ur * wr * sphi) / denominator) + collision_data.v_com
    else
        particles_electron.v = SVector{3, Float64}(ctheta * ur, stheta * cphi * ur, stheta * sphi * ur) + collision_data.v_com
    end
end

"""
    scatter_ionization_electrons!(rng, collision_data, particles_electron, i1, i2)

Scatter two electrons off of a neutral using VHS (isotropic) scattering.

# Positional arguments
* `rng`: the random number generator
* `collision_data`: the `CollisionData` instance to which stores the
    center-of-mass velocity the magnitude of the pre-collisional relative velocity
    of the electron and the neutral, and the post-collisional magnitudes of the velocities
    of the electrons
* `particles_electron`: the vector of electron particles to scatter off of the neutral particle
* `i1`: the index of the first electron particle to scatter off of the neutral
* `i2`: the index of the second electron particle to scatter off of the neutral
"""
function scatter_ionization_electrons!(rng, collision_data, particles_electron, i1, i2)
    scatter_electron_vhs!(rng, collision_data, particles_electron[i1], collision_data.g_new_1)
    scatter_electron_vhs!(rng, collision_data, particles_electron[i2], collision_data.g_new_2)
end