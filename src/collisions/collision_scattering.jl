@muladd begin

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
    scatter_electron_vhs!(rng, collision_data, particles_electron, g_new)

Scatter an electron using VHS (isotropic) scattering and
re-scale its relative velocity to `g_new`. This **DOES NOT** add
the velocity of the center of mass to the electron.

# Positional arguments
* `rng`: the random number generator
* `particles_electron`: the electron particle to scatter off of the neutral particle
* `g_new`: the magnitude of the post-collisional relative velocity
"""
@inline function scatter_electron_vhs!(rng, particles_electron, g_new)
    ϕ = twopi * rand(rng, Float64)
    cphi = cos(ϕ)
    sphi = sin(ϕ)

    ctheta = 2.0 * rand(rng, Float64) - 1.0
    stheta = sqrt(1.0 - ctheta^2)

    particles_electron.v = SVector{3, Float64}(ctheta * g_new, stheta * cphi * g_new, stheta * sphi * g_new)
end

"""
    scatter_ionization_electrons_and_ion!(rng, collision_data, particles_electron, particles_ion, i1, i2, k1, mass_ratio)

Scatter electrons and ion after an ionization reaction using VHS (isotropic) scattering.

# Positional arguments
* `rng`: the random number generator
* `collision_data`: the `CollisionData` instance to which stores the
    center-of-mass velocity the magnitude of the pre-collisional relative velocity
    of the electron and the neutral, and the post-collisional magnitudes of the velocities
    of the electrons
* `particles_electron`: the vector of electron particles
* `particles_electron`: the vector of ion particles
* `i1`: the index of the first electron particle to scatter off of the neutral
* `i2`: the index of the second electron particle to scatter off of the neutral
* `k1`: the index of the ion produced in the ionization reaction
* `mass_ratio`: ratio of the electron mass to the ion mass

# References
* K. Nanbu, Eqns. (47)-(53b), [IEEE Trans. Plasma. Sci., 2000](https://doi.org/10.1109/27.887765)
"""
function scatter_ionization_electrons_and_ion!(rng, collision_data, particles_electron, particles_ion, i1, i2, k1, mass_ratio)
    scatter_electron_vhs!(rng, particles_electron[i1], collision_data.g_new_1)
    scatter_electron_vhs!(rng, particles_electron[i2], collision_data.g_new_2)
    particles_ion[k1].v = -mass_ratio * (particles_electron[i1].v + particles_electron[i2].v) + collision_data.v_com
    
    particles_electron[i1].v = particles_electron[i1].v + collision_data.v_com
    particles_electron[i2].v = particles_electron[i2].v + collision_data.v_com
end

end