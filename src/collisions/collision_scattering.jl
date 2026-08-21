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
    sphi, cphi = sincos(ϕ)

    ctheta = 2.0 * rand(rng, Float64) - 1.0
    stheta = sqrt(1.0 - ctheta^2)

    g_vec_new = collision_data.g * SVector{3,Float64}(stheta * cphi, stheta * sphi, ctheta)

    v_com = collision_data.v_com

    p1.v = v_com + interaction.μ2 * g_vec_new
    p2.v = v_com - interaction.μ1 * g_vec_new
    collision_data.g_vec_new = g_vec_new
end

"""
    scatter_vss!(rng, collision_data, interaction, p1, p2)

Scatter two particles using VSS (anisotropic) scattering.
The cosine of the deflection angle ``\\chi`` (the angle between the pre- and post-collisional
relative velocities) is sampled as ``\\cos\\chi = 2 R^{1/\\alpha} - 1``, where ``R`` is a uniformly
distributed random number and ``\\alpha`` is the VSS exponent of the interaction; the azimuthal
angle ``\\varepsilon`` is sampled uniformly in ``[0, 2\\pi)``. The post-collisional relative velocity
is then computed by rotating the pre-collisional one, which requires the pre-collisional relative
velocity vector stored in `collision_data` (as computed by [`Merzbild.compute_g!`](@ref)).

The power ``R^{1/\\alpha}`` is evaluated as `exp2(vss_inv_alpha * log2(R))` for the same reason
as in [`Merzbild.sigma_vhs`](@ref).

If the relative velocity is nearly aligned with the x-axis, the rotation is performed in the
``(y,z)`` plane instead, as the general expressions become ill-conditioned.

# Positional arguments
* `rng`: the random number generator
* `collision_data`: the `CollisionData` instance which stores the center-of-mass velocity,
    the pre-collisional relative velocity of the particles (vector and magnitude), and to which
    the new post-collisional velocity will be written
* `interaction`: the `Interaction` instance for the colliding particles
* `p1`: the first colliding particle
* `p2`: the second colliding particle

# References
* G.A. Bird, Eq. (2.22), Molecular gas dynamics and the direct simulation of gas flows,
    [Clarendon Press, Oxford, 1994](https://doi.org/10.1093/oso/9780198561958.001.0001).
* K. Koura, H. Matsumoto, Variable soft sphere molecular model for inverse-power-law or Lennard-Jones
    potential. [Phys. Fluids A, 1991](https://doi.org/10.1063/1.857792).
"""
@inline function scatter_vss!(rng, collision_data, interaction, p1, p2)
    cchi = 2.0 * exp2(interaction.vss_inv_alpha * log2(rand(rng, Float64))) - 1.0
    schi = sqrt(1.0 - cchi^2)

    ε = twopi * rand(rng, Float64)
    seps, ceps = sincos(ε)

    g = collision_data.g
    g_vec = collision_data.g_vec
    gyz = sqrt(g_vec[2]^2 + g_vec[3]^2)

    if (gyz > 1e-5 * g)
        inv_gyz = 1.0 / gyz
        g_vec_new = SVector{3,Float64}(cchi * g_vec[1] + schi * seps * gyz,
                                       cchi * g_vec[2] + schi * (g * g_vec[3] * ceps - g_vec[1] * g_vec[2] * seps) * inv_gyz,
                                       cchi * g_vec[3] - schi * (g * g_vec[2] * ceps + g_vec[1] * g_vec[3] * seps) * inv_gyz)
    else
        g_vec_new = SVector{3,Float64}(cchi * g_vec[1], schi * g * ceps, schi * g * seps)
    end

    v_com = collision_data.v_com

    p1.v = v_com + interaction.μ2 * g_vec_new
    p2.v = v_com - interaction.μ1 * g_vec_new
    collision_data.g_vec_new = g_vec_new
end

"""
    scatter!(rng, model, collision_data, interaction, p1, p2)

Scatter two particles using the scattering law of the elastic scattering `model`:
isotropic scattering ([`Merzbild.scatter_vhs!`](@ref)) for the [`VHS`](@ref) model,
VSS scattering ([`Merzbild.scatter_vss!`](@ref)) for the [`VSS`](@ref) model.

# Positional arguments
* `rng`: the random number generator
* `model`: the `AbstractScatteringModel` singleton tag of the scattering model
* `collision_data`: the `CollisionData` instance which stores the center-of-mass velocity and
    the pre-collisional relative velocity of the particles, and to which
    the new post-collisional velocity will be written
* `interaction`: the `Interaction` instance for the colliding particles
* `p1`: the first colliding particle
* `p2`: the second colliding particle
"""
@inline scatter!(rng, ::VHS, collision_data, interaction, p1, p2) = scatter_vhs!(rng, collision_data, interaction, p1, p2)
@inline scatter!(rng, ::VSS, collision_data, interaction, p1, p2) = scatter_vss!(rng, collision_data, interaction, p1, p2)

"""
    scatter_electron_vhs!(rng, particle_electron, g_new)

Scatter an electron using VHS (isotropic) scattering and
re-scale its relative velocity to `g_new`. This **DOES NOT** add
the velocity of the center of mass to the electron.

# Positional arguments
* `rng`: the random number generator
* `particle_electron`: the electron particle to scatter off of the neutral particle
* `g_new`: the magnitude of the post-collisional relative velocity
"""
@inline function scatter_electron_vhs!(rng, particle_electron, g_new)
    ϕ = twopi * rand(rng, Float64)
    sphi, cphi = sincos(ϕ)

    ctheta = 2.0 * rand(rng, Float64) - 1.0
    stheta = sqrt(1.0 - ctheta^2)

    particle_electron.v = g_new * SVector{3, Float64}(ctheta, stheta * cphi, stheta * sphi)
end

"""
    scatter_ionization_electrons_and_ion!(rng, collision_data, p_e1, p_e2, p_ion, mass_ratio)

Scatter electrons and ion after an ionization reaction using VHS (isotropic) scattering.

# Positional arguments
* `rng`: the random number generator
* `collision_data`: the `CollisionData` instance to which stores the
    center-of-mass velocity the magnitude of the pre-collisional relative velocity
    of the electron and the neutral, and the post-collisional magnitudes of the velocities
    of the electrons
* `p_e1`: the first electron particle to scatter off of the neutral
* `p_e2`: the second electron particle to scatter off of the neutral
* `p_ion`: the ion produced in the ionization reaction
* `mass_ratio`: ratio of the electron mass to the ion mass

# References
* K. Nanbu, Eqns. (47)-(53b), [IEEE Trans. Plasma. Sci., 2000](https://doi.org/10.1109/27.887765)
"""
function scatter_ionization_electrons_and_ion!(rng, collision_data, p_e1, p_e2, p_ion, mass_ratio)
    scatter_electron_vhs!(rng, p_e1, collision_data.g_new_1)
    scatter_electron_vhs!(rng, p_e2, collision_data.g_new_2)
    p_ion.v = -mass_ratio * (p_e1.v + p_e2.v) + collision_data.v_com
    
    p_e1.v = p_e1.v + collision_data.v_com
    p_e2.v = p_e2.v + collision_data.v_com
end

end