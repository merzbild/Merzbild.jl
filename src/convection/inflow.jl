using SpecialFunctions: erf

@muladd begin

"""
    AbstractInflowSampler

Abstract type for the samplers of the Maxwellian inflow distribution used by the inflow boundary
conditions.

Each concrete sampler is a zero-field singleton type ([`InflowNoDrift`](@ref),
[`InflowEnvelope4`](@ref)) used as a compile-time tag, and [`Merzbild.sample_inflow_z`](@ref) is
dispatched on it. The tag is stored in an [`InflowBC1D`](@ref) instances.

Which sampler is applicable is determined by the speed ratio ``a``, see [`InflowBC1D`](@ref).
"""
abstract type AbstractInflowSampler end

"""
    InflowNoDrift

Singleton type tag for the sampler of the Maxwellian inflow distribution in the absence of a
streaming velocity normal to the surface (``a = 0``), for which the distribution reduces to a
Rayleigh distribution and is sampled by direct inversion.
"""
struct InflowNoDrift <: AbstractInflowSampler end

"""
    InflowEnvelope4

Singleton type tag for the sampler of the Maxwellian inflow distribution based on Envelope 4 of
Garcia and Wagner, an exact acceptance-rejection scheme efficient over a wide range of speed
ratios. Only the ``a > 0`` branch is implemented, i.e. the streaming velocity may not point
out of the domain; ``a = 0`` is handled by [`InflowNoDrift`](@ref).

# References
* A.L. Garcia, W. Wagner, Generation of the Maxwellian inflow distribution.
    [J. Comput. Phys., 2006](https://doi.org/10.1016/j.jcp.2006.01.025), section 5.2.4 and Table 2.
"""
struct InflowEnvelope4 <: AbstractInflowSampler end

"""
    maxwellian_inflow_normalization(a)

Compute the normalization factor ``m(a) = \\exp(-a^2) + a \\sqrt{\\pi} [1 + \\mathrm{erf}(a)]``
of the Maxwellian inflow distribution.

# Positional arguments
* `a`: the speed ratio (normal streaming velocity divided by thermal velocity)

# Returns
The value of ``m(a)``.
"""
@inline function maxwellian_inflow_normalization(a)
    return exp(-a * a) + a * sqrt_pi * (1.0 + erf(a))
end

"""
    compute_inflow_speed_ratio(m, T, v, normal)

Compute the speed ratio ``a = dot(v, n) / \\sqrt{2kT/m}`` of the streaming velocity along the inward
surface normal.

# Positional arguments
* `m`: the molecular mass of the species
* `T`: the temperature of the reservoir
* `v`: the streaming velocity vector of the reservoir
* `normal`: the inward surface normal

# Returns
A tuple of the thermal velocity and the speed ratio.
"""
function compute_inflow_speed_ratio(m, T, v, normal)
    v_thermal = compute_thermal_velocity(m, T)
    @inbounds v_drift_normal = v[1] * normal[1] + v[2] * normal[2] + v[3] * normal[3]
    return v_thermal, v_drift_normal / v_thermal
end

"""
    compute_inflow_number_flux(species, species_data, T, v, normal, ndens)

Compute the number flux (number of particles crossing a unit surface area per unit time) of an
equilibrium reservoir into the domain,
``\\Gamma = n \\sqrt{2kT/m} \\, m(a) / (2\\sqrt{\\pi})``,
where ``m(a)`` is given by [`Merzbild.maxwellian_inflow_normalization`](@ref) and ``a`` is the
speed ratio of the streaming velocity along the inward surface normal (i.e., normal streaming velocity divided by thermal velocity).
For a reservoir at rest this reduces to the effusion flux ``n \\sqrt{kT/(2 \\pi m)}``.

# Positional arguments
* `species`: index of the species for which the flux is computed
* `species_data`: the list of `Species` data
* `T`: the temperature of the reservoir
* `v`: the streaming velocity vector of the reservoir
* `normal`: the inward surface normal
* `ndens`: the number density of the reservoir

# Returns
The number flux.

# References
* A.L. Garcia, W. Wagner, Generation of the Maxwellian inflow distribution.
    [J. Comput. Phys., 2006](https://doi.org/10.1016/j.jcp.2006.01.025), section 7.
"""
function compute_inflow_number_flux(species, species_data, T, v, normal, ndens)
    v_thermal, a = compute_inflow_speed_ratio(species_data[species].mass, T, v, normal)
    return ndens * v_thermal * maxwellian_inflow_normalization(a) / (2 * sqrt_pi)
end

"""
    check_inflow_sampler(sampler, a)

Check that a sampler of the Maxwellian inflow distribution is applicable for a given speed ratio.

# Positional arguments
* `sampler`: the `AbstractInflowSampler` tag
* `a`: the speed ratio

# Throws
`ArgumentError` if the sampler cannot be used for the speed ratio in question.
"""
function check_inflow_sampler(::InflowNoDrift, a)
    if a != 0.0
        throw(ArgumentError("InflowNoDrift requires a zero speed ratio, got a = " * string(a)))
    end
    return nothing
end

function check_inflow_sampler(::InflowEnvelope4, a)
    if a < 0.0
        throw(ArgumentError("InflowEnvelope4 is implemented for positive speed ratios only, got a = " * string(a) *
                            "; sampling a streaming velocity pointing out of the domain requires Envelope 2"))
    elseif !(0.0 < a < Inf)
        throw(ArgumentError("InflowEnvelope4 requires a positive finite speed ratio (InflowNoDrift is to be used for a = 0), got a = " *
                            string(a)))
    end
    return nothing
end

"""
    InflowBC1D{S} <: AbstractBC

A struct to hold information about an inflow surface orthogonal to the x-axis, across which an
equilibrium reservoir with a given number density, temperature and streaming velocity flows into
the domain. This is species-specific, as the thermal velocity depends on the species' mass.

Particles of the domain reaching such a surface leave the domain and are deleted, as for an
[`OutflowBC`](@ref).

The type parameter `S` is the [`AbstractInflowSampler`](@ref) tag used to sample the velocity
component normal to the surface. It is selected from the speed ratio `a` by the constructor and is
therefore known only at run-time: instances are to be created during set-up and not in a hot loop.

# Fields
* `sampler`: the `AbstractInflowSampler` tag
* `T`: temperature of the reservoir
* `v`: streaming velocity vector of the reservoir
* `normal`: the inward surface normal
* `ndens`: number density of the reservoir
* `ndens_flux`: pre-computed number flux across the surface
* `v_thermal`: pre-computed thermal velocity of the species at the reservoir temperature
* `v_tangential_scale`: pre-computed standard deviation of the tangential velocity components
* `v_drift_normal`: pre-computed projection of the streaming velocity on the inward surface normal
* `a`: pre-computed speed ratio
* `inv_a`: pre-computed inverse of the speed ratio (``0`` if the speed ratio is zero)
* `envelope4_prob`: pre-computed mixture probability ``1 / (2 a \\sqrt{\\pi} + 1)`` of Envelope 4
"""
struct InflowBC1D{S<:AbstractInflowSampler} <: AbstractBC
    sampler::S
    T::Float64
    v::SVector{3,Float64}
    normal::SVector{3,Float64}
    ndens::Float64
    ndens_flux::Float64
    v_thermal::Float64
    v_tangential_scale::Float64
    v_drift_normal::Float64
    a::Float64
    inv_a::Float64
    envelope4_prob::Float64

    @doc """
        InflowBC1D(sampler::AbstractInflowSampler, species, species_data, T::Float64, v, normal, ndens::Float64)

    Construct an `InflowBC1D` instance for a given species using an explicitly chosen sampler of the
    Maxwellian inflow distribution.

    # Positional arguments
    * `sampler`: the `AbstractInflowSampler` tag
    * `species`: index of the species for which to create the BC
    * `species_data`: the list of `Species` data
    * `T`: temperature of the reservoir
    * `v`: streaming velocity vector of the reservoir
    * `normal`: the inward surface normal, `(1, 0, 0)` or `(-1, 0, 0)`
    * `ndens`: number density of the reservoir

    # Throws
    `ArgumentError` if `normal` is not `(1, 0, 0)` or `(-1, 0, 0)`, or if the sampler is not applicable
    for the resulting speed ratio.
    """
    function InflowBC1D(sampler::AbstractInflowSampler, species, species_data, T::Float64, v, normal, ndens::Float64)
        v = SVector{3,Float64}(v)
        normal = SVector{3,Float64}(normal)

        if !(abs(normal[1]) == 1.0 && normal[2] == 0.0 && normal[3] == 0.0)
            throw(ArgumentError("InflowBC1D requires an inward surface normal of (1, 0, 0) or (-1, 0, 0), got " *
                                string(normal)))
        end

        v_thermal, a = compute_inflow_speed_ratio(species_data[species].mass, T, v, normal)
        check_inflow_sampler(sampler, a)

        inv_a = a == 0.0 ? 0.0 : 1.0 / a

        return new{typeof(sampler)}(sampler, T, v, normal, ndens,
                                    compute_inflow_number_flux(species, species_data, T, v, normal, ndens),
                                    v_thermal, inv_sqrt_2 * v_thermal, a * v_thermal, a, inv_a,
                                    1.0 / (2 * a * sqrt_pi + 1.0))
    end

    @doc """
        InflowBC1D(species, species_data, T::Float64, v, normal, ndens::Float64)

    Construct an `InflowBC1D` instance for a given species, choosing the sampler of the Maxwellian
    inflow distribution based on the speed ratio: [`InflowNoDrift`](@ref) if the streaming velocity has
    no component along the surface normal, [`InflowEnvelope4`](@ref) otherwise.
    Since the type of the returned instance depends on the value of `v`, this constructor is not type-stable:
    in performance-critical code, either use the constructor with an explicitly specified sampler,
    or pass the instance (or the `Tuple` of boundary conditions containing it) through a function barrier.

    # Positional arguments
    * `species`: index of the species for which to create the BC
    * `species_data`: the list of `Species` data
    * `T`: temperature of the reservoir
    * `v`: streaming velocity vector of the reservoir
    * `normal`: the inward surface normal, `(1, 0, 0)` or `(-1, 0, 0)`
    * `ndens`: number density of the reservoir

    # Throws
    `ArgumentError` if `normal` is not `(1, 0, 0)` or `(-1, 0, 0)`, or if the streaming velocity points
    out of the domain.
    """
    function InflowBC1D(species, species_data, T::Float64, v, normal, ndens::Float64)
        v = SVector{3,Float64}(v)
        normal = SVector{3,Float64}(normal)
        _, a = compute_inflow_speed_ratio(species_data[species].mass, T, v, normal)

        if a == 0.0
            return InflowBC1D(InflowNoDrift(), species, species_data, T, v, normal, ndens)
        else
            return InflowBC1D(InflowEnvelope4(), species, species_data, T, v, normal, ndens)
        end
    end
end

"""
    sample_inflow_z(rng, sampler, inflowbc::InflowBC1D)

Sample the scaled velocity component ``z`` from the Maxwellian inflow distribution
``p_a(z) = 2 m(a)^{-1} (a - z) \\exp(-z^2)``, ``z \\in (-\\infty, a)``, where ``a`` is the speed
ratio of the boundary condition. The velocity component along the inward surface normal is
recovered as ``(a - z) \\sqrt{2kT/m}``.

The [`InflowNoDrift`](@ref) sampler inverts the ``a = 0`` distribution directly; the
[`InflowEnvelope4`](@ref) sampler uses acceptance-rejection with Envelope 4.

# Positional arguments
* `rng`: the random number generator
* `sampler`: the `AbstractInflowSampler` tag
* `inflowbc`: the `InflowBC1D` instance

# Returns
The sampled value of ``z``.

# References
* A.L. Garcia, W. Wagner, Generation of the Maxwellian inflow distribution.
    [J. Comput. Phys., 2006](https://doi.org/10.1016/j.jcp.2006.01.025), equation (3.6) and Table 2.
"""
@inline function sample_inflow_z(rng, ::InflowNoDrift, inflowbc::InflowBC1D)
    return -sqrt(randexp(rng))
end

function sample_inflow_z(rng, ::InflowEnvelope4, inflowbc::InflowBC1D)
    a = inflowbc.a
    inv_a = inflowbc.inv_a

    while true
        if inflowbc.envelope4_prob > rand(rng)
            z = -sqrt(randexp(rng))
        else
            z = inv_sqrt_2 * randn(rng)
        end

        if (a - z) * inv_a > rand(rng)
            return z
        end
    end
end

"""
    sample_inflow_velocity(rng, inflowbc::InflowBC1D)

Sample the velocity of a particle entering the domain across an inflow surface orthogonal to the
x-axis. The component normal to the surface is drawn from the Maxwellian inflow distribution and
always points into the domain; the tangential components are Gaussians centered on the tangential
streaming velocity.

# Positional arguments
* `rng`: the random number generator
* `inflowbc`: the `InflowBC1D` instance

# Returns
The sampled velocity vector.
"""
@inline function sample_inflow_velocity(rng, inflowbc::InflowBC1D)
    z = sample_inflow_z(rng, inflowbc.sampler, inflowbc)

    v_normal = inflowbc.v_drift_normal - z * inflowbc.v_thermal
    v_tangential_scale = inflowbc.v_tangential_scale

    @inbounds return SVector{3,Float64}(inflowbc.normal[1] * v_normal,
                                        inflowbc.v[2] + v_tangential_scale * randn(rng),
                                        inflowbc.v[3] + v_tangential_scale * randn(rng))
end

end
