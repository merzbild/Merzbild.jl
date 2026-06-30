@muladd begin

"""
    MaxwellWallBC1D <: AbstractBC

A struct to hold information about a Maxwell reflecting wall (mixture of specular and diffuse scattering) orthogonal to the x-axis. This is species-specific,
as `reflection_velocity_sq` is dependent on the species' mass.

# Fields
* `T`: temperature
* `v`: wall velocity vector
* `accommodation`: accommodation coefficient (a value of 0 corresponds to specular reflection,
    a value of 1 corresponds to purely diffuse reflection) - these are implemented in separate more efficient functions as well
* `reflection_velocities_sq`: pre-computed squared thermal reflection velocities
"""
struct MaxwellWallBC1D <: AbstractBC
    T::Float64
    v::SVector{3,Float64}
    accommodation::Float64  # accommodation coefficient, = 0: specular, = 1: purely diffuse
    reflection_velocity_sq::Float64

    @doc """
        MaxwellWallBC1D(species_data, species, T::Float64, v, accommodation::Float64)

    Construct a `MaxwellWallBC1D` instance for a given species.

    # Fields
    * `species_data`: the list of `Species` data
    * `species`: index of the species for which to create the BC
    * `T`: wall temperature
    * `v`: wall velocity vector
    * `accommodation`: accommodation coefficient
    """
    function MaxwellWallBC1D(species_data, species, T::Float64, v, accommodation::Float64)
        return new(T, v, accommodation,
                   2 * k_B * T / species_data[species].mass)
    end
end

"""
    FullyDiffuseBC1D <: AbstractBC

A struct to hold information about a fully diffuse reflecting wall orthogonal to the x-axis. This is species-specific,
as `reflection_velocity_sq` is dependent on the species' mass.

# Fields
* `T`: temperature
* `v`: wall velocity vector
* `reflection_velocities_sq`: pre-computed squared thermal reflection velocities
"""
struct FullyDiffuseBC1D <: AbstractBC
    T::Float64
    v::SVector{3,Float64}
    reflection_velocity_sq::Float64

    @doc """
        FullyDiffuseBC1D(species_data, species, T::Float64, v)

    Construct a `FullyDiffuseBC1D` instance for a given species.

    # Fields
    * `species_data`: the list of `Species` data
    * `species`: index of the species for which to create the BC
    * `T`: wall temperature
    * `v`: wall velocity vector
    """
    function FullyDiffuseBC1D(species_data, species, T::Float64, v)
        return new(T, v, 2 * k_B * T / species_data[species].mass)
    end
end

"""
    FullyDiffuseBC1D <: AbstractBC

A struct to hold information about a fully specularly reflecting wall orthogonal to the x-axis. Since such reflection simply
flips the sign of the particle's x-velocity, no data is stored in the struct.
"""
struct FullySpecularBC1D <: AbstractBC
end

"""
    specular_reflection_x!(particle, Δt)

Perform specular reflection of a particle in the x direction.

# Positional arguments
* `particle`: the `Particle` instance for which the velocity is reflected
* `Δt`: the time left for the particle to move

# Returns
Unchanged value of `Δt`.
"""
@inline function specular_reflection_x!(particle, Δt)
    @inbounds particle.v = SVector{3, Float64}(-particle.v[1], particle.v[2], particle.v[3])
    return Δt
end

"""
    diffuse_reflection_x!(rng, particle, wall_reflection_v_sq, wall_normal_sign, wall_v, Δt)

Perform diffuse reflection of a particle, assuming the wall is orthogonal to the x axis.

# Positional arguments
* `rng`: the random number generator
* `particle`: the `Particle` instance for which the velocity is reflected
* `wall_reflection_v_sq`: the squared thermal velocity of the species reflected at the wall temperature
* `wall_normal_sign`: sign of the wall normal
* `wall_v`: wall velocity vector
* `Δt`: the time left for the particle to move

# Returns
Unchanged value of `Δt`.
"""
function diffuse_reflection_x!(rng, particle, wall_reflection_v_sq, wall_normal_sign, wall_v, Δt)
    # Note: inlining this function via @inline slows the code down!
    R = max(1e-50, rand(rng))
    v_normal = wall_normal_sign * sqrt(-wall_reflection_v_sq * log(R))
    
    R = max(1e-50, rand(rng))
    v_tang = sqrt(-wall_reflection_v_sq * log(R))

    R = twopi * rand(rng)
    v_tang1, v_tang2 = sincos(R)

    @inbounds particle.v = SVector{3, Float64}(v_normal + wall_v[1], v_tang1 * v_tang + wall_v[2], v_tang2 * v_tang + wall_v[3])
    return Δt
end

"""
    apply_bc!(rng, particle::Particle{D}, wallbc::MaxwellWallBC1D, Δt)

Apply a Maxwell 1D boundary condition.

# Positional arguments
* `rng`: the random number generator
* `particle`: the `Particle` instance for which the velocity is reflected
* `wallbc`: the `MaxwellWallBC1D` instance
* `surface_normal`: the vector of the surface normal
* `Δt`: the time left for the particle to move

# Returns
Unchanged value of `Δt`.
"""
@inline function apply_bc!(rng, particle::Particle{D}, wallbc::MaxwellWallBC1D, surface_normal::SVector{3,Float64}, Δt) where D
    wa = wallbc.accommodation
    if wa == 0.0
        return specular_reflection_x!(particle, Δt)
    elseif wa == 1.0
        @inbounds return diffuse_reflection_x!(rng, particle, wallbc.reflection_velocity_sq, surface_normal[1], wallbc.v, Δt)
    else
        R = rand(rng)
        if R < wa
            @inbounds return diffuse_reflection_x!(rng, particle, wallbc.reflection_velocity_sq, surface_normal[1], wallbc.v, Δt)
        else
            return specular_reflection_x!(particle, Δt)
        end
    end
end

"""
    apply_bc!(rng, particle::Particle{D}, wallbc::FullyDiffuseBC1D, surface_normal::SVector{3,Float64}, Δt)

Apply a fully diffuse 1D boundary condition.

# Positional arguments
* `rng`: the random number generator
* `particle`: the `Particle` instance to which the boundary condition is applied
* `wallbc`: the `FullyDiffuseBC1D` instance
* `surface_normal`: the vector of the surface normal
* `Δt`: the time left for the particle to move

# Returns
Unchanged value of `Δt`.
"""
@inline function apply_bc!(rng, particle::Particle{D}, wallbc::FullyDiffuseBC1D, surface_normal::SVector{3,Float64}, Δt) where D
    @inbounds return diffuse_reflection_x!(rng, particle, wallbc.reflection_velocity_sq, surface_normal[1], wallbc.v, Δt)
end

"""
    apply_bc!(rng, particle::Particle{D}, wallbc::FullySpecularBC1D, Δt)

Apply a fully specular 1D boundary condition.

# Positional arguments
* `rng`: the random number generator
* `particle`: the `Particle` instance to which the boundary condition is applied
* `wallbc`: the `FullySpecularBC1D` instance
* `surface_normal`: the vector of the surface normal
* `Δt`: the time left for the particle to move

# Returns
Unchanged value of `Δt`.
"""
@inline function apply_bc!(rng, particle::Particle{D}, wallbc::FullySpecularBC1D, surface_normal::SVector{3,Float64}, Δt) where D
    return specular_reflection_x!(particle, Δt)
end

"""
    apply_bc_dispatched!(rng, particle::Particle{D}, bcs::Tuple, i::Int, surface_normal::SVector{3,Float64}, Δt::Float64)

Generates a function that dispatches to the correct boundary condition function based on the index `i`. For use in a particle
convection routine.

# Positional arguments:
* `rng`: the random number generator
* `particle`: the `Particle` instance for which the boundary condition is applied
* `bcs`: Tuple of boundary conditions
* `i`: index of boundary condition to apply
* `surface_normal`: the vector of the surface normal
* `Δt`: the time left for the particle to move

# Returns
Potentially changed value of `Δt`.
"""
@generated function apply_bc_dispatched!(rng, particle::Particle{D}, bcs::Tuple, i::Int, surface_normal::SVector{3,Float64}, Δt::Float64) where D
    N_bcs = length(bcs.parameters)
    ex = quote end
    
    # Dynamically build an if-elseif chain at compile time
    for j in 1:N_bcs
        push!(ex.args, :(if i == $j
            # Because $j is a literal constant here, 
            # bcs[$j] is perfectly type-stable to the compiler
            return apply_bc!(rng, particle, bcs[$j], surface_normal, Δt)
        end))
    end
    return ex
end

end