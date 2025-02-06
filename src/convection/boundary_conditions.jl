"""
    MaxwellWallBC

A struct to hold information about a diffuse reflecting wall.

# Fields
* `T`: temperature
* `v`: wall velocity vector
* `accommodation`: accommodation coefficient (a value of 0 corresponds to specular reflection,
    a value of 1 corresponds to purely diffuse reflection)
"""
struct MaxwellWallBC
    T::Float64
    v::SVector{3,Float64}
    accommodation::Float64  # accommodation coefficient, = 0: specular, = 1: purely diffuse
end

"""
    MaxwellWalls1D

A struct to hold information about all diffuse reflecting walls in the simulation.

# Fields
* `boundaries`: a vector of `MaxwellWallBC` walls
* `reflection_velocities_sq`: array of pre-computed squared thermal reflection velocities (n_walls x n_species)
"""
struct MaxwellWalls1D
    boundaries::Vector{MaxwellWallBC}
    reflection_velocities_sq::Array{Float64,2} # per-wall, per-species
end

@doc """
    MaxwellWalls1D(species_data, T_l::Float64, T_r::Float64, vy_l::Float64, vy_r::Float64, accomodation_l::Float64, accomodation_r::Float64)

Create a `MaxwellWalls1D` struct for a 1-D simulation with 2 walls ("left" and "right") with a velocity in
    the y-direction only.

# Positional arguments
* `species_data`: vector of `Species` data
* `T_l`: temperature of left wall
* `T_r`: temperature of right wall
* `vy_l`: y-velocity of left wall
* `vy_r`: y-velocity of right wall
* `accomodation_l`: accommodation coefficient of left wall
* `accomodation_r`: accommodation coefficient of right wall
"""
MaxwellWalls1D(species_data, T_l::Float64, T_r::Float64, vy_l::Float64, vy_r::Float64, accomodation_l::Float64, accomodation_r::Float64) =
MaxwellWalls1D([MaxwellWallBC(T_l, [0.0, vy_l, 0.0], accomodation_l),
               MaxwellWallBC(T_r, [0.0, vy_r, 0.0], accomodation_r)],
               vcat([2 * k_B * T_l / sd.mass for sd in species_data]',
                    [2 * k_B * T_r / sd.mass for sd in species_data]'))

"""
    specular_reflection_x!(particle)

Perform specular reflection of a particle in the x direction.

# Positional arguments
* `particle`: the `Particle` instance for which the velocity is reflected
"""
function specular_reflection_x!(particle)
    particle.v = SVector{3, Float64}(-particle.v[1], particle.v[2], particle.v[3])
end

"""
    specular_reflection_x!(particle)

Perform diffuse reflection of a particle, assuming the wall is orthogonal to the x axis.

# Positional arguments
* `rng`: the random number generator
* `particle`: the `Particle` instance for which the velocity is reflected
* `wall_reflection_v_sq`: the squared thermal velocity of the species reflected at the wall temperature
* `wall_normal_sign`: sign of the wall normal
* `wall_v`: wall velocity vector
"""
function diffuse_reflection_x!(rng, particle, wall_reflection_v_sq, wall_normal_sign, wall_v)
    R = max(1e-50, rand(rng))

    v_normal = wall_normal_sign * sqrt(-wall_reflection_v_sq * log(R))
    
    R = max(1e-50, rand(rng))
    v_tang = sqrt(-wall_reflection_v_sq * log(R))

    R = rand(rng)
    v_tang1 = sin(twopi * R) * v_tang 
    v_tang2 = cos(twopi * R) * v_tang

    particle.v = SVector{3, Float64}(v_normal + wall_v[1], v_tang1 + wall_v[2], v_tang2 + wall_v[3])
end

"""
    reflect_particle_x!(rng, particle, wall_reflection_v_sq, wall_normal_sign, wall_v, wall_accommodation)

Reflect particle from a Maxwell wall orthogonal to the x axis.

# Positional arguments
* `rng`: the random number generator
* `particle`: the `Particle` instance for which the velocity is reflected
* `wall_reflection_v_sq`: the squared thermal velocity of the species reflected at the wall temperature
* `wall_normal_sign`: sign of the wall normal
* `wall_v`: wall velocity vector
* `wall_accommodation`: wall accommodation coefficient
"""
function reflect_particle_x!(rng, particle, wall_reflection_v_sq, wall_normal_sign, wall_v, wall_accommodation)
    if wall_accommodation == 0.0
        specular_reflection_x!(particle)
    elseif wall_accommodation == 1.0
        diffuse_reflection_x!(rng, particle, wall_reflection_v_sq, wall_normal_sign, wall_v)
    else
        R = rand(rng)
        if R < wall_accommodation
            diffuse_reflection_x!(rng, particle, wall_reflection_v_sq, wall_normal_sign, wall_v)
        else
            specular_reflection_x!(particle)
        end
    end
end