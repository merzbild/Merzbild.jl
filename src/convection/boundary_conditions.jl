"""
A struct to hold information about a diffuse reflecting wall
"""
struct MaxwellWallBC
    T::Float64
    v::SVector{3,Float64}
    accomodation::Float64  # accomodation coefficient, = 0: specular, = 1: purely diffuse
end

"""
A struct to hold information about all diffuse reflecting walls in the simulation
"""
struct MaxwellWalls1D
    boundaries::Vector{MaxwellWallBC}
    reflection_velocities_sq::Array{Float64,2} # per-wall, per-species
end

"""
Create a MaxwellWalls struct for a 1-D simulation with 2 walls ("left" and "right)
"""
MaxwellWalls1D(species_data, T_l::Float64, T_r::Float64, vy_l::Float64, vy_r::Float64, accomodation_l::Float64, accomodation_r::Float64) =
MaxwellWalls1D([MaxwellWallBC(T_l, [0.0, vy_l, 0.0], accomodation_l),
               MaxwellWallBC(T_r, [0.0, vy_r, 0.0], accomodation_r)],
               vcat([2 * k_B * T_l / sd.mass for sd in species_data]',
                    [2 * k_B * T_r / sd.mass for sd in species_data]'))

"""
Perform specular reflection of a particle in the x direction
"""
function specular_reflection_x!(particle)
    particle.v = SVector{3, Float64}(-particle.v[1], particle.v[2], particle.v[3])
end

"""
Perform diffuse reflection of a particle in the x direction
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
Reflect particle in the x direction 
"""
function reflect_particle_x!(rng, particle, wall_reflection_v_sq, wall_normal_sign, wall_v, wall_accomodation)
    if wall_accomodation == 0.0
        specular_reflection_x!(particle)
    elseif wall_accomodation == 1.0
        diffuse_reflection_x!(rng, particle, wall_reflection_v_sq, wall_normal_sign, wall_v)
    else
        R = rand(rng)
        if R < wall_accomodation
            diffuse_reflection_x!(rng, particle, wall_reflection_v_sq, wall_normal_sign, wall_v)
        else
            specular_reflection_x!(particle)
        end
    end
end