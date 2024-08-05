struct MaxwellWallBC
    T::Float64
    v::SVector{3,Float64}
    accomodation::Float64  # accomodation coefficient, = 0: specular, = 1: purely diffuse
end

struct MaxwellWalls
    boundaries::Vector{MaxwellWallBC}
    # reflection_velocities::Array{Float64,2} per-wall, per-species
end

function create_1D_boundaries(T_l, T_r, vy_l, vy_r, accomodation_l, accomodation_r)
    return MaxwellWalls([MaxwellWallBC(T_l, [0.0, vy_l, 0.0], accomodation_l),
                         MaxwellWallBC(T_r, [0.0, vy_r, 0.0], accomodation_r)])
end

function specular_reflection_x!(particle)
    particle.v = SVector{3, Float64}(-particle.v[1], particle.v[2], particle.v[3])
end

function diffuse_reflection_x!(rng, particle, species_data, wall_normal_sign, wall_T, wall_v)
    R = max(1e-50, rand(rng))

    wall_char_v = 2 * k_B * wall_T / species_data.mass # TODO: Precompute!
    v_normal = wall_normal_sign * sqrt(-wall_char_v * log(R))
    
    R = max(1e-50, rand(rng))
    v_tang = sqrt(-wall_char_v * log(R))

    R = rand(rng)
    v_tang1 = sin(twopi * R) * v_tang 
    v_tang2 = cos(twopi * R) * v_tang

    particle.v = SVector{3, Float64}(v_normal + wall_v[1], v_tang1 + wall_v[2], v_tang2 + wall_v[3])
end

function reflect_particle_x!(rng, particle, species_data, wall_normal_sign, wall_T, wall_v, wall_accomodation)
    if wall_accomodation == 0.0
        specular_reflection_x!(particle)
    elseif wall_accomodation == 1.0
        diffuse_reflection_x!(rng, particle, species_data, wall_normal_sign, wall_T, wall_v)
    else
        R = rand(rng)
        if R > wall_accomodation
            diffuse_reflection_x!(rng, particle, species_data, wall_normal_sign, wall_T, wall_v)
        else
            specular_reflection_x!(particle)
        end
    end
end