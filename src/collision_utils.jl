mutable struct CollisionData
    v_com::MVector{3,Float64}
    g::Float64
    g_vec::MVector{3,Float64}
    g_vec_new::MVector{3,Float64}
end

struct Interaction
    m_r::Float64
    μ1::Float64  # m1 / (m1 + m2)
    μ2::Float64  # m2 / (m1 + m2)
    vhs_d::Float64
    vhs_o::Float64
    vhs_factor::Float64
end

function compute_com(collision_data::CollisionData, interaction::Interaction, p1, p2)
    collision_data.v_com .= interaction.μ1 * p1.v .+ interaction.μ2 * p2.v
end

function compute_g(collision_data::CollisionData, p1, p2)
    collision_data.g_vec = p1.v .- p2.v
    collision_data.g = norm(collision_data.g_vec)
end