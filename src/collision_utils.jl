using SpecialFunctions: gamma
using TOML
using StaticArrays

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
    vhs_Tref::Float64
    vhs_d::Float64
    vhs_o::Float64
    vhs_factor::Float64

    # vhs_factor = π * vhs_d^2 * (2 * vhs_Tref/m_r)^(vhs_o - 0.5) / gamma(2.5 - vhs_o)
end

function load_interaction_data(interactions_filename, species_names, species_list)
    interactions_data = TOML.parsefile(species_filename)

    interactions_list = Array{Interaction,2}(undef,length(species_names),length(species_names))

    for (i, species_name1) in enumerate(species_names)
        for (k, species_name2) in enumerate(species_names)
            if i >= k
                m_r = species_list[i].mass * species_list[k].mass / (species_list[i].mass + species_list[k].mass)

                μ1 = species_list[i].mass / (species_list[i].mass + species_list[k].mass) 
                μ2 = species_list[i].mass / (species_list[i].mass + species_list[k].mass) 
                # interactions_list[i,k] = Interaction(m_r, μ1, μ2)
            end
            # push!(species_list, Species(species_name, species_data[species_name]["mass"]))
        end
    end
end

function create_collision_data()
    nothing
end

function compute_com(collision_data::CollisionData, interaction::Interaction, p1, p2)
    collision_data.v_com = interaction.μ1 * p1.v + interaction.μ2 * p2.v
end

function compute_g(collision_data::CollisionData, p1, p2)
    collision_data.g_vec = p1.v - p2.v
    collision_data.g = norm(collision_data.g_vec)
end