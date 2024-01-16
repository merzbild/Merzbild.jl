using SpecialFunctions: gamma
using TOML
using StaticArrays

mutable struct CollisionData
    v_com::SVector{3,Float64}
    g::Float64
    g_vec::SVector{3,Float64}
    g_vec_new::SVector{3,Float64}
end

struct Interaction
    m_r::Float64
    μ1::Float64  # m1 / (m1 + m2)
    μ2::Float64  # m2 / (m1 + m2)
    # vhs_Tref::Float64
    vhs_d::Float64
    vhs_o::Float64
    vhs_factor::Float64
    # vhs_factor = π * vhs_d^2 * (2 * vhs_Tref/m_r)^(vhs_o - 0.5) / gamma(2.5 - vhs_o)
end

function compute_vhs_factor(vhs_Tref, vhs_d, vhs_o, m_r)
    # println(m_r, ", ", vhs_o, ", ", vhs_Tref, ", ", (2 * k_B * vhs_Tref/m_r)^(vhs_o - 0.5), ", ",
    # π * vhs_d^2 * ((2 * k_B * vhs_Tref/m_r)^(vhs_o - 0.5)) / gamma(2.5 - vhs_o), ", ", π * vhs_d^2)
    return π * vhs_d^2 * (2 * k_B * vhs_Tref/m_r)^(vhs_o - 0.5) / gamma(2.5 - vhs_o)
end

function load_interaction_data(interactions_filename, species_list)
    interactions_data = TOML.parsefile(interactions_filename)

    species_names = [species.name for species in species_list]

    interactions_list = Array{Interaction,2}(undef,length(species_names),length(species_names))

    for (i, species_name1) in enumerate(species_names)
        for (k, species_name2) in enumerate(species_names)
            if i >= k

                s1s2 = species_name1 * "," * species_name2

                m_r = species_list[i].mass * species_list[k].mass / (species_list[i].mass + species_list[k].mass)

                μ1 = species_list[i].mass / (species_list[i].mass + species_list[k].mass) 
                μ2 = species_list[k].mass / (species_list[i].mass + species_list[k].mass) 

                interaction_s1s2 = nothing
                try
                    interaction_s1s2 = interactions_data[s1s2]
                catch KeyError
                    s2s1 = species_name2 * "," * species_name1
                    interaction_s1s2 = interactions_data[s2s1]
                end
                println(i, ", ", k)
                interactions_list[i,k] = Interaction(m_r, μ1, μ2, # interaction_s1s2["vhs_Tref"],
                interaction_s1s2["vhs_d"], interaction_s1s2["vhs_o"],
                compute_vhs_factor(interaction_s1s2["vhs_Tref"], interaction_s1s2["vhs_d"], interaction_s1s2["vhs_o"],
                m_r))
                interactions_list[k,i] = interactions_list[i,k]
            end
            # push!(species_list, Species(species_name, species_data[species_name]["mass"]))
        end
    end

    return interactions_list
end

function create_collision_data()
    return CollisionData(SVector{3,Float64}(0.0, 0.0, 0.0), 0.0,
                         SVector{3,Float64}(0.0, 0.0, 0.0),
                         SVector{3,Float64}(0.0, 0.0, 0.0))
end

function compute_com!(collision_data::CollisionData, interaction::Interaction, p1, p2)
    collision_data.v_com = interaction.μ1 * p1.v + interaction.μ2 * p2.v
end

function compute_g!(collision_data::CollisionData, p1, p2)
    collision_data.g_vec = p1.v - p2.v
    collision_data.g = norm(collision_data.g_vec)
end

function estimate_sigma_g_w_max(interaction, species, T, Fnum; mult_factor=1.0)
    g_thermal = sqrt(2 * T * k_B / species.mass)
    return mult_factor * sigma_vhs(interaction, g_thermal) * g_thermal * Fnum
end

function estimate_sigma_g_w_max(interaction, species1, species2, T1, T2, Fnum; mult_factor=1.0)
    g_thermal1 = sqrt(2 * T1 * k_B / species1.mass)
    g_thermal2 = sqrt(2 * T2 * k_B / species2.mass)
    g_thermal = 0.5 * (g_thermal1 + g_thermal2)
    return mult_factor * sigma_vhs(interaction, g_thermal) * g_thermal * Fnum
end

function estimate_sigma_g_w_max!(collision_factors, interactions, species_list, T_list, Fnum; mult_factor=1.0)
    for (k, species2) in enumerate(species_list)
        for (i, species1) in enumerate(species_list)
            g_thermal1 = sqrt(2 * T_list[1] * k_B / species1.mass)
            g_thermal2 = sqrt(2 * T_list[2] * k_B / species2.mass)
            g_thermal = 0.5 * (g_thermal1 + g_thermal2)
            collision_factors[i,k].sigma_g_w_max = mult_factor * sigma_vhs(interactions[i,k], g_thermal) * g_thermal * Fnum
        end
    end
end