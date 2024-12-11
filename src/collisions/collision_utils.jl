using SpecialFunctions: gamma
using TOML
using StaticArrays

"""
Used to store temporary collision data for a specific collision
"""
mutable struct CollisionData
    v_com::SVector{3,Float64}
    g::Float64
    E_coll::Float64  # in case we also need to store the collision energy
    E_coll_eV::Float64  # in case we also need to store the collision energy in eV
    E_coll_electron_eV::Float64  # collision energy of impacting electron in eV
    g_vec::SVector{3,Float64}
    g_vec_new::SVector{3,Float64}
    g_new_1::Float64
    g_new_2::Float64  # to store post-collision energies of multiple particles (i.e. during chemical reactions)
end

"""
Used to store Fokker-Planck collision data
"""
mutable struct CollisionDataFP
    vel_ave::SVector{3,Float64}
    mean::SVector{3,Float64}
    stddev::SVector{3,Float64}
end

"""
Used to store interaction parameters for a 2-species interaction
"""
struct Interaction
    m_r::Float64
    μ1::Float64  # m1 / (m1 + m2)
    μ2::Float64  # m2 / (m1 + m2)
    vhs_d::Float64
    vhs_o::Float64
    vhs_Tref::Float64
    vhs_muref::Float64
    vhs_factor::Float64 # = π * vhs_d^2 * (2 * vhs_Tref/m_r)^(vhs_o - 0.5) / gamma(2.5 - vhs_o)
end

"""
Compute interaction-specific factor ``\\pi D_{VHS}^2 (2 T_{ref,VHS}/m_r)^{(\\omega_{VHS} - 0.5)} \\frac{1}{\\Gamma(2.5 - \\omega_{VHS})}``
"""
function compute_vhs_factor(vhs_Tref, vhs_d, vhs_o, m_r)
    # π * vhs_d^2 * ((2 * k_B * vhs_Tref/m_r)^(vhs_o - 0.5)) / gamma(2.5 - vhs_o), ", ", π * vhs_d^2)
    return π * vhs_d^2 * (2 * k_B * vhs_Tref/m_r)^(vhs_o - 0.5) / gamma(2.5 - vhs_o)
end

"""
Create empty CollisionData instance
"""
CollisionData() = CollisionData(SVector{3,Float64}(0.0, 0.0, 0.0), 0.0, 0.0, 0.0, 0.0,
                                SVector{3,Float64}(0.0, 0.0, 0.0),
                                SVector{3,Float64}(0.0, 0.0, 0.0), 0.0, 0.0)

"""
Create empty CollisionDataFP instance
"""
CollisionDataFP() = CollisionDataFP(SVector{3,Float64}(0.0, 0.0, 0.0),
                                    SVector{3,Float64}(0.0, 0.0, 0.0),
                                    SVector{3,Float64}(0.0, 0.0, 0.0))

"""
Load interaction data
"""
function load_interaction_data(interactions_filename, species_data)
    interactions_data = TOML.parsefile(interactions_filename)

    species_names = [species.name for species in species_data]

    interactions_list = Array{Interaction,2}(undef,length(species_names),length(species_names))

    for (i, species_name1) in enumerate(species_names)
        for (k, species_name2) in enumerate(species_names)
            if i >= k

                s1s2 = species_name1 * "," * species_name2

                m_r = species_data[i].mass * species_data[k].mass / (species_data[i].mass + species_data[k].mass)

                μ1 = species_data[i].mass / (species_data[i].mass + species_data[k].mass) 
                μ2 = species_data[k].mass / (species_data[i].mass + species_data[k].mass) 

                interaction_s1s2 = nothing
                try
                    interaction_s1s2 = interactions_data[s1s2]
                catch KeyError
                    s2s1 = species_name2 * "," * species_name1
                    interaction_s1s2 = interactions_data[s2s1]
                end

                interactions_list[i,k] = Interaction(m_r, μ1, μ2,
                interaction_s1s2["vhs_d"], interaction_s1s2["vhs_o"],
                interaction_s1s2["vhs_Tref"],
                compute_mu_ref(0.5*(species_data[i].mass + species_data[k].mass), interaction_s1s2["vhs_o"], interaction_s1s2["vhs_Tref"], interaction_s1s2["vhs_d"]),
                compute_vhs_factor(interaction_s1s2["vhs_Tref"], interaction_s1s2["vhs_d"], interaction_s1s2["vhs_o"],
                m_r))
                interactions_list[k,i] = Interaction(m_r, μ2, μ1,
                interaction_s1s2["vhs_d"], interaction_s1s2["vhs_o"],
                interaction_s1s2["vhs_Tref"],
                compute_mu_ref(0.5*(species_data[i].mass + species_data[k].mass), interaction_s1s2["vhs_o"], interaction_s1s2["vhs_Tref"], interaction_s1s2["vhs_d"]),
                compute_vhs_factor(interaction_s1s2["vhs_Tref"], interaction_s1s2["vhs_d"], interaction_s1s2["vhs_o"],
                m_r))
            end
        end
    end

    return interactions_list
end

"""
Compute reference viscosity
"""
function compute_mu_ref(mass, omega, Tref, diameter)
    numerator = 30.0 * sqrt(mass * k_B * Tref)
    denumerator = 4.0 * sqrt(π) * (5.0 - 2.0 * omega) * (7.0 - 2.0 * omega) * diameter * diameter

    return numerator / denumerator;
end

"""
Load interaction data filling with dummy data if needed
"""
function load_interaction_data(interactions_filename, species_data, dummy_vhs_d, dummy_vhs_o, dummy_vhs_Tref)
    interactions_data = TOML.parsefile(interactions_filename)

    species_names = [species.name for species in species_data]

    interactions_list = Array{Interaction,2}(undef,length(species_names),length(species_names))

    for (i, species_name1) in enumerate(species_names)
        for (k, species_name2) in enumerate(species_names)
            if i >= k

                s1s2 = species_name1 * "," * species_name2

                m_r = species_data[i].mass * species_data[k].mass / (species_data[i].mass + species_data[k].mass)

                μ1 = species_data[i].mass / (species_data[i].mass + species_data[k].mass) 
                μ2 = species_data[k].mass / (species_data[i].mass + species_data[k].mass) 

                interaction_s1s2 = nothing
                try
                    interaction_s1s2 = interactions_data[s1s2]
                catch KeyError
                    s2s1 = species_name2 * "," * species_name1
                    try
                        interaction_s1s2 = interactions_data[s2s1]
                    catch KeyError
                        nothing
                    end
                end
                
                if interaction_s1s2 !== nothing
                    interactions_list[i,k] = Interaction(m_r, μ1, μ2,
                    interaction_s1s2["vhs_d"], interaction_s1s2["vhs_o"],
                    interaction_s1s2["vhs_Tref"],
                    compute_mu_ref(0.5*(species_data[i].mass + species_data[k].mass), interaction_s1s2["vhs_o"], interaction_s1s2["vhs_Tref"], interaction_s1s2["vhs_d"]),
                    compute_vhs_factor(interaction_s1s2["vhs_Tref"], interaction_s1s2["vhs_d"], interaction_s1s2["vhs_o"],
                    m_r))
                    
                    interactions_list[k,i] = Interaction(m_r, μ2, μ1,
                    interaction_s1s2["vhs_d"], interaction_s1s2["vhs_o"],
                    interaction_s1s2["vhs_Tref"],
                    compute_mu_ref(0.5*(species_data[i].mass + species_data[k].mass), interaction_s1s2["vhs_o"], interaction_s1s2["vhs_Tref"], interaction_s1s2["vhs_d"]),
                    compute_vhs_factor(interaction_s1s2["vhs_Tref"], interaction_s1s2["vhs_d"], interaction_s1s2["vhs_o"],
                    m_r))
                else
                    interactions_list[i,k] = Interaction(m_r, μ1, μ2,
                    dummy_vhs_d, dummy_vhs_o,
                    dummy_vhs_Tref,
                    compute_mu_ref(0.5*(species_data[i].mass + species_data[k].mass), dummy_vhs_o, dummy_vhs_Tref, dummy_vhs_d),
                    compute_vhs_factor(dummy_vhs_Tref, dummy_vhs_d, dummy_vhs_o, m_r))
                    
                    interactions_list[k,i] = Interaction(m_r, μ2, μ1,
                    dummy_vhs_d, dummy_vhs_o,
                    dummy_vhs_Tref,
                    compute_mu_ref(0.5*(species_data[i].mass + species_data[k].mass), dummy_vhs_o, dummy_vhs_Tref, dummy_vhs_d),
                    compute_vhs_factor(dummy_vhs_Tref, dummy_vhs_d, dummy_vhs_o, m_r))
                end
            end
        end
    end

    return interactions_list
end

"""
Load interaction data filling with dummy data if needed
"""
function load_interaction_data_with_dummy(interactions_filename, species_data)
    return load_interaction_data(interactions_filename, species_data, 1e-10, 1.0, 273.0)
end

"""
Compute center of mass velocity
"""
function compute_com!(collision_data::CollisionData, interaction::Interaction, p1, p2)
    collision_data.v_com = interaction.μ1 * p1.v + interaction.μ2 * p2.v
end

"""
Compute relative velocity
"""
function compute_g!(collision_data::CollisionData, p1, p2)
    collision_data.g_vec = p1.v - p2.v
    collision_data.g = norm(collision_data.g_vec)
end

"""
Estimate ``(\\sigma g w)_{max}`` for a two-species interaction
"""
function estimate_sigma_g_w_max(interaction, species1, species2, T1, T2, Fnum; mult_factor=1.0)
    g_thermal1 = sqrt(2 * T1 * k_B / species1.mass)
    g_thermal2 = sqrt(2 * T2 * k_B / species2.mass)
    g_thermal = 0.5 * (g_thermal1 + g_thermal2)
    return mult_factor * sigma_vhs(interaction, g_thermal) * g_thermal * Fnum
end

"""
Estimate ``(\\sigma g w)_{max}`` for a single-species interaction
"""
function estimate_sigma_g_w_max(interaction, species, T, Fnum; mult_factor=1.0)
    return estimate_sigma_g_w_max(interaction, species, species, T, T, Fnum, mult_factor=mult_factor)
end

"""
Estimate ``(\\sigma g w)_{max}`` for all species
"""
function estimate_sigma_g_w_max!(collision_factors, interactions, species_data, T_list, Fnum; mult_factor=1.0)
    for (k, species2) in enumerate(species_data)
        for (i, species1) in enumerate(species_data)
            sigma_g_w_max = estimate_sigma_g_w_max(interactions[i,k], species1, species2, T_list[i], T_list[k], Fnum, mult_factor=mult_factor)

            for cell in 1:length(collision_factors[i,k,:])
                collision_factors[i,k,cell].sigma_g_w_max = sigma_g_w_max
            end
        end
    end
end

"""
Compute post-ionization relative velocity
"""
function compute_g_new_ionization!(coll_data, interaction, E_i, energy_splitting)
    E_new_coll = coll_data.E_coll_electron_eV - E_i

    if (energy_splitting == ElectronEnergySplitEqual)
        coll_data.g_new_1 = sqrt(E_new_coll / e_mass_div_electron_volt)
        coll_data.g_new_2 = coll_data.g_new_1
    else
        coll_data.g_new_1 = sqrt(2 * E_new_coll / e_mass_div_electron_volt)
        coll_data.g_new_2 = 0.0
    end
end