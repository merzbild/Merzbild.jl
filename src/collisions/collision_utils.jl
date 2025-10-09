using SpecialFunctions: gamma
using TOML
using StaticArrays

@muladd begin

"""
    CollisionData

Structure to store temporary collision data for a specific collision (relative velocity, collision energy, post-collision velocities, etc.)

# Fields
* `v_com`: vector of center-of-mass velocity
* `g`: magnitude of relative velocity
* `E_coll`: relative translational energy of the colliding particles
* `E_coll_eV`: relative translational energy of the colliding particles in electron-volt
* `E_coll_electron_eV`: collisional energy of the an electron in electron-volt for electron-neutral collisions
* `g_vec`: vector of pre-collisional relative velocity
* `g_vec_new`: vector of post-collisional relative velocity
* `g_new_1`: magnitude of post-collisional relative velocity of the first particle in the collision pair
* `g_new_2`: magnitude of post-collisional relative velocity of the second particle in the collision pair
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
    CollisionData

Structure to store temporary collision data for the particle Fokker-Planck approach.

# Fields
* `vel_ave`: average velocity of the particles in the cell
* `mean`: mean value of the sampled velocities
* `stddev`: standard deviation of the sampled velocities
* `xvel_rand`: pre-allocated storage for sampled x-velocity components
* `yvel_rand`: pre-allocated storage for sampled y-velocity components
* `zvel_rand`: pre-allocated storage for sampled z-velocity components
"""
mutable struct CollisionDataFP
    vel_ave::SVector{3,Float64}
    mean::SVector{3,Float64}
    stddev::SVector{3,Float64}
    xvel_rand::Vector{Float64}
    yvel_rand::Vector{Float64}
    zvel_rand::Vector{Float64}
end

"""
    Interaction

Structure to store interaction parameters for a 2-species interaction.
The VHS model uses the following power law: ``\\sigma_{VHS} = C g^(1 - 2 \\omega_{VHS})``, where
``omega`` is the exponent of the VHS potential, and ``C`` is the pre-computed factor:
``C = \\pi D_{VHS}^2 (2 T_{ref,VHS}/m_r)^{(\\omega_{VHS} - 0.5)} \\frac{1}{\\Gamma(2.5 - \\omega_{VHS})}``.

# Fields
* `m_r`: collision-reduced mass
* `μ1`: relative mass of the first species
* `μ2`: relative mass of the second species
* `vhs_d`: diameter for the VHS potential
* `vhs_o`: exponent for the VHS potential
* `vhs_Tref`: reference temperature for the VHS potential
* `vhs_muref`: reference viscosity for the VHS potential
* `vhs_factor`: pre-computed factor for calculation of the VHS cross-section
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
    compute_vhs_factor(vhs_Tref, vhs_d, vhs_o, m_r)

Compute interaction-specific factor ``\\pi D_{VHS}^2 (2 T_{ref,VHS}/m_r)^{(\\omega_{VHS} - 0.5)} \\frac{1}{\\Gamma(2.5 - \\omega_{VHS})}``.

# Positional arguments
* `vhs_Tref`: reference temperature for the VHS potential
* `vhs_d`: diameter for the VHS potential
* `vhs_o`: exponent for the VHS potential
* `m_r`: collision-reduced mass

# Returns
* computed VHS cross-section factor
"""
function compute_vhs_factor(vhs_Tref, vhs_d, vhs_o, m_r)
    # π * vhs_d^2 * ((2 * k_B * vhs_Tref/m_r)^(vhs_o - 0.5)) / gamma(2.5 - vhs_o), ", ", π * vhs_d^2)
    return π * vhs_d^2 * (2 * k_B * vhs_Tref/m_r)^(vhs_o - 0.5) / gamma(2.5 - vhs_o)
end

"""
    CollisionData() 

Create an empty CollisionData instance.
"""
CollisionData() = CollisionData(SVector{3,Float64}(0.0, 0.0, 0.0), 0.0, 0.0, 0.0, 0.0,
                                SVector{3,Float64}(0.0, 0.0, 0.0),
                                SVector{3,Float64}(0.0, 0.0, 0.0), 0.0, 0.0)

"""
    CollisionDataFP()
    
Create an empty CollisionDataFP instance.
"""
CollisionDataFP() = CollisionDataFP(SVector{3,Float64}(0.0, 0.0, 0.0),
                                    SVector{3,Float64}(0.0, 0.0, 0.0),
                                    SVector{3,Float64}(0.0, 0.0, 0.0),
                                    [0.0], [0.0], [0.0])

"""
    CollisionDataFP(n_particles_in_cell)
    
Create an empty CollisionDataFP instance, pre-allocating the arrays for sampled normal variables
for `n_particles_in_cell`.

# Positional arguments
* `n_particles_in_cell`: estimate of expected maximum number of particles in cell
"""
CollisionDataFP(n_particles_in_cell) = CollisionDataFP(SVector{3,Float64}(0.0, 0.0, 0.0),
                                    SVector{3,Float64}(0.0, 0.0, 0.0),
                                    SVector{3,Float64}(0.0, 0.0, 0.0),
                                    zeros(n_particles_in_cell),
                                    zeros(n_particles_in_cell),
                                    zeros(n_particles_in_cell))

"""
    load_interaction_data(interactions_filename, species_data)

Load interaction data from a TOML file given a list of species' data (list of `Species` instances).
It will load interaction data for all possible pair-wise interactions of the species in the list.

The resulting 2-D array has the interaction data for `Species[i]` with `Species[k]` in position `[i,k]`.
It is not symmetric, as the relative collision masses `μ1` and `μ2` are swapped when
comparing the `Interaction` instances in positions `[i,k]` and `[k,i]`. If no data is found,
the function throws an error.

# Positional arguments
* `interactions_filename`: the path to the TOML file containing the data
* `species_data`: list of `Species` instances for which to search for the interaction data

# Returns
* 2-dimensional array of `Interaction` instances of size `(n_species, n_species)`

# Throws
`KeyError` if interaction data not found in the file.
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
    compute_mu_ref(m_r, vhs_o, vhs_Tref, vhs_d)
    
Compute reference viscosity for the VHS model.

# Positional arguments
* `m_r`: collision-reduced mass
* `vhs_o`: exponent for the VHS potential
* `vhs_Tref`: reference temperature for the VHS potential
* `vhs_d`: diameter for the VHS potential

# Returns
* reference viscosity
"""
function compute_mu_ref(m_r, vhs_o, vhs_Tref, vhs_d)
    numerator = 30.0 * sqrt(m_r * k_B * vhs_Tref)
    denumerator = 4.0 * sqrt(π) * (5.0 - 2.0 * vhs_o) * (7.0 - 2.0 * vhs_o) * vhs_d * vhs_d

    return numerator / denumerator;
end

"""
    load_interaction_data(interactions_filename, species_data)

Load interaction data from a TOML file given a list of species' data (list of `Species` instances),
filling in dummy VHS data in case no entry is found in the TOML file. Useful for interactions where
the VHS model doesn't make sense, for example electron-neutral interactions.

It will load interaction data for all possible pair-wise interactions of the species in the list.
The resulting 2-D array has the interaction data for `Species[i]` with `Species[k]` in position `[i,k]`.
It is not symmetric, as the relative collision masses `μ1` and `μ2` are swapped when
comparing the `Interaction` instances in positions `[i,k]` and `[k,i]`.

# Positional arguments
* `interactions_filename`: the path to the TOML file containing the data
* `species_data`: list of `Species` instances for which to search for the interaction data
* `dummy_vhs_d`: value to use for the VHS diameter if no interaction data found in the file
* `dummy_vhs_o`: value to use for the VHS exponent if no interaction data found in the file
* `dummy_vhs_Tref`: value to use for the VHS reference temperature if no interaction data found in the file

# Returns
* 2-dimensional array of `Interaction` instances of size `(n_species, n_species)`
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
    load_interaction_data_with_dummy(interactions_filename, species_data)

Load interaction data from a TOML file given a list of species' data (list of `Species` instances),
filling in dummy VHS data in case no entry is found in the TOML file. Useful for interactions where
the VHS model doesn't make sense, for example electron-neutral interactions.
Uses a value of `1e-10` for the dummy VHS diameter, `1.0` for the dummy VHS exponent, and `273.0` for the
dummy VHS reference temperature.

It will load interaction data for all possible pair-wise interactions of the species in the list.
The resulting 2-D array has the interaction data for `Species[i]` with `Species[k]` in position `[i,k]`.
It is not symmetric, as the relative collision masses `μ1` and `μ2` are swapped when
comparing the `Interaction` instances in positions `[i,k]` and `[k,i]`.
    
# Positional arguments
* `interactions_filename`: the path to the TOML file containing the data
* `species_data`: list of `Species` instances for which to search for the interaction data

# Returns
* 2-dimensional array of `Interaction` instances of size `(n_species, n_species)`
"""
function load_interaction_data_with_dummy(interactions_filename, species_data)
    return load_interaction_data(interactions_filename, species_data, 1e-10, 1.0, 273.0)
end


"""
    load_species_and_interaction_data(species_filename, interactions_filename, species_names; fill_dummy=true)

Given a list of species' names, load the species and interaction data (filling with dummy data if needed).

# Positional arguments
* `species_filename`: the path to the TOML file containing the species' data
* `interactions_filename`: the path to the TOML file containing the interaction data
* `species_name`: the name of the species for which to load the data

# Keyword arguments
* `fill_dummy`: if `true`, fill interaction data with computed dummy values if no entry found for a species pair in the interaction file

# Returns
* `Vector` of `Species` instances
* 2-dimensional array of `Interaction` instances of size `(n_species, n_species)`
"""
function load_species_and_interaction_data(species_filename, interactions_filename, species_names; fill_dummy=true)
    species_data = load_species_data(species_filename, species_names)
    if fill_dummy
        return species_data, load_interaction_data_with_dummy(interactions_filename, species_data)
    else
        return species_data, load_interaction_data(interactions_filename, species_data)
    end
end

"""
    compute_com!(collision_data::CollisionData, interaction::Interaction, p1, p2)

Compute center of mass velocity of two particles.

# Positional arguments
* `collision_data`: the `CollisionData` instance where the computed velocity will be stored
* `interaction`: the `Interaction` instance for the species of the colliding particles
* `p1`: the first particle
* `p2`: the second particle
"""
@inline function compute_com!(collision_data::CollisionData, interaction::Interaction, p1, p2)
    collision_data.v_com = interaction.μ1 * p1.v + interaction.μ2 * p2.v
end

"""
    compute_g!(collision_data::CollisionData, p1, p2)

Compute relative velocity (vector and its magnitude) of two particles.

# Positional arguments
* `collision_data`: the `CollisionData` instance where the computed velocity will be stored
* `p1`: the first particle
* `p2`: the second particle
"""
@inline function compute_g!(collision_data::CollisionData, p1, p2)
    collision_data.g_vec = p1.v - p2.v
    collision_data.g = norm(collision_data.g_vec)
end

"""
    estimate_sigma_g_w_max(interaction, species1, species2, T1, T2, Fnum; mult_factor=1.0)

Estimate ``(\\sigma g w)_{max}`` for a two-species interaction, assuming a constant particle computational weight `Fnum` and
a VHS cross-section.
The relative velocity ``g`` is estimated as ``g = 0.5 (\\sqrt{2T_1 k_B / m_1} + \\sqrt{2T_2 k_B / m_2})``,
where ``T_1`` and ``m_1`` are the temperature and mass of the first species (`species1`), and 
``T_2`` and ``m_2`` are the temperature and mass of the second species (`species2`).
This relative velocity estimate is then plugged into the VHS cross-section model to compute ``\\sigma``.
The result is then multiplied by `Fnum` and an (optional) factor `mult_factor`.

# Positional arguments
* `interaction`: the `Interaction` instance for the interacting species
* `species1`: the `Species` instance of the first interacting species
* `species2`: the `Species` instance of the second interacting species
* `T1`: the temperature of the first interacting species
* `T2`: the temperature of the second interacting species
* `Fnum`: the constant computational weight of the particles

# Keyword arguments
* `mult_factor`: a factor by which to multiply the result (default value is 1.0)

# Returns
The estimate of ``(\\sigma g w)_{max}``.
"""
function estimate_sigma_g_w_max(interaction, species1, species2, T1, T2, Fnum; mult_factor=1.0)
    g_thermal1 = sqrt(2 * T1 * k_B / species1.mass)
    g_thermal2 = sqrt(2 * T2 * k_B / species2.mass)
    g_thermal = 0.5 * (g_thermal1 + g_thermal2)
    return mult_factor * sigma_vhs(interaction, g_thermal) * g_thermal * Fnum
end

"""
    estimate_sigma_g_w_max(interaction, species, T, Fnum; mult_factor=1.0)

Estimate ``(\\sigma g w)_{max}`` for a single-species interaction, assuming a constant particle computational weight `Fnum` and
a VHS cross-section. Uses the same methodology as the estimate for a two-species interaction.

# Positional arguments
* `interaction`: the `Interaction` instance for the interacting species
* `species`: the `Species` instance of the interacting species
* `T`: the temperature of the interacting species
* `Fnum`: the constant computational weight of the particles

# Keyword arguments
* `mult_factor`: a factor by which to multiply the result (default value is 1.0)

# Returns
* the estimate of ``(\\sigma g w)_{max}``
"""
function estimate_sigma_g_w_max(interaction, species, T, Fnum; mult_factor=1.0)
    return estimate_sigma_g_w_max(interaction, species, species, T, T, Fnum, mult_factor=mult_factor)
end

"""
    estimate_sigma_g_w_max!(collision_factors, interactions, species_data, T_list, Fnum; mult_factor=1.0)

Estimate ``(\\sigma g w)_{max}`` for all species in all cells, assuming
a constant particle computational weight `Fnum`, a VHS cross-section, and that each species' temperature is constant across all cells.
Uses the same methodology as the estimate for a two-species interaction.

# Positional arguments
* `collision_factors`: 3-dimensional array of `CollisionFactors` of shape `(n_species, n_species, n_cells)`
* `interactions`: the 2-dimensional array of `Interaction` instances (of shape `(n_species, n_species)`) of all the pair-wise interactions
* `species_data`: the vector of `Species` instances of the species in the flow 
* `T_list`: the list of temperatures of the species
* `Fnum`: the constant computational weight of the particles

# Keyword arguments
* `mult_factor`: a factor by which to multiply the result (default value is 1.0)
"""
function estimate_sigma_g_w_max!(collision_factors::Array{CollisionFactors}, interactions, species_data, T_list, Fnum::Number; mult_factor=1.0)
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
    estimate_sigma_g_w_max!(collision_factors, interactions, species_data, T_list, Fnum; mult_factor=1.0)

Estimate ``(\\sigma g w)_{max}`` for all species in all cells, assuming
species-specific computational weights `Fnum`, a VHS cross-section, and that each species' temperature is constant across all cells.
Uses the same methodology as the estimate for a two-species interaction.

# Positional arguments
* `collision_factors`: 3-dimensional array of `CollisionFactors` of shape `(n_species, n_species, n_cells)`
* `interactions`: the 2-dimensional array of `Interaction` instances (of shape `(n_species, n_species)`) of all the pair-wise interactions
* `species_data`: the vector of `Species` instances of the species in the flow 
* `T_list`: the list of temperatures of the species
* `Fnum`: a list of the computational weights of the species

# Keyword arguments
* `mult_factor`: a factor by which to multiply the result (default value is 1.0)
"""
function estimate_sigma_g_w_max!(collision_factors::Array{CollisionFactors}, interactions, species_data, T_list, Fnum::Vector; mult_factor=1.0)
    for (k, species2) in enumerate(species_data)
        for (i, species1) in enumerate(species_data)
            sigma_g_w_max = estimate_sigma_g_w_max(interactions[i,k], species1, species2, T_list[i], T_list[k],
                                                   max(Fnum[i], Fnum[k]), mult_factor=mult_factor)

            for cell in 1:length(collision_factors[i,k,:])
                collision_factors[i,k,cell].sigma_g_w_max = sigma_g_w_max
            end
        end
    end
end

"""
    estimate_sigma_g_max!(collision_factors::Array{CollisionFactorsSWPM}, interactions, species_data, T_list, Fnum; mult_factor=1.0)

Estimate ``(\\sigma g)_{max}`` for all species in all cells, assuming
a constant particle computational weight `Fnum`, a VHS cross-section, and that each species' temperature is constant across all cells.
Uses the same methodology as the estimate for a two-species interaction.

# Positional arguments
* `collision_factors`: 3-dimensional array of `CollisionFactorsSWPM` of shape `(n_species, n_species, n_cells)`
* `interactions`: the 2-dimensional array of `Interaction` instances (of shape `(n_species, n_species)`) of all the pair-wise interactions
* `species_data`: the vector of `Species` instances of the species in the flow 
* `T_list`: the list of temperatures of the species
* `Fnum`: the constant computational weight of the particles

# Keyword arguments
* `mult_factor`: a factor by which to multiply the result (default value is 1.0)
"""
function estimate_sigma_g_max!(collision_factors::Array{CollisionFactorsSWPM}, interactions, species_data, T_list; mult_factor=1.0)
    for (k, species2) in enumerate(species_data)
        for (i, species1) in enumerate(species_data)
            sigma_g_max = estimate_sigma_g_w_max(interactions[i,k], species1, species2, T_list[i], T_list[k], 1.0, mult_factor=mult_factor)

            for cell in 1:length(collision_factors[i,k,:])
                collision_factors[i,k,cell].sigma_g_max = sigma_g_max
            end
        end
    end
end

"""
    compute_g_new_ionization!(collision_data, interaction, E_i, energy_splitting)

Compute post-ionization magnitudes of the relative velocities of the impacting and newly created electrons.

# Positional arguments
* `collision_data`: the `CollisionData` instance where the computed velocity will be stored
* `interaction`: the `Interaction` instance for the electron-neutral interaction which produced the ion
* `E_i`: ionization energy of the neutral species in electron-volt
* `energy_splitting`: how the energy is divided among the electrons: if `ElectronEnergySplitEqual`, the energy is split equally,
    if `ElectronEnergySplitZeroE`, one electron takes all of the energy
"""
function compute_g_new_ionization!(collision_data, interaction, E_i, energy_splitting)
    E_new_coll = collision_data.E_coll_electron_eV - E_i

    if (energy_splitting == ElectronEnergySplitEqual)
        collision_data.g_new_1 = sqrt(E_new_coll / e_mass_div_electron_volt)
        collision_data.g_new_2 = collision_data.g_new_1
    else
        collision_data.g_new_1 = sqrt(2 * E_new_coll / e_mass_div_electron_volt)
        collision_data.g_new_2 = 0.0
    end
end

end