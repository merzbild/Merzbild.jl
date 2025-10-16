"""
    mean_free_path(interaction, species, n, T)

Compute elastic VHS mean free path for single-species collisions.

# Positional arguments
* `interactions`: the 2-dimensional array of `Interaction` instances (of shape `(n_species, n_species)`) of all the pair-wise interactions
* `species`: the species for which to compute the mean free path
* `n`: the number density
* `T`: the temperature

# Returns
Mean free path.

# References
* Eqn. (4.65) in "Molecular Gas Dynamics and the Direct Simulation of Gas Flows" or (D.21) in "Nonequilibrium Gas Dynamics and Molecular Simulation".
"""
function mean_free_path(interaction, species, n, T)
    @inbounds λ = sqrt(2.0) * π * (interaction[species, species].vhs_d^2) * n
    @inbounds λ *= (interaction[species, species].vhs_Tref / T)^(interaction[species, species].vhs_o - 0.5)
    return 1.0 / λ
end

"""
    mean_collision_frequency(interaction, species_data, species, n, T)

Compute elastic VHS mean collision frequency for single-species collisions.

# Positional arguments
* `interactions`: the 2-dimensional array of `Interaction` instances (of shape `(n_species, n_species)`) of all the pair-wise interactions
* `species_data`: the vector of `Species` instances of the species in the flow 
* `species`: the species for which to compute the mean free path
* `n`: the number density
* `T`: the temperature

# Returns
Mean collision frequency.

# References
* Eqn. (4.64) in "Molecular Gas Dynamics and the Direct Simulation of Gas Flows"
"""
function mean_collision_frequency(interaction, species_data, species, n, T)
    @inbounds nu = 4.0 * (interaction[species, species].vhs_d^2) * n
    @inbounds nu *= sqrt(π * k_B * interaction[species, species].vhs_Tref / species_data[species].mass)
    @inbounds nu *= (T / interaction[species, species].vhs_Tref)^(1.0 - interaction[species, species].vhs_o)
    return nu
end