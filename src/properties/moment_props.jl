using LinearAlgebra

@muladd begin

"""
    compute_moment_scaling!(moment_scaling, moment_powers, species, species_data, Tref)

Fill the vector `moment_scaling` with the appropriate scaling based on the Maxwell-Boltzmann distribution
for a specific species.

# Positional arguments
* `moment_scaling`: vector to be filled with the moment scaling factors (length `n_moments`)
* `moment_powers`: vector of moment powers (length `n_moments`)
* `species`: the species for which the mixed moment is being computed
* `species_data`: the `Vector` of `SpeciesData`
* `Tref`: reference temperature used to compute the total moments of a Maxwellian distribution which
    are used to scale the total moments
"""
function compute_moment_scaling!(moment_scaling, moment_powers, species, species_data, Tref)
    moment_factor = 4 * π * (species_data[species].mass / (twopi * k_B * Tref))^(1.5) * 0.5
    moment_vref = (species_data[species].mass / (2 * k_B * Tref))^0.5
    
    l_mom = length(moment_powers)

    @inbounds for n_mom in 1:l_mom
        m = moment_powers[n_mom]
        moment_scaling[n_mom] = moment_factor * moment_vref^(-(3 + m)) * gamma((3 + m) / 2)
    end
end

"""
    compute_moments!(moment_values, moment_scaling, moment_powers, particles, pia, cell, species, species_data, phys_props)

Fill the vector `moment_values` with the computed scaled moment values for a single species in a single cell.
The values of `n` (number density or number of physical particles) and `v` (mean velocity) are taken from `phys_props`.
It's important that function `compute_props!` must be called before this function, so that the correct `n`` and `v`` values are available.

# Positional arguments
* `moment_values`: vector to be filled with the computed moment values (length `n_moments`)
* `moment_scaling`: vector of precomputed moment scaling factors (length `n_moments`)
* `moment_powers`: vector of moment powers (length `n_moments`)
* `particles`: the `Vector` of `ParticleVector`s containing all the particles in a simulation
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the cell in which the moment is being computed
* `species`: the species for which the mixed moment is being computed
* `species_data`: the `Vector` of `SpeciesData`
* `phys_props`: the `PhysProps` instance in which the computed physical properties are stored
"""
function compute_moments!(moment_values, moment_scaling, moment_powers, particles, pia, cell, species, species_data, phys_props)
    # Initialize moment values to zero
    fill!(moment_values, 0.0)
    
    # Get the mean velocity from phys_props
    @inbounds v_mean = SVector{3,Float64}(phys_props.v[1, cell, species],
                                          phys_props.v[2, cell, species],
                                          phys_props.v[3, cell, species])
    
    # Get number density / number of physical particles from phys_props
    @inbounds n = phys_props.n[cell, species]
    @inbounds indexer = pia.indexer[cell, species]
    
    if (n > 0.0)
        # Accumulate moments from group1 particles
        s1 = indexer.start1
        e1 = indexer.end1

        @inbounds for i in s1:e1
            normv = norm(particles[species][i].v - v_mean)
            for (n_mom, m) in enumerate(moment_powers)
                moment_values[n_mom] += particles[species][i].w * normv^m
            end
        end
        
        # Accumulate moments from group2 particles if present
        if indexer.n_group2 > 0
            s2 = indexer.start2
            e2 = indexer.end2
            @inbounds for i in s2:e2
                normv = norm(particles[species][i].v - v_mean)
                for (n_mom, m) in enumerate(moment_powers)
                    moment_values[n_mom] += particles[species][i].w * normv^m
                end
            end
        end
        
        # Normalize by scaling and n
        @inbounds for n_mom in eachindex(moment_values)
            moment_values[n_mom] /= (moment_scaling[n_mom] * n)
        end
    end
end

end

