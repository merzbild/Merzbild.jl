@muladd begin

"""
    CollisionFactorsSWPM

Structure to store SWPM-related collision factors for collisions between particles of two species
in a given cell.

# Fields
* `n1`: the number of particles of the first species in the cell
* `n2`: the number of particles of the second species in the cell
* `sigma_g_max`: estimate of the ``(\\sigma g)_{max}`` (``\\sigma`` is the total collision cross-section,
    ``g`` is the relative collision velocity
* `n_coll`: number of collisions to be tested
* `n_coll_performed`: number of collisions actually performed
"""
mutable struct CollisionFactorsSWPM
    n1::Int64
    n2::Int64
    sigma_g_max::Float64  # (σ(g) * g)_max
    n_coll::Int64  # number of collisions tested
    n_coll_performed::Int64  # number of collisions performed
end

"""
    CollisionFactorsSWPM()

Create an empty `CollisionFactorsSWPM` instance (all values set to 0).
"""
CollisionFactorsSWPM() = CollisionFactorsSWPM(0, 0, 0.0, 0, 0)

"""
    create_collision_factors_swpm_array(n_species)

Create a 3-dimensional array of SWPM collision factors for all interaction pairs for a 0-D case (1 spatial cell),
with shape `(n_species,n_species,1)`.

# Positional arguments
* `n_species`: number of species in the flow

# Returns
3-dimensional array of `CollisionFactorsSWPM` instances with shape `(n_species,n_species,1)`.
"""
function create_collision_factors_swpm_array(n_species)
    coll_factor_array = Array{CollisionFactorsSWPM, 3}(undef, (n_species, n_species, 1))
    for k in 1:n_species
        for i in 1:n_species
            coll_factor_array[i,k,1] = CollisionFactorsSWPM()
        end
    end
    return coll_factor_array
end

"""
    create_collision_factors_swpm_array(n_species, n_cells)

Create a 3-dimensional array of SWPM collision factors for all interaction pairs for all cells
in the simulation, with shape `(n_species,n_species,n_cells)`.

# Positional arguments
* `n_species`: number of species in the flow
* `n_cells`: number of cells in the simulation

# Returns
3-dimensional array of `CollisionFactorsSWPM` instances with shape `(n_species,n_species,n_cells)`.
"""
function create_collision_factors_swpm_array(n_species, n_cells)
    coll_factor_array = Array{CollisionFactorsSWPM, 3}(undef, (n_species, n_species, n_cells))
    for k in 1:n_cells
        for j in 1:n_species
            for i in 1:n_species
                coll_factor_array[i,j,k] = CollisionFactorsSWPM()
            end
        end
    end
    return coll_factor_array
end

"""
    create_collision_factors_swpm_array(pia)

Create a 3-dimensional array of SWPM collision factors for all interaction pairs for all cells
in the simulation, with shape `(n_species,n_species,n_cells)`.

# Positional arguments
* `pia`: the ParticleIndexerArray instance

# Returns
3-dimensional array of `CollisionFactorsSWPM` instances with shape `(n_species,n_species,n_cells)`.
"""
function create_collision_factors_swpm_array(pia::ParticleIndexerArray)
    return create_collision_factors_swpm_array(size(pia.indexer)[2], size(pia.indexer)[1])
end

"""
    create_collision_factors_swpm_array(pia, interactions, species_data, T_list; mult_factor=1.0)

Create a 3-dimensional array of SWPM collision factors for all interaction pairs for all cells
in the simulation, with shape `(n_species,n_species,n_cells)`.
This will fill the array with the estimates ``(\\sigma g)_{max}`` for all species in all cells, assuming
a VHS cross-section, and that the temperature of each
species is constant across all cells.

# Positional arguments
* `pia`: the ParticleIndexerArray instance
* `interactions`: the 2-dimensional array of `Interaction` instances (of shape `(n_species, n_species)`) of all the pair-wise interactions
* `species_data`: the vector of `Species` instances of the species in the flow 
* `T_list`: the list of temperatures of the species

# Keyword arguments
* `mult_factor`: a factor by which to multiply the result (default value is 1.0)

# Returns
3-dimensional array of `CollisionFactorsSWPM` instances with shape `(n_species,n_species,n_cells)` filled with estimated
values of ``(\\sigma g)_{max}``.
"""
function create_collision_factors_swpm_array(pia, interactions, species_data, T_list; mult_factor=1.0)
    coll_factor_array = create_collision_factors_swpm_array(pia)

    estimate_sigma_g_max!(coll_factor_array, interactions, species_data, T_list; mult_factor=mult_factor)
    return coll_factor_array
end

"""
    create_collision_factors_swpm_array(pia, interactions, species_data, T::Real; mult_factor=1.0)

Create a 3-dimensional array of SWPM collision factors for all interaction pairs for all cells
in the simulation, with shape `(n_species,n_species,n_cells)`.
This will fill the array with the estimates ``(\\sigma g)_{max}`` for all species in all cells, assuming
a VHS cross-section, and that all species have a single temperature
that is constant across all cells.

# Positional arguments
* `pia`: the ParticleIndexerArray instance
* `interactions`: the 2-dimensional array of `Interaction` instances (of shape `(n_species, n_species)`) of all the pair-wise interactions
* `species_data`: the vector of `Species` instances of the species in the flow 
* `T`: the temperatures of the flow

# Keyword arguments
* `mult_factor`: a factor by which to multiply the result (default value is 1.0)

# Returns
3-dimensional array of `CollisionFactorsSWPM` instances with shape `(n_species,n_species,n_cells)` filled with estimated
values of ``(\\sigma g)_{max}``.
"""
function create_collision_factors_swpm_array(pia, interactions, species_data, T::Real; mult_factor=1.0)
    coll_factor_array = create_collision_factors_swpm_array(pia)
    estimate_sigma_g_max!(coll_factor_array, interactions, species_data, repeat([T], length(species_data)); mult_factor=mult_factor)
    return coll_factor_array
end

"""
    compute_n_coll_single_species(rng, collision_factors_swpm, np, w_max, G, Δt, V)

Compute the non-integer number of collisions between particles of same species for use in the SWPM method.

# Positional arguments
* `rng`: the random number generator
* `collision_factors`: the `CollisionFactorsSWPM` instance holding the estimate of ``(\\sigma g)_{max}``
    for the species in question in the cell
* `np`: number of particles in the cell
* `w_max`: the maximum computational weight of the particles in the cell
* `G`: non-negative value defining the weight transfer function
* `Δt`: timestep
* `V`: cell volume

# Returns
The non-integer number of collisions.
"""
function compute_n_coll_single_species(rng, collision_factors_swpm::CollisionFactorsSWPM, np, w_max, G, Δt, V)
    return 0.5 * Δt * np * (np - 1) * collision_factors_swpm.sigma_g_max * w_max * (G+1) / V +
        rand(rng, Float64)
end

"""
    swpm!(rng, collision_factors_swpm, collision_data, interaction, particles, pia,
          cell, species, G, Δt, V)

Perform elastic collisions between variable-weight particles of same species using the SWPM algorithm
and the VHS cross-section model.
During a collision of particles with weights ``w_i``, ``w_j``, the weights are depleted by
``\\min(w_i, w_j) / (1+G)``, where ``G \\geq 0`` is a user-defined parameter.

# Positional arguments
* `rng`: the random number generator
* `collision_factors_swpm`: the `CollisionFactorsSWPM` for the species in question in the cell
* `collision_data`: `CollisionData` instance used for storing collisional quantities
* `interaction`: 2-dimensional array of `Interaction` instances for all possible species pairs
* `particles`: `ParticleVector` of the particles being collided
* `pia`: the `ParticleIndexerArray`
* `cell`: the index of the cell in which collisions are performed
* `species`: the index of the species for which collisions are performed
* `G`: non-negative value defining the weight transfer function
* `Δt`: timestep
* `V`: cell volume

# References
* S. Rjasanow, W.Wagner, Stochastic numerics for the Boltzmann equation.
    [Springer Berlin, Heidelberg, 2005](https://doi.org/10.1007/3-540-27689-0).
"""
function swpm!(rng, collision_factors_swpm, collision_data, interaction, particles, pia,
               cell, species, G, Δt, V)
    # single-species swpm
    # find w_max
    w_max = 0.0
    @inbounds s1 = pia.indexer[cell, species].start1
    @inbounds e1 = pia.indexer[cell, species].end1
    @inbounds for i in s1:e1
        w_max = max(w_max, particles[i].w)
    end

    if pia.indexer[cell, species].n_group2 > 0
        @inbounds s2 = pia.indexer[cell, species].start2
        @inbounds e2 = pia.indexer[cell, species].end2
        @inbounds for i in s2:e2
            w_max = max(w_max, particles[i].w)
        end
    end

    wtf = 1.0 / (1.0 + G)
    inv_w_max = 1.0 / w_max

    # compute ncoll
    # loop over particles
    # update sigma_g_w_max
    # collide
    @inbounds collision_factors_swpm.n1 = pia.indexer[cell, species].n_local
    @inbounds collision_factors_swpm.n2 = pia.indexer[cell, species].n_local
    @inbounds n_coll_float = compute_n_coll_single_species(rng, collision_factors_swpm, pia.indexer[cell, species].n_local,
                                                           w_max, G, Δt, V)
    n_coll_int = floor(Int64, n_coll_float)

    collision_factors_swpm.n_coll = n_coll_int
    collision_factors_swpm.n_coll_performed = 0

    for _ in 1:n_coll_int
        # TODO: check if 0 or 1 !!!
        @inbounds i = floor(Int64, rand(rng, Float64) * pia.indexer[cell, species].n_local)
        @inbounds k = floor(Int64, rand(rng, Float64) * pia.indexer[cell, species].n_local)

        while (i == k)
            @inbounds k = floor(Int64, rand(rng, Float64) * pia.indexer[cell, species].n_local)
        end

        # example: bounds from [1,4], [7,9]; n_total = 7
        # n_group1 = 4, n_group2 = 3
        # i = 0,1,2,3 - [1,4]
        # i = 4,5,6 - [7,9]
        @inbounds i = map_cont_index(pia.indexer[cell, species], i)
        @inbounds k = map_cont_index(pia.indexer[cell, species], k)
        
        @inbounds compute_g!(collision_data, particles[i], particles[k])
        # println("NTC: ", particle_indexer.n_total, ", ", length(particles))
        if (collision_data.g > eps())

            @inbounds sigma = sigma_vhs(interaction[species, species], collision_data.g)
            @inbounds sigma_g_max = sigma * collision_data.g

            # update (σ g)_max if needed
            collision_factors_swpm.sigma_g_max = max(sigma_g_max, collision_factors_swpm.sigma_g_max)
            
            if (rand(rng, Float64) < sigma_g_max * max(particles[i].w, particles[k].w) * inv_w_max / collision_factors_swpm.sigma_g_max)
                collision_factors_swpm.n_coll_performed += 1
                @inbounds compute_com!(collision_data, interaction[species, species], particles[i], particles[k])
                # do collision
                if (length(particles) <= pia.n_total[species])
                    resize!(particles, length(particles)+DELTA_PARTICLES)
                end

                @inbounds Δw = min(particles[i].w, particles[k].w) * wtf

                @inbounds particles[i].w -= Δw
                @inbounds particles[k].w -= Δw

                update_buffer_index_new_particle!(particles, pia, cell, species)
                @inbounds particles[pia.n_total[species]] = Particle(Δw, particles[i].v, particles[i].x)
                update_buffer_index_new_particle!(particles, pia, cell, species)
                @inbounds particles[pia.n_total[species]] = Particle(Δw, particles[k].v, particles[k].x)
                @inbounds scatter_vhs!(rng, collision_data, interaction[species, species],
                                       particles[pia.n_total[species]-1], particles[pia.n_total[species]])
            end
        end
    end
end

end