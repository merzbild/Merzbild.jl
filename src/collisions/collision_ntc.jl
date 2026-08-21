@muladd begin

"""
    CollisionFactors

Structure to store NTC-related collision factors for collisions between particles of two species
in a given cell.

# Fields
* `n1`: the number of particles of the first species in the cell
* `n2`: the number of particles of the second species in the cell
* `sigma_g_w_max`: estimate of the ``(\\sigma g w)_{max}`` (``\\sigma`` is the total collision cross-section,
    ``g`` is the relative collision velocity, ``w`` is the computational weight of the particles)
* `n_coll`: number of collisions to be tested
* `n_coll_performed`: number of collisions actually performed
* `n_eq_w_coll_performed`: number of collisions between particles with equal weights actually performed
"""
mutable struct CollisionFactors
    n1::Int64
    n2::Int64
    sigma_g_w_max::Float64  # (σ(g) * g * w)_max
    n_coll::Int64  # number of collisions tested
    n_coll_performed::Int64  # number of collisions performed
    n_eq_w_coll_performed::Int64  # number of collisions of particles with equal weights (no splitting), for debugging/etc
end

"""
    CollisionFactors()

Create an empty `CollisionFactors` instance (all values set to 0).
"""
CollisionFactors() = CollisionFactors(0, 0, 0.0, 0, 0, 0)

"""
    create_collision_factors_array(n_species)

Create a 3-dimensional array of collision factors for all interaction pairs for a 0-D case (1 spatial cell),
with shape `(n_species,n_species,1)`.

# Positional arguments
* `n_species`: number of species in the flow

# Returns
3-dimensional array of `CollisionFactors` instances with shape `(n_species,n_species,1)`.
"""
function create_collision_factors_array(n_species)
    coll_factor_array = Array{CollisionFactors, 3}(undef, (n_species, n_species, 1))
    for k in 1:n_species
        for i in 1:n_species
            coll_factor_array[i,k,1] = CollisionFactors()
        end
    end
    return coll_factor_array
end

"""
    create_collision_factors_array(n_species, n_cells)

Create a 3-dimensional array of collision factors for all interaction pairs for all cells
in the simulation, with shape `(n_species,n_species,n_cells)`.

# Positional arguments
* `n_species`: number of species in the flow
* `n_cells`: number of cells in the simulation

# Returns
3-dimensional array of `CollisionFactors` instances with shape `(n_species,n_species,n_cells)`.
"""
function create_collision_factors_array(n_species, n_cells)
    coll_factor_array = Array{CollisionFactors, 3}(undef, (n_species, n_species, n_cells))
    for k in 1:n_cells
        for j in 1:n_species
            for i in 1:n_species
                coll_factor_array[i,j,k] = CollisionFactors()
            end
        end
    end
    return coll_factor_array
end

"""
    create_collision_factors_array(pia)

Create a 3-dimensional array of collision factors for all interaction pairs for all cells
in the simulation, with shape `(n_species,n_species,n_cells)`.

# Positional arguments
* `pia`: the ParticleIndexerArray instance

# Returns
3-dimensional array of `CollisionFactors` instances with shape `(n_species,n_species,n_cells)`.
"""
function create_collision_factors_array(pia::ParticleIndexerArray)
    return create_collision_factors_array(pia.n_species, pia.n_cells)
end

"""
    create_collision_factors_array(pia, interactions, species_data, T_list, Fnum::Real; mult_factor=1.0)

Create a 3-dimensional array of collision factors for all interaction pairs for all cells
in the simulation, with shape `(n_species,n_species,n_cells)`.
This will fill the array with the estimates ``(\\sigma g w)_{max}`` for all species in all cells, assuming
a constant particle computational weight `Fnum`, a VHS cross-section, and that the temperature of each
species is constant across all cells.

# Positional arguments
* `pia`: the ParticleIndexerArray instance
* `interactions`: the 2-dimensional array of `Interaction` instances (of shape `(n_species, n_species)`) of all the pair-wise interactions
* `species_data`: the vector of `Species` instances of the species in the flow 
* `T_list`: the list of temperatures of the species
* `Fnum`: the constant computational weight of the particles

# Keyword arguments
* `mult_factor`: a factor by which to multiply the result (default value is 1.0)

# Returns
3-dimensional array of `CollisionFactors` instances with shape `(n_species,n_species,n_cells)` filled with estimated
values of ``(\\sigma g w)_{max}``.
"""
function create_collision_factors_array(pia, interactions, species_data, T_list, Fnum::Real; mult_factor=1.0)
    coll_factor_array = create_collision_factors_array(pia)
    estimate_sigma_g_w_max!(coll_factor_array, interactions, species_data, T_list,
                            Fnum; mult_factor=mult_factor)
    return coll_factor_array
end

"""
    create_collision_factors_array(pia, interactions, species_data, T::Real, Fnum::Real; mult_factor=1.0)

Create a 3-dimensional array of collision factors for all interaction pairs for all cells
in the simulation, with shape `(n_species,n_species,n_cells)`.
This will fill the array with the estimates ``(\\sigma g w)_{max}`` for all species in all cells, assuming
a constant particle computational weight `Fnum`, a VHS cross-section, and that all species have a single temperature
that is constant across all cells.

# Positional arguments
* `pia`: the ParticleIndexerArray instance
* `interactions`: the 2-dimensional array of `Interaction` instances (of shape `(n_species, n_species)`) of all the pair-wise interactions
* `species_data`: the vector of `Species` instances of the species in the flow 
* `T`: the temperatures of the flow
* `Fnum`: the constant computational weight of the particles

# Keyword arguments
* `mult_factor`: a factor by which to multiply the result (default value is 1.0)

# Returns
3-dimensional array of `CollisionFactors` instances with shape `(n_species,n_species,n_cells)` filled with estimated
values of ``(\\sigma g w)_{max}``.
"""
function create_collision_factors_array(pia, interactions, species_data, T::Real, Fnum::Real; mult_factor=1.0)
    coll_factor_array = create_collision_factors_array(pia)
    estimate_sigma_g_w_max!(coll_factor_array, interactions, species_data, repeat([T], length(species_data)),
                            Fnum; mult_factor=mult_factor)
    return coll_factor_array
end

"""
    compute_n_coll_single_species(rng, collision_factors, np, Δt, V)

Compute the non-integer number of collisions between particles of same species.

# Positional arguments
* `rng`: the random number generator
* `collision_factors`: the `CollisionFactors` holding the estimate of ``(\\sigma g w)_{max}``
    for the species in question in the cell
* `np`: number of particles in the cell
* `Δt`: timestep
* `V`: cell volume

# Returns
The non-integer number of collisions.
"""
function compute_n_coll_single_species(rng, collision_factors, np, Δt, V)
    return 0.5 * Δt * np * (np - 1) * collision_factors.sigma_g_w_max / V +
        rand(rng, Float64)
end

"""
    compute_n_coll_two_species(rng, collision_factors, np1, np2, Δt, V)

Compute number of collisions between particles of different species

# Positional arguments
* `rng`: the random number generator
* `collision_factors`: the `CollisionFactors` holding the estimate of ``(\\sigma g w)_{max}``
    for the species in question in the cell
* `np1`: number of particles of the first species in the cell
* `np2`: number of particles of the second species in the cell
* `Δt`: timestep
* `V`: cell volume

# Returns
The non-integer number of collisions.
"""
function compute_n_coll_two_species(rng, collision_factors, np1, np2, Δt, V)
    return Δt * np1 * np2 * collision_factors.sigma_g_w_max / V + rand(rng, Float64)
end

"""
    collide_2particles!(rng, model, collision_data, collision_factors, interaction, pa_i::Particle{D}, pa_k::Particle{D},
                        particles_1::ParticleVector{D}, particles_2::ParticleVector{D}, pia, cell, species1, species2; dw_tol=1e-16) where D

Collide two particles elastically using the elastic scattering model `model`.
Particles can be of same or different species.
If particles' weights differ by less than `dw_tol`, an equal-weight collision is performed and no particles are split.

# Positional arguments
* `rng`: the random number generator
* `model`: the `AbstractScatteringModel` singleton tag of the scattering model
* `collision_data`: `CollisionData` instance used for storing collisional quantities
* `collision_factors`: the `CollisionFactors` holding the estimate of ``(\\sigma g w)_{max}``
    for the species in question in the cell
* `interaction`: the `Interaction` instance for the colliding species
* `particles_1`: `ParticleVector` of the particles of the first species being collided
* `particles_2`: `ParticleVector` of the particles of the second species being collided
* `pia`: the `ParticleIndexerArray`
* `cell`: the index of the cell in which collisions are performed
* `species1`: the index of the first species for which collisions are performed
* `species2`: the index of the second species for which collisions are performed

# Keyword arguments
* `dw_tol`: if weights of particles differ by less than this amount, an equal-weight collision is assumed
and no particle splitting is performed
"""
@inline function collide_2particles!(rng, model::AbstractScatteringModel, collision_data, collision_factors, interaction,
                                     pa_i::Particle{D}, pa_k::Particle{D},
                                     particles_1::ParticleVector{D}, particles_2::ParticleVector{D}, pia, cell, species1, species2; dw_tol=1e-16) where D
    sigma_coll = sigma(model, interaction, collision_data.g)
    sigma_g_w_max = sigma_coll * collision_data.g * max(pa_i.w, pa_k.w)

    # update (σ g w)_max if needed
    collision_factors.sigma_g_w_max = max(sigma_g_w_max, collision_factors.sigma_g_w_max)

    @inbounds if (rand(rng, Float64) * collision_factors.sigma_g_w_max < sigma_g_w_max)
        collision_factors.n_coll_performed += 1
        compute_com!(collision_data, interaction, pa_i, pa_k)
        # do collision
        if (abs(pa_i.w - pa_k.w) < dw_tol)
            collision_factors.n_eq_w_coll_performed += 1
        elseif (pa_i.w > pa_k.w)
            # we split particle i, update velocity of i and k (split part remains unchanged)

            # first need to grow particle array
            if (length(particles_1) <= pia.index_last[species1])
                resize!(particles_1, length(particles_1)+DELTA_PARTICLES)
            end

            # first need to update the particle indexer struct
            update_buffer_index_new_particle!(particles_1, pia, cell, species1)

            Δw = pa_i.w - pa_k.w
            pa_i.w = pa_k.w

            p1_new = particles_1[pia.index_last[species1]]

            p1_new.w = Δw
            p1_new.v = pa_i.v
            p1_new.x = pa_i.x
        else  # (particles[k].w > particles[i].w)
            if (length(particles_2) <= pia.index_last[species2])
                resize!(particles_2, length(particles_2)+DELTA_PARTICLES)
            end

            update_buffer_index_new_particle!(particles_2, pia, cell, species2)

            Δw = pa_k.w - pa_i.w
            pa_k.w = pa_i.w

            p2_new = particles_2[pia.index_last[species2]]

            p2_new.w = Δw
            p2_new.v = pa_k.v
            p2_new.x = pa_k.x
        end
        scatter!(rng, model, collision_data, interaction, pa_i, pa_k)
    end
end

"""
    collide_2particles_equal_weight!(rng, model, collision_data, collision_factors, interaction, pa_i::Particle{D}, pa_k::Particle{D}) where D

Collide two particles elastically using the elastic scattering model `model`, assuming equal weights -
no particle splitting is performed even if weights are unequal.
Particles can be of same or different species.

# Positional arguments
* `rng`: the random number generator
* `model`: the `AbstractScatteringModel` singleton tag of the scattering model
* `collision_data`: `CollisionData` instance used for storing collisional quantities
* `collision_factors`: the `CollisionFactors` holding the estimate of ``(\\sigma g w)_{max}``
    for the species in question in the cell
* `interaction`: the `Interaction` instance for the colliding species
* `pa_i`: the first particle being collided
* `pa_k`: the second particle being collided
"""
@inline function collide_2particles_equal_weight!(rng, model::AbstractScatteringModel, collision_data, collision_factors, interaction,
                                                  pa_i::Particle{D}, pa_k::Particle{D}) where D
    sigma_coll = sigma(model, interaction, collision_data.g)
    sigma_g_w_max = sigma_coll * collision_data.g * max(pa_i.w, pa_k.w)

    # update (σ g w)_max if needed
    collision_factors.sigma_g_w_max = max(sigma_g_w_max, collision_factors.sigma_g_w_max)

    @inbounds if (rand(rng, Float64) * collision_factors.sigma_g_w_max < sigma_g_w_max)
        collision_factors.n_coll_performed += 1
        collision_factors.n_eq_w_coll_performed += 1
        compute_com!(collision_data, interaction, pa_i, pa_k)
        # do collision
        scatter!(rng, model, collision_data, interaction, pa_i, pa_k)
    end
end

"""
    ntc!(rng, collision_factors, collision_data, interaction, particles::ParticleVector{D}, pia,
         cell, species, Δt, V; dw_tol=1e-16) where D

Perform elastic collisions between particles of same species using the NTC algorithm.
The elastic scattering model is taken from the `Interaction` instance of the species pair.

# Positional arguments
* `rng`: the random number generator
* `collision_factors`: the `CollisionFactors` for the species in question in the cell
* `collision_data`: `CollisionData` instance used for storing collisional quantities
* `interaction`: 2-dimensional array of `Interaction` instances for all possible species pairs
* `particles`: `ParticleVector` of the particles being collided
* `pia`: the `ParticleIndexerArray`
* `cell`: the index of the cell in which collisions are performed
* `species`: the index of the species for which collisions are performed
* `Δt`: timestep
* `V`: cell volume

# Keyword arguments
* `dw_tol`: if weights of particles differ by less than this amount, an equal-weight collision is assumed
and no particle splitting is performed

# References
* D.P. Schmidt, C.J. Rutland, A New Droplet Collision Algorithm.
    [J. Comput. Phys, 2000](https://doi.org/10.1006/jcph.2000.6568).
"""
function ntc!(rng, collision_factors, collision_data, interaction, particles::ParticleVector{D}, pia,
              cell, species, Δt, V; dw_tol=1e-16) where D
    @inbounds model = interaction[species, species].model

    @scattering_barrier model ntc!(rng, collision_factors, collision_data, interaction, particles, pia,
                                   cell, species, Δt, V; dw_tol=dw_tol)
end

"""
    ntc!(rng, model::AbstractScatteringModel, collision_factors, collision_data, interaction,
         particles::ParticleVector{D}, pia, cell, species, Δt, V; dw_tol=1e-16) where D

Perform elastic collisions between particles of same species using the NTC algorithm
and the elastic scattering model `model`, overriding the model stored in the `Interaction` instance
of the species pair.

# Positional arguments
* `rng`: the random number generator
* `model`: the `AbstractScatteringModel` singleton tag of the scattering model
* `collision_factors`: the `CollisionFactors` for the species in question in the cell
* `collision_data`: `CollisionData` instance used for storing collisional quantities
* `interaction`: 2-dimensional array of `Interaction` instances for all possible species pairs
* `particles`: `ParticleVector` of the particles being collided
* `pia`: the `ParticleIndexerArray`
* `cell`: the index of the cell in which collisions are performed
* `species`: the index of the species for which collisions are performed
* `Δt`: timestep
* `V`: cell volume

# Keyword arguments
* `dw_tol`: if weights of particles differ by less than this amount, an equal-weight collision is assumed
and no particle splitting is performed

# References
* D.P. Schmidt, C.J. Rutland, A New Droplet Collision Algorithm.
    [J. Comput. Phys, 2000](https://doi.org/10.1006/jcph.2000.6568).
"""
function ntc!(rng, model::AbstractScatteringModel, collision_factors, collision_data, interaction,
              particles::ParticleVector{D}, pia,
              cell, species, Δt, V; dw_tol=1e-16) where D
    # single-species ntc
    # compute ncoll
    # loop over particles
    # update sigma_g_w_max
    # collide

    @inbounds indexer = pia.indexer[cell, species]

    @inbounds collision_factors.n1 = indexer.n_local
    @inbounds collision_factors.n2 = indexer.n_local
    @inbounds n_coll_float = compute_n_coll_single_species(rng, collision_factors, indexer.n_local, Δt, V)
    n_coll_int = trunc(Int64, n_coll_float)
    # println(n_coll_float, ", ", n_coll_int)

    collision_factors.n_coll = n_coll_int
    collision_factors.n_coll_performed = 0
    collision_factors.n_eq_w_coll_performed = 0

    @inbounds interaction_l = interaction[species, species]

    @inbounds for _ in 1:n_coll_int

        n_loc = indexer.n_local  # can change due to splitting!
        i = trunc(Int64, rand(rng, Float64) * n_loc)
        k = trunc(Int64, rand(rng, Float64) * n_loc)

        while (i == k)
            k = trunc(Int64, rand(rng, Float64) * n_loc)
        end

        # example: bounds from [1,4], [7,9]; n_total = 7
        # n_group1 = 4, n_group2 = 3
        # i = 0,1,2,3 - [1,4]
        # i = 4,5,6 - [7,9]
        i = map_cont_index(indexer, i)
        k = map_cont_index(indexer, k)
        pa_i = particles[i]
        pa_k = particles[k]
        
        compute_g!(collision_data, pa_i, pa_k)
        if (collision_data.g > eps())
            collide_2particles!(rng, model, collision_data, collision_factors, interaction_l, pa_i, pa_k,
                                particles, particles, pia, cell, species, species; dw_tol=dw_tol)
        end
    end
end

"""
    ntc!(rng, collision_factors, collision_data, interaction,
         particles_1::ParticleVector{D}, particles_2::ParticleVector{D}, pia,
         cell, species1, species2, Δt, V; dw_tol=1e-16) where D

Perform elastic collisions between particles of different species using the NTC algorithm.
The elastic scattering model is taken from the `Interaction` instance of the species pair.

# Positional arguments
* `rng`: the random number generator
* `collision_factors`: the `CollisionFactors` for the species in question in the cell
* `collision_data`: `CollisionData` instance used for storing collisional quantities
* `interaction`: 2-dimensional array of `Interaction` instances for all possible species pairs
* `particles_1`: `ParticleVector` of the particles of the first species being collided
* `particles_2`: `ParticleVector` of the particles of the second species being collided
* `pia`: the `ParticleIndexerArray`
* `cell`: the index of the cell in which collisions are performed
* `species1`: the index of the first species for which collisions are performed
* `species2`: the index of the second species for which collisions are performed
* `Δt`: timestep
* `V`: cell volume

# Keyword arguments
* `dw_tol`: if weights of particles differ by less than this amount, an equal-weight collision is assumed
and no particle splitting is performed

# References
* D.P. Schmidt, C.J. Rutland, A New Droplet Collision Algorithm.
    [J. Comput. Phys, 2000](https://doi.org/10.1006/jcph.2000.6568).
"""
function ntc!(rng, collision_factors, collision_data, interaction,
              particles_1::ParticleVector{D}, particles_2::ParticleVector{D}, pia,
              cell, species1, species2, Δt, V; dw_tol=1e-16) where D
    @inbounds model = interaction[species1, species2].model

    @scattering_barrier model ntc!(rng, collision_factors, collision_data, interaction,
                                   particles_1, particles_2, pia,
                                   cell, species1, species2, Δt, V; dw_tol=dw_tol)
end

"""
    ntc!(rng, model::AbstractScatteringModel, collision_factors, collision_data, interaction,
         particles_1::ParticleVector{D}, particles_2::ParticleVector{D}, pia,
         cell, species1, species2, Δt, V; dw_tol=1e-16) where D

Perform elastic collisions between particles of different species using the NTC algorithm
and the elastic scattering model `model`, overriding the model stored in the `Interaction` instance
of the species pair.

# Positional arguments
* `rng`: the random number generator
* `model`: the `AbstractScatteringModel` singleton tag of the scattering model
* `collision_factors`: the `CollisionFactors` for the species in question in the cell
* `collision_data`: `CollisionData` instance used for storing collisional quantities
* `interaction`: 2-dimensional array of `Interaction` instances for all possible species pairs
* `particles_1`: `ParticleVector` of the particles of the first species being collided
* `particles_2`: `ParticleVector` of the particles of the second species being collided
* `pia`: the `ParticleIndexerArray`
* `cell`: the index of the cell in which collisions are performed
* `species1`: the index of the first species for which collisions are performed
* `species2`: the index of the second species for which collisions are performed
* `Δt`: timestep
* `V`: cell volume

# Keyword arguments
* `dw_tol`: if weights of particles differ by less than this amount, an equal-weight collision is assumed
and no particle splitting is performed

# References
* D.P. Schmidt, C.J. Rutland, A New Droplet Collision Algorithm.
    [J. Comput. Phys, 2000](https://doi.org/10.1006/jcph.2000.6568).
"""
function ntc!(rng, model::AbstractScatteringModel, collision_factors, collision_data, interaction,
              particles_1::ParticleVector{D}, particles_2::ParticleVector{D}, pia,
              cell, species1, species2, Δt, V; dw_tol=1e-16) where D
    # compute ncoll
    # loop over particles
    # update sigma_g_w_max
    # collide

    @inbounds indexer1 = pia.indexer[cell, species1]
    @inbounds indexer2 = pia.indexer[cell, species2]

    @inbounds collision_factors.n1 = indexer1.n_local
    @inbounds collision_factors.n2 = indexer2.n_local
    @inbounds n_coll_float = compute_n_coll_two_species(rng, collision_factors,
                                                        indexer1.n_local, indexer2.n_local, Δt, V)
    n_coll_int = trunc(Int64, n_coll_float)
    # println(n_coll_float, ", ", n_coll_int)

    collision_factors.n_coll = n_coll_int
    collision_factors.n_coll_performed = 0
    collision_factors.n_eq_w_coll_performed = 0

    @inbounds interaction_l = interaction[species1, species2]

    @inbounds for _ in 1:n_coll_int

        i = trunc(Int64, rand(rng, Float64) * indexer1.n_local)
        k = trunc(Int64, rand(rng, Float64) * indexer2.n_local)

        # example: bounds from [1,4], [7,9]; n_total = 7
        # n_group1 = 4, n_group2 = 3
        # i = 0,1,2,3 - [1,4]
        # i = 4,5,6 - [7,9]
        i = map_cont_index(indexer1, i)
        k = map_cont_index(indexer2, k)

        pa_i = particles_1[i]
        pa_k = particles_2[k]
        
        compute_g!(collision_data, pa_i, pa_k)

        if (collision_data.g > eps())
            collide_2particles!(rng, model, collision_data, collision_factors, interaction_l, pa_i, pa_k,
                                particles_1, particles_2, pia, cell, species1, species2; dw_tol=dw_tol)
        end
    end
end

"""
    ntc_equal_weight!(rng, collision_factors, collision_data, interaction, particles::ParticleVector{D}, pia,
                      cell, species, Δt, V) where D

Perform elastic collisions between particles of same species using the NTC algorithm.
Particle weights are assumed to be equal, and no weight checks/splitting is performed.
The elastic scattering model is taken from the `Interaction` instance of the species pair.

# Positional arguments
* `rng`: the random number generator
* `collision_factors`: the `CollisionFactors` for the species in question in the cell
* `collision_data`: `CollisionData` instance used for storing collisional quantities
* `interaction`: 2-dimensional array of `Interaction` instances for all possible species pairs
* `particles`: `ParticleVector` of the particles being collided
* `pia`: the `ParticleIndexerArray`
* `cell`: the index of the cell in which collisions are performed
* `species`: the index of the species for which collisions are performed
* `Δt`: timestep
* `V`: cell volume

# References
* G.A. Bird, Molecular gas dynamics and the direct simulation of gas flows,
    [Clarendon Press, Oxford, 1994](https://doi.org/10.1093/oso/9780198561958.001.0001).
"""
function ntc_equal_weight!(rng, collision_factors, collision_data, interaction, particles::ParticleVector{D}, pia,
                           cell, species, Δt, V) where D
    @inbounds model = interaction[species, species].model

    @scattering_barrier model ntc_equal_weight!(rng, collision_factors, collision_data, interaction, particles, pia,
                                                 cell, species, Δt, V)
end

"""
    ntc_equal_weight!(rng, model::AbstractScatteringModel, collision_factors, collision_data, interaction,
                      particles::ParticleVector{D}, pia, cell, species, Δt, V) where D

Perform elastic collisions between particles of same species using the NTC algorithm
and the elastic scattering model `model`, overriding the model stored in the `Interaction` instance
of the species pair. Particle weights are assumed to be equal,
and no weight checks/splitting is performed.

# Positional arguments
* `rng`: the random number generator
* `model`: the `AbstractScatteringModel` singleton tag of the scattering model
* `collision_factors`: the `CollisionFactors` for the species in question in the cell
* `collision_data`: `CollisionData` instance used for storing collisional quantities
* `interaction`: 2-dimensional array of `Interaction` instances for all possible species pairs
* `particles`: `ParticleVector` of the particles being collided
* `pia`: the `ParticleIndexerArray`
* `cell`: the index of the cell in which collisions are performed
* `species`: the index of the species for which collisions are performed
* `Δt`: timestep
* `V`: cell volume

# References
* G.A. Bird, Molecular gas dynamics and the direct simulation of gas flows,
    [Clarendon Press, Oxford, 1994](https://doi.org/10.1093/oso/9780198561958.001.0001).
"""
function ntc_equal_weight!(rng, model::AbstractScatteringModel, collision_factors, collision_data, interaction,
                           particles::ParticleVector{D}, pia,
                           cell, species, Δt, V) where D
    # single-species ntc
    # compute ncoll
    # loop over particles
    # update sigma_g_w_max
    # collide
    @inbounds indexer = pia.indexer[cell, species]

    @inbounds collision_factors.n1 = indexer.n_local
    @inbounds collision_factors.n2 = indexer.n_local
    @inbounds n_coll_float = compute_n_coll_single_species(rng, collision_factors, indexer.n_local, Δt, V)
    n_coll_int = floor(Int64, n_coll_float)
    # println(n_coll_float, ", ", n_coll_int)

    collision_factors.n_coll = n_coll_int
    collision_factors.n_coll_performed = 0
    collision_factors.n_eq_w_coll_performed = 0

    @inbounds interaction_l = interaction[species, species]

    @inbounds for _ in 1:n_coll_int
        i = floor(Int64, rand(rng, Float64) * indexer.n_local)
        k = floor(Int64, rand(rng, Float64) * indexer.n_local)

        while (i == k)
            k = floor(Int64, rand(rng, Float64) * indexer.n_local)
        end

        # example: bounds from [1,4], [7,9]; n_total = 7
        # n_group1 = 4, n_group2 = 3
        # i = 0,1,2,3 - [1,4]
        # i = 4,5,6 - [7,9]
        i = map_cont_index(indexer, i)
        k = map_cont_index(indexer, k)
        pa_i = particles[i]
        pa_k = particles[k]
        
        compute_g!(collision_data, pa_i, pa_k)
        if (collision_data.g > eps())
            collide_2particles_equal_weight!(rng, model, collision_data, collision_factors, interaction_l, pa_i, pa_k)
        end
    end
end

"""
    ntc_equal_weight!(rng, collision_factors, collision_data, interaction,
                      particles_1, particles_2, pia,
                      cell, species1, species2, Δt, V)

Perform elastic collisions between particles of different species using the NTC algorithm.
Particle weights are assumed to be equal, and no weight checks/splitting is performed.
The elastic scattering model is taken from the `Interaction` instance of the species pair.

# Positional arguments
* `rng`: the random number generator
* `collision_factors`: the `CollisionFactors` for the species in question in the cell
* `collision_data`: `CollisionData` instance used for storing collisional quantities
* `interaction`: 2-dimensional array of `Interaction` instances for all possible species pairs
* `particles_1`: `ParticleVector` of the particles of the first species being collided
* `particles_2`: `ParticleVector` of the particles of the second species being collided
* `pia`: the `ParticleIndexerArray`
* `cell`: the index of the cell in which collisions are performed
* `species1`: the index of the first species for which collisions are performed
* `species2`: the index of the second species for which collisions are performed
* `Δt`: timestep
* `V`: cell volume

# References
* G.A. Bird, Molecular gas dynamics and the direct simulation of gas flows,
    [Clarendon Press, Oxford, 1994](https://doi.org/10.1093/oso/9780198561958.001.0001).
"""
function ntc_equal_weight!(rng, collision_factors, collision_data, interaction,
                           particles_1, particles_2, pia,
                           cell, species1, species2, Δt, V)
    @inbounds model = interaction[species1, species2].model

    @scattering_barrier model ntc_equal_weight!(rng, collision_factors, collision_data, interaction,
                                                 particles_1, particles_2, pia,
                                                 cell, species1, species2, Δt, V)
end

"""
    ntc_equal_weight!(rng, model::AbstractScatteringModel, collision_factors, collision_data, interaction,
                      particles_1, particles_2, pia,
                      cell, species1, species2, Δt, V)

Perform elastic collisions between particles of different species using the NTC algorithm
and the elastic scattering model `model`, overriding the model stored in the `Interaction` instance
of the species pair. Particle weights are assumed to be equal,
and no weight checks/splitting is performed.

# Positional arguments
* `rng`: the random number generator
* `model`: the `AbstractScatteringModel` singleton tag of the scattering model
* `collision_factors`: the `CollisionFactors` for the species in question in the cell
* `collision_data`: `CollisionData` instance used for storing collisional quantities
* `interaction`: 2-dimensional array of `Interaction` instances for all possible species pairs
* `particles_1`: `ParticleVector` of the particles of the first species being collided
* `particles_2`: `ParticleVector` of the particles of the second species being collided
* `pia`: the `ParticleIndexerArray`
* `cell`: the index of the cell in which collisions are performed
* `species1`: the index of the first species for which collisions are performed
* `species2`: the index of the second species for which collisions are performed
* `Δt`: timestep
* `V`: cell volume

# References
* G.A. Bird, Molecular gas dynamics and the direct simulation of gas flows,
    [Clarendon Press, Oxford, 1994](https://doi.org/10.1093/oso/9780198561958.001.0001).
"""
function ntc_equal_weight!(rng, model::AbstractScatteringModel, collision_factors, collision_data, interaction,
                           particles_1, particles_2, pia,
                           cell, species1, species2, Δt, V)
    # compute ncoll
    # loop over particles
    # update sigma_g_w_max
    # collide
    @inbounds indexer1 = pia.indexer[cell, species1]
    @inbounds indexer2 = pia.indexer[cell, species2]

    @inbounds collision_factors.n1 = indexer1.n_local
    @inbounds collision_factors.n2 = indexer2.n_local
    @inbounds n_coll_float = compute_n_coll_two_species(rng, collision_factors,
                                                        indexer1.n_local, indexer2.n_local, Δt, V)
    n_coll_int = floor(Int64, n_coll_float)
    # println(n_coll_float, ", ", n_coll_int)

    collision_factors.n_coll = n_coll_int
    collision_factors.n_coll_performed = 0
    collision_factors.n_eq_w_coll_performed = 0

    @inbounds interaction_l = interaction[species1, species2]

    @inbounds for _ in 1:n_coll_int
        i = floor(Int64, rand(rng, Float64) * indexer1.n_local)
        k = floor(Int64, rand(rng, Float64) * indexer2.n_local)

        # example: bounds from [1,4], [7,9]; n_total = 7
        # n_group1 = 4, n_group2 = 3
        # i = 0,1,2,3 - [1,4]
        # i = 4,5,6 - [7,9]
        i = map_cont_index(indexer1, i)
        k = map_cont_index(indexer2, k)

        pa_i = particles_1[i]
        pa_k = particles_2[k]
        
        compute_g!(collision_data, pa_i, pa_k)

        if (collision_data.g > eps())
            collide_2particles_equal_weight!(rng, model, collision_data, collision_factors, interaction_l, pa_i, pa_k)
        end
    end
end

"""
    estimate_sigma_g_w_max_ntc_n_e!(rng, collision_factors, collision_data, interaction,
                                    n_e_interactions, n_e_cs, particles_n, particles_e,
                                    pia, cell, species_n, species_e, Δt, V; min_coll=5, n_loops=3)

Estimate ``(\\sigma g w)_{max}`` for an electron-neutral interaction by stochastically choosing particle pairs
multiple times and computing ``(\\sigma g w)`` for each pair. The number of collisions is computed
using the standard variable-weight NTC formula, the value of `min_coll` is added to this number,
and particles are randomly sampled. The whole procedure is repeated `n_loops` times, so that
an increased value ``(\\sigma g w)_{max}`` can have an impact on the computed number of pairs to select
during the next loop iteration.

# Positional arguments
* `rng`: the random number generator
* `collision_factors`: the `CollisionFactors` for the species in question in the cell
* `collision_data`: `CollisionData` instance used for storing collisional quantities
* `interaction`: 2-dimensional array of `Interaction` instances for all possible species pairs
* `n_e_interactions`: the `ElectronNeutralInteractions` instance
* `n_e_cs`: a vector of `ComputedCrossSections` instances
* `particles_n`: `ParticleVector` of the particles of neutral species
* `particles_e`: `ParticleVector` of the particles of the electron species
* `pia`: the `ParticleIndexerArray`
* `cell`: the index of the cell in which collisions are performed
* `species_n`: the index of the neutral species
* `species_e`: the index of the electron species
* `Δt`: timestep
* `V`: cell volume

# Keyword arguments:
* `min_coll`: the minimum number of pairs to test
* `n_loops`: the number of loops to perform (in each loop the number of collisions is computed using the
    estimated value of  ``(\\sigma g w)_{max}``)
* `extend`: enum of `CSExtend` type that sets how out-of-range energy values are treated when computing cross-sections
"""
function estimate_sigma_g_w_max_ntc_n_e!(rng, collision_factors, collision_data, interaction,
                                         n_e_interactions, n_e_cs, particles_n, particles_e,
                                         pia, cell, species_n, species_e, Δt, V; min_coll=5, n_loops=3, extend::CSExtend=CSExtendConstant)
    # compute ncoll
    # loop over particles
    # update sigma_g_w_max
    # collide
    @inbounds indexer_n = pia.indexer[cell, species_n]
    @inbounds indexer_e = pia.indexer[cell, species_e]

    @inbounds collision_factors.n1 = indexer_n.n_local
    @inbounds collision_factors.n2 = indexer_e.n_local

    # need to index into n_e_cs
    @inbounds species_n_en_i = n_e_interactions.neutral_indexer[species_n]

    @inbounds interaction_l = interaction[species_n, species_e]

    @inbounds for _ in 1:n_loops
        n_coll_float = compute_n_coll_two_species(rng, collision_factors,
                                                  indexer_n.n_local,
                                                  indexer_e.n_local, Δt, V) + min_coll
        n_coll_int = floor(Int64, n_coll_float)

        collision_factors.n_coll = n_coll_int
        collision_factors.n_coll_performed = 0
        collision_factors.n_eq_w_coll_performed = 0

        for _ in 1:n_coll_int
            i = floor(Int64, rand(rng, Float64) * indexer_n.n_local)
            k = floor(Int64, rand(rng, Float64) * indexer_e.n_local)

            # example: bounds from [1,4], [7,9]; n_total = 7
            # n_group1 = 4, n_group2 = 3
            # i = 0,1,2,3 - [1,4]
            # i = 4,5,6 - [7,9]
            i = map_cont_index(indexer_n, i)
            k = map_cont_index(indexer_e, k)
            
            compute_g!(collision_data, particles_n[i], particles_e[k])

            if (collision_data.g > eps())
                collision_data.E_coll_eV = compute_cross_sections!(n_e_cs, interaction_l, collision_data.g, n_e_interactions, species_n;
                                                                             extend=extend)
                sigma_g_w_max = n_e_cs[species_n_en_i].cs_total * collision_data.g * max(particles_n[i].w, particles_e[k].w)

                # update (σ g w)_max if needed
                collision_factors.sigma_g_w_max = max(sigma_g_w_max, collision_factors.sigma_g_w_max)
            end
        end
    end
end



"""
    ntc_n_e!(rng, collision_factors, collision_data, interaction,
             n_e_interactions, n_e_cs, particles_n, particles_e, particles_ion,
             pia, cell, species_n, species_e, species_ion, Δt, V; extend::CSExtend=CSExtendConstant, dw_tol=1e-16)

Perform electron-neutral elastic scattering and electron-impact ionization collisions.

# Positional arguments
* `rng`: the random number generator
* `collision_factors`: the `CollisionFactors` for the species in question in the cell
* `collision_data`: `CollisionData` instance used for storing collisional quantities
* `interaction`: 2-dimensional array of `Interaction` instances for all possible species pairs
* `n_e_interactions`: the `ElectronNeutralInteractions` instance
* `n_e_cs`: a vector of `ComputedCrossSections` instances
* `particles_n`: `ParticleVector` of the particles of neutral species
* `particles_e`: `ParticleVector` of the particles of the electron species
* `particles_ion`: `ParticleVector` of the particles of the ion species
* `pia`: the `ParticleIndexerArray`
* `cell`: the index of the cell in which collisions are performed
* `species_n`: the index of the neutral species
* `species_e`: the index of the electron species
* `species_ion`: the index of the ion species
* `Δt`: timestep
* `V`: cell volume

# Keyword arguments
* `extend`: enum of `CSExtend` type that sets how out-of-range energy values are treated when computing cross-sections
* `dw_tol`: if weights of particles differ by less than this amount, an equal-weight collision is assumed
and no particle splitting is performed

# References
* D.P. Schmidt, C.J. Rutland, A New Droplet Collision Algorithm.
    [J. Comput. Phys, 2000](https://doi.org/10.1006/jcph.2000.6568).
"""
function ntc_n_e!(rng, collision_factors, collision_data, interaction,
                  n_e_interactions, n_e_cs, particles_n, particles_e, particles_ion,
                  pia, cell, species_n, species_e, species_ion, Δt, V; extend::CSExtend=CSExtendConstant, dw_tol=1e-16)
    # no event splitting
    # compute ncoll
    # loop over particles
    # update sigma_g_w_max
    # collide
    @inbounds indexer_n = pia.indexer[cell, species_n]
    @inbounds indexer_e = pia.indexer[cell, species_e]

    @inbounds collision_factors.n1 = indexer_n.n_local
    @inbounds collision_factors.n2 = pia.indexer[cell, species_e].n_local
    @inbounds n_coll_float = compute_n_coll_two_species(rng, collision_factors,
                                              indexer_n.n_local, pia.indexer[cell, species_e].n_local, Δt, V)

    @inbounds interaction_l = interaction[species_n, species_e]

    # need to index into n_e_cs
    @inbounds species_n_en_i = n_e_interactions.neutral_indexer[species_n]
    
    @inbounds mass_ratio = n_e_interactions.mass_ratios[species_n_en_i]

    n_coll_int = floor(Int64, n_coll_float)
    # println(n_coll_float, ", ", n_coll_int)

    collision_factors.n_coll = n_coll_int
    collision_factors.n_coll_performed = 0
    collision_factors.n_eq_w_coll_performed = 0

    @inbounds for _ in 1:n_coll_int
        i = floor(Int64, rand(rng, Float64) * indexer_n.n_local)
        k = floor(Int64, rand(rng, Float64) * indexer_e.n_local)

        # example: bounds from [1,4], [7,9]; n_total = 7
        # n_group1 = 4, n_group2 = 3
        # i = 0,1,2,3 - [1,4]
        # i = 4,5,6 - [7,9]
        i = map_cont_index(indexer_n, i)
        k = map_cont_index(indexer_e, k)

        particles_n_i = particles_n[i]
        particles_e_k = particles_e[k]
        
        compute_g!(collision_data, particles_n_i, particles_e_k)

        if (collision_data.g > eps())

            collision_data.E_coll_eV = compute_cross_sections!(n_e_cs, interaction_l,
                                                               collision_data.g, n_e_interactions, species_n;
                                                               extend=extend)
            sigma_g_w_max = get_cs_total(n_e_interactions, n_e_cs, species_n) * collision_data.g * max(particles_n_i.w, particles_e_k.w)

            # update (σ g w)_max if needed
            collision_factors.sigma_g_w_max = max(sigma_g_w_max, collision_factors.sigma_g_w_max)

            if (rand(rng, Float64) * collision_factors.sigma_g_w_max < sigma_g_w_max)
                collision_factors.n_coll_performed += 1
                compute_com!(collision_data, interaction_l, particles_n_i, particles_e_k)

                # do collision

                if (particles_n_i.w > particles_e_k.w)
                    # we split particle i, update velocity of i and k (split part remains unchanged)
                    if (length(particles_n) <= pia.index_last[species_n])
                        resize!(particles_n, length(particles_n)+DELTA_PARTICLES)
                    end
                    update_buffer_index_new_particle!(particles_n, pia, cell, species_n)

                    Δw = particles_n_i.w - particles_e_k.w
                    particles_n_i.w = particles_e_k.w

                    p_n_new = particles_n[pia.index_last[species_n]]

                    p_n_new.w = Δw
                    p_n_new.v = particles_n_i.v
                    p_n_new.x = particles_n_i.x
                elseif (abs(particles_n_i.w - particles_e_k.w) < dw_tol)
                    collision_factors.n_eq_w_coll_performed += 1
                else  # (particles[k].w > particles[i].w)
                    if (length(particles_e) <= pia.index_last[species_e])
                        resize!(particles_e, length(particles_e)+DELTA_PARTICLES)
                    end
                    update_buffer_index_new_particle!(particles_e, pia, cell, species_e)

                    Δw = particles_e_k.w - particles_n_i.w
                    particles_e_k.w = particles_n_i.w

                    p_e_new = particles_e[pia.index_last[species_e]]

                    p_e_new.w = Δw
                    p_e_new.v = particles_e_k.v
                    p_e_new.x = particles_e_k.x
                end

                # now we collide the 2 equal-weight particles
                R = rand(rng, Float64)  # decide which process we do

                if (R < n_e_cs[species_n_en_i].prob_vec[1]) 
                    # elastic collision
                    scatter_vhs!(rng, collision_data, interaction_l, particles_n_i, particles_e_k)
                else
                    # perform the ionization: 
                    if (length(particles_ion) <= pia.index_last[species_ion])
                        resize!(particles_ion, length(particles_ion)+DELTA_PARTICLES)
                    end
                    # create the ion particle
                    update_buffer_index_new_particle!(particles_ion, pia, cell, species_ion)

                    p_i_new = particles_ion[pia.index_last[species_ion]]

                    p_i_new.w = particles_n_i.w
                    p_i_new.v = particles_n_i.v
                    p_i_new.x = particles_n_i.x

                    # add a second electron
                    if (length(particles_e) <= pia.index_last[species_e])
                        resize!(particles_e, length(particles_e)+DELTA_PARTICLES)
                    end
                    update_buffer_index_new_particle!(particles_e, pia, cell, species_e)

                    p_e_new = particles_e[pia.index_last[species_e]]

                    p_e_new.w = particles_n_i.w
                    p_e_new.v = particles_e_k.v
                    p_e_new.x = particles_e_k.x

                    # set neutral particle weight to 0
                    particles_n_i.w = 0.0

                    # compute energy split across the primare and secondary electrons
                    compute_g_new_ionization!(collision_data, interaction_l,
                                              get_ionization_threshold(n_e_interactions, species_n), get_electron_energy_split(n_e_interactions, species_n))

                    scatter_ionization_electrons_and_ion!(rng, collision_data, particles_e_k, p_e_new,
                                                          p_i_new, mass_ratio)
                end
            end
        end
    end
end

"""
    ntc_n_e_es!(rng, collision_factors, collision_data, interaction,
             n_e_interactions, n_e_cs, particles_n, particles_e, particles_ion,
             pia, cell, species_n, species_e, species_ion, Δt, V; extend::CSExtend=CSExtendConstant, dw_tol=1e-16)

Perform electron-neutral elastic scattering and electron-impact ionization collisions
using the event splitting method.

# Positional arguments
* `rng`: the random number generator
* `collision_factors`: the `CollisionFactors` for the species in question in the cell
* `collision_data`: `CollisionData` instance used for storing collisional quantities
* `interaction`: 2-dimensional array of `Interaction` instances for all possible species pairs
* `n_e_interactions`: the `ElectronNeutralInteractions` instance
* `n_e_cs`: a vector of `ComputedCrossSections` instances
* `particles_n`: `ParticleVector` of the particles of neutral species
* `particles_e`: `ParticleVector` of the particles of the electron species
* `particles_ion`: `ParticleVector` of the particles of the ion species
* `pia`: the `ParticleIndexerArray`
* `cell`: the index of the cell in which collisions are performed
* `species_n`: the index of the neutral species
* `species_e`: the index of the electron species
* `species_ion`: the index of the ion species
* `Δt`: timestep
* `V`: cell volume

# Keyword arguments
* `extend`: enum of `CSExtend` type that sets how out-of-range energy values are treated when computing cross-sections
* `dw_tol`: if weights of particles differ by less than this amount, an equal-weight collision is assumed
and no particle splitting is performed

# References
* G. Oblapenko, D. Goldstein, P. Varghese, C. Moore, Hedging direct simulation Monte Carlo bets via event splitting.
    [J. Comput. Phys, 2022](https://doi.org/10.1016/j.jcp.2022.111390).
* D.P. Schmidt, C.J. Rutland, A New Droplet Collision Algorithm.
    [J. Comput. Phys, 2000](https://doi.org/10.1006/jcph.2000.6568).
"""
function ntc_n_e_es!(rng, collision_factors, collision_data, interaction,
    n_e_interactions, n_e_cs, particles_n, particles_e, particles_ion, 
    pia, cell, species_n, species_e, species_ion, Δt, V; extend::CSExtend=CSExtendConstant, dw_tol=1e-16)
    # event splitting
    # compute ncoll
    # loop over particles
    # update sigma_g_w_max
    # collide
    @inbounds indexer_n = pia.indexer[cell, species_n]
    @inbounds indexer_e = pia.indexer[cell, species_e]

    @inbounds collision_factors.n1 = indexer_n.n_local
    @inbounds collision_factors.n2 = indexer_e.n_local
    @inbounds n_coll_float = compute_n_coll_two_species(rng, collision_factors, indexer_n.n_local, indexer_e.n_local, Δt, V)
    n_coll_int = floor(Int64, n_coll_float)
    # println(n_coll_float, ", ", n_coll_int)

    collision_factors.n_coll = n_coll_int
    collision_factors.n_coll_performed = 0
    collision_factors.n_eq_w_coll_performed = 0

    # need to index into n_e_cs
    @inbounds species_n_en_i = n_e_interactions.neutral_indexer[species_n]
    
    @inbounds mass_ratio = n_e_interactions.mass_ratios[species_n_en_i]

    @inbounds interaction_l = interaction[species_n, species_e]

    @inbounds for _ in 1:n_coll_int
        i = floor(Int64, rand(rng, Float64) * indexer_n.n_local)
        k = floor(Int64, rand(rng, Float64) * indexer_e.n_local)

        # example: bounds from [1,4], [7,9]; n_total = 7
        # n_group1 = 4, n_group2 = 3
        # i = 0,1,2,3 - [1,4]
        # i = 4,5,6 - [7,9]
        i = map_cont_index(indexer_n, i)
        k = map_cont_index(indexer_e, k)

        particles_n_i = particles_n[i]
        particles_e_k = particles_e[k]

        compute_g!(collision_data, particles_n_i, particles_e_k)

        if (collision_data.g > eps())

            collision_data.E_coll_eV = compute_cross_sections!(n_e_cs, interaction_l, collision_data.g, n_e_interactions, species_n;
                                                               extend=extend)
            sigma_g_w_max = get_cs_total(n_e_interactions, n_e_cs, species_n) * collision_data.g * max(particles_n_i.w, particles_e_k.w)

            # update (σ g w)_max if needed
            collision_factors.sigma_g_w_max = max(sigma_g_w_max, collision_factors.sigma_g_w_max)

            if (rand(rng, Float64) * collision_factors.sigma_g_w_max < sigma_g_w_max)
                collision_factors.n_coll_performed += 1
                compute_com!(collision_data, interaction_l, particles_n_i, particles_e_k)

                # do collision

                if (particles_n_i.w > particles_e_k.w)
                    # we split particle i, update velocity of i and k (split part remains unchanged)
                    if (length(particles_n) <= pia.index_last[species_n])
                        resize!(particles_n, length(particles_n)+DELTA_PARTICLES)
                    end
                    # first need to update the particle indexer struct
                    update_buffer_index_new_particle!(particles_n, pia, cell, species_n)

                    Δw = particles_n_i.w - particles_e_k.w
                    particles_n_i.w = particles_e_k.w

                    p_n_new = particles_n[pia.index_last[species_n]]

                    p_n_new.w = Δw
                    p_n_new.v = particles_n_i.v
                    p_n_new.x = particles_n_i.x
                elseif (abs(particles_n_i.w - particles_e_k.w) < dw_tol)
                    collision_factors.n_eq_w_coll_performed += 1
                else  # (particles[k].w > particles[i].w)
                    if (length(particles_e) <= pia.index_last[species_e])
                        resize!(particles_e, length(particles_e)+DELTA_PARTICLES)
                    end
                    update_buffer_index_new_particle!(particles_e, pia, cell, species_e)

                    Δw = particles_e_k.w - particles_n_i.w
                    particles_e_k.w = particles_n_i.w

                    p_e_new = particles_e[pia.index_last[species_e]]

                    p_e_new.w = Δw
                    p_e_new.v = particles_e_k.v
                    p_e_new.x = particles_e_k.x
                end

                # now we collide the 2 equal-weight particles
                
                if (n_e_cs[species_n_en_i].prob_vec[2] == 0.0)
                    scatter_vhs!(rng, collision_data, interaction_l, particles_n_i, particles_e_k)
                else
                    w_ionized = particles_e_k.w * n_e_cs[species_n_en_i].prob_vec[2]
    
                    particles_n_i.w -= w_ionized
                    particles_e_k.w -= w_ionized

                    if (length(particles_ion) <= pia.index_last[species_ion])
                        resize!(particles_ion, length(particles_ion)+DELTA_PARTICLES)
                    end
                    # create the ion particle
                    update_buffer_index_new_particle!(particles_ion, pia, cell, species_ion)

                    # velocity will be set later in scattering function
                    p_ion_new = particles_ion[pia.index_last[species_ion]]
                    p_ion_new.w = w_ionized
                    p_ion_new.x = particles_n_i.x

                    # add 2 electrons (split + secondary)
                    if (length(particles_e) < pia.index_last[species_e] + 2)
                        resize!(particles_e, length(particles_e)+DELTA_PARTICLES)
                    end
                    update_buffer_index_new_particle!(particles_e, pia, cell, species_e)

                    p_e_new1 = particles_e[pia.index_last[species_e]]
                    p_e_new1.w = w_ionized
                    p_e_new1.v = particles_e_k.v
                    p_e_new1.x = particles_e_k.x

                    update_buffer_index_new_particle!(particles_e, pia, cell, species_e)

                    p_e_new2 = particles_e[pia.index_last[species_e]]
                    p_e_new2.w = w_ionized
                    p_e_new2.v = particles_e_k.v
                    p_e_new2.x = particles_e_k.x

                    # elastic scattering
                    scatter_vhs!(rng, collision_data, interaction_l, particles_n_i, particles_e_k)

                    compute_g_new_ionization!(collision_data, interaction_l,
                                              get_ionization_threshold(n_e_interactions, species_n), get_electron_energy_split(n_e_interactions, species_n))

                    scatter_ionization_electrons_and_ion!(rng, collision_data, p_e_new1, p_e_new2,
                                                          p_ion_new, mass_ratio)
                end
            end
        end
    end
end

end