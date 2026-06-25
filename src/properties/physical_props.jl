using LinearAlgebra

@muladd begin

"""
    PhysProps

Structure to store computed physical properties in a physical cell.

# Fields
* `ndens_not_Np`: whether the `n` field stores number density (if `true`) and not the number of physical particles in a cell (if `false`)
* `n_cells`: number of physical cells
* `n_species`: number of species
* `lpa`: length of the particle array (vector of length `n_species`)
* `np`: number of particles (array of shape `(n_cells, n_species)`)
* `n`: number density or number of physical particles in a cell (array of shape `(n_cells, n_species)`)
* `v`: per-species flow velocity in a cell (array of shape `(3, n_cells, n_species)`)
* `T`: per-species temperature in a cell (array of shape `(n_cells, n_species)`)
* `Tref`: reference temperature used to scale moments
"""
mutable struct PhysProps
    ndens_not_Np::Bool
    n_cells::Int64
    n_species::Int64
    lpa::Vector{Float64}  # length of particle array: species
    np::Array{Float64,2}  # number of particles: cells x species
    n::Array{Float64,2}  # number density: cells x species
    v::Array{Float64,3}  # velocity: velocity component x cells x species
    T::Array{Float64,2}  # temperature: cells x species
end

"""
    PhysProps(n_cells, n_species; ndens_not_Np=false)

Construct physical properties given the number of cells and species.

# Positional arguments
* `n_cells`: number of cells
* `n_species`: number of species

# Keyword arguments
* `ndens_not_Np`: whether the `n` field stores number density (if `true`) and not the number of physical particles in a cell (if `false`)

"""
PhysProps(n_cells, n_species; ndens_not_Np=false, Tref=300.0) = PhysProps(ndens_not_Np, n_cells, n_species,
                                                                    zeros(n_species), zeros(n_cells, n_species),
                                                                    zeros(n_cells, n_species), zeros(3, n_cells, n_species), zeros(n_cells, n_species))

"""
    PhysProps(pia::ParticleIndexerArray; ndens_not_Np=false)

Construct physical properties given a `ParticleIndexerArray` instance,.

# Positional arguments
* `pia`: the `ParticleIndexerArray` instance

# Keyword arguments
* `ndens_not_Np`: whether the `n` field stores number density (if `true`) and not the number of physical particles in a cell (if `false`)
"""
PhysProps(pia::ParticleIndexerArray; ndens_not_Np=false) = PhysProps(size(pia.indexer)[1], size(pia.indexer)[2], ndens_not_Np=ndens_not_Np)

"""
    compute_props!(particles, pia, species_data, phys_props)

Compute the physical properties of all species in all cells and store the result in a `PhysProps` instance.
This function does not compute the total moments, even if `phys_props.n_moments > 0`.

# Positional arguments
* `particles`: the `Vector` of `ParticleVector`s containing all the particles in a simulation
* `pia`: the `ParticleIndexerArray` instance
* `species_data`: the `Vector` of `SpeciesData`
* `phys_props`: the `PhysProps` instance in which the computed physical properties are stored
"""
function compute_props!(particles, pia, species_data, phys_props)
    @inbounds for species in 1:phys_props.n_species
        for cell in 1:phys_props.n_cells
            np = 0
            n = 0.0
            E = 0.0
            T = 0.0
            v = SVector{3,Float64}(0.0, 0.0, 0.0)

            for i in pia.indexer[cell,species].start1:pia.indexer[cell,species].end1
                n += particles[species][i].w
                v = v + particles[species][i].v * particles[species][i].w
                np += 1
            end

            if pia.indexer[cell,species].n_group2 > 0
                for i in pia.indexer[cell,species].start2:pia.indexer[cell,species].end2
                    n += particles[species][i].w
                    v = v + particles[species][i].v * particles[species][i].w
                    np += 1
                end
            end

            if (n > 0.0)
                v /= n
                for i in pia.indexer[cell,species].start1:pia.indexer[cell,species].end1
                    E = E + particles[species][i].w * ((particles[species][i].v[1] - v[1])^2
                                                  + (particles[species][i].v[2] - v[2])^2
                                                  + (particles[species][i].v[3] - v[3])^2)
                end
            
                if pia.indexer[cell,species].n_group2 > 0
                    for i in pia.indexer[cell,species].start2:pia.indexer[cell,species].end2
                        E = E + particles[species][i].w * ((particles[species][i].v[1] - v[1])^2
                                                      + (particles[species][i].v[2] - v[2])^2
                                                      + (particles[species][i].v[3] - v[3])^2)
                    end
                end
                
                E *= 0.5 * species_data[species].mass / (n * k_B)
                T = (2.0/3.0) * E
            end

            phys_props.lpa[species] = length(particles[species])
            phys_props.np[cell,species] = np
            phys_props.n[cell,species] = n
            phys_props.v[:,cell,species] = v
            phys_props.T[cell,species] = T
        end
    end
end

"""
    clear_props!(phys_props::PhysProps)

Clear all data from PhysProps, for use when physical properties are averaged over timesteps
and averaging over a new set of timesteps needs to be started.

# Positional arguments
* `phys_props`: the `PhysProps` instance to be cleared
"""
function clear_props!(phys_props::PhysProps)
    fill!(phys_props.lpa, 0)
    fill!(phys_props.np, 0)
    fill!(phys_props.n, 0.0)
    fill!(phys_props.v, 0.0)
    fill!(phys_props.T, 0.0)
end

"""
    avg_props!(phys_props_avg::PhysProps, phys_props::PhysProps, n_avg_timesteps)

Used to time-average computed physical properties, not including the total moments.
For each instantaneous value of a property computed and stored in `phys_props`,
it is divided by `n_avg_timesteps` and added to `phys_props_avg`.

# Positional arguments
* `phys_props_avg`: the `PhysProps` instance used to store the time-averaged properties
* `phys_props`: the `PhysProps` instance holding the current values of the properties
    to be used for the averaging at the current timestep
* `n_avg_timesteps`: the number of timesteps over which the averaging is performed

# Throws
`ErrorException` if `phys_props_avg` computes number density and `phys_props` computes number of physical
particles, or vice versa.
"""
function avg_props!(phys_props_avg::PhysProps, phys_props::PhysProps, n_avg_timesteps)
    if (phys_props_avg.ndens_not_Np != phys_props.ndens_not_Np)
        throw(ErrorException("Inconsistent computation of ndens/number of physical particles in cell"))
    end

    inv_nt_avg = 1.0 / n_avg_timesteps
    n_cells = phys_props.n_cells
    n_species = phys_props.n_species

    avg_lpa = phys_props_avg.lpa
    curr_lpa = phys_props.lpa
    
    avg_np = phys_props_avg.np
    curr_np = phys_props.np
    
    avg_n = phys_props_avg.n
    curr_n = phys_props.n
    
    avg_v = phys_props_avg.v
    curr_v = phys_props.v
    
    avg_T = phys_props_avg.T
    curr_T = phys_props.T

    @inbounds for species in 1:n_species
        avg_lpa[species] += curr_lpa[species] * inv_nt_avg

        @simd for cell in 1:n_cells
            avg_np[cell, species] += curr_np[cell, species] * inv_nt_avg
            avg_n[cell, species]  += curr_n[cell, species] * inv_nt_avg
            
            avg_v[1, cell, species] += curr_v[1, cell, species] * inv_nt_avg
            avg_v[2, cell, species] += curr_v[2, cell, species] * inv_nt_avg
            avg_v[3, cell, species] += curr_v[3, cell, species] * inv_nt_avg
            
            avg_T[cell, species]  += curr_T[cell, species] * inv_nt_avg
        end
    end
end

"""
    compute_props_sorted!(particles, pia, species_data, phys_props, cell_chunk)

Compute the physical properties of all species in a
subset of cells and store the result in a `PhysProps` instance,
assuming the particles are sorted.
This function does not compute the total moments, even if `phys_props.n_moments > 0`. Currently this does not
compute the length of the particle array.

# Positional arguments
* `particles`: the `Vector` of `ParticleVector`s containing all the particles in a simulation
* `pia`: the `ParticleIndexerArray` instance
* `species_data`: the `Vector` of `SpeciesData`
* `phys_props`: the `PhysProps` instance in which the computed physical properties are stored
* `cell_chunk`: the list of cell indices or range of cell indices in which to compute the properties
"""
function compute_props_sorted!(particles, pia, species_data, phys_props, cell_chunk)
    indexer = pia.indexer
    pp_np = phys_props.np
    pp_n  = phys_props.n
    pp_v  = phys_props.v
    pp_T  = phys_props.T

    @inbounds for species in 1:phys_props.n_species
        species_particles = particles[species]
        species_mass = species_data[species].mass

        for cell in cell_chunk
            n = 0.0
            np = 0.0
            E = 0.0
            T = 0.0
            v = SVector{3,Float64}(0.0, 0.0, 0.0)

            idx = indexer[cell, species]
            s1 = idx.start1
            e1 = idx.end1
            for i in s1:e1
                particle = species_particles[i]

                n += particle.w
                v = v + particle.v * particle.w
            end

            np = e1 >= s1 ? e1-s1 + 1.0 : 0.0

            if (n > 0.0)
                v /= n
                for i in s1:e1
                    particle = species_particles[i]
                    
                    E = E + particle.w * ((particle.v[1] - v[1])^2
                                            + (particle.v[2] - v[2])^2
                                            + (particle.v[3] - v[3])^2)
                end
                E *= 0.5 * species_mass / (n * k_B)
                T = (2.0/3.0) * E
            end
            
            pp_np[cell, species] = np
            pp_n[cell, species] = n
            pp_v[1, cell, species] = v[1]
            pp_v[2, cell, species] = v[2]
            pp_v[3, cell, species] = v[3]
            pp_T[cell, species] = T
        end
    end
end

"""
    compute_props_sorted!(particles, pia, species_data, phys_props)

Compute the physical properties of all species in all cells and store the result in a `PhysProps` instance,
assuming the particles are sorted.
This function does not compute the total moments, even if `phys_props.n_moments > 0`. Currently this does not
compute the length of the particle array.

# Positional arguments
* `particles`: the `Vector` of `ParticleVector`s containing all the particles in a simulation
* `pia`: the `ParticleIndexerArray` instance
* `species_data`: the `Vector` of `SpeciesData`
* `phys_props`: the `PhysProps` instance in which the computed physical properties are stored
"""
@inline function compute_props_sorted!(particles, pia, species_data, phys_props)
    compute_props_sorted!(particles, pia, species_data, phys_props, 1:phys_props.n_cells)
end

"""
    compute_props_sorted!(particles::Vector{ParticleVector{D}}, pia, species_data, phys_props, grid::G, cell_chunk) where {G<:AbstractGrid,D}

Compute the physical properties of all species in a
subset of cells and store the result in a `PhysProps` instance,
assuming the particles are sorted.
This function does not compute the total moments, even if `phys_props.n_moments > 0`.
If `ndens_not_Np` is `true`, the number density will be computed based on the volumes of the grid cells;
otherwise, the number of physical particles in each cell will be computed. Currently this does not
compute the length of the particle array.

# Positional arguments
* `particles`: the `Vector` of `ParticleVector`s containing all the particles in a simulation
* `pia`: the `ParticleIndexerArray` instance
* `species_data`: the `Vector` of `SpeciesData`
* `phys_props`: the `PhysProps` instance in which the computed physical properties are stored
* `grid`: the physical grid
* `cell_chunk`: the list of cell indices or range of cell indices in which to compute the properties
"""
function compute_props_sorted!(particles::Vector{ParticleVector{D}}, pia, species_data, phys_props, grid::G, cell_chunk) where {G<:AbstractGrid,D}
    if !phys_props.ndens_not_Np
        compute_props_sorted!(particles, pia, species_data, phys_props)
    else
        @inbounds for species in 1:phys_props.n_species
            for cell in cell_chunk
                n = 0.0
                np = 0.0
                E = 0.0
                T = 0.0
                v = SVector{3,Float64}(0.0, 0.0, 0.0)

                s1 = pia.indexer[cell,species].start1
                e1 = pia.indexer[cell,species].end1
                for i in s1:e1
                    particle = particles[species][i]

                    n += particle.w
                    v = v + particle.v * particle.w
                    np += 1.0
                end

                if (n > 0.0)
                    v /= n
                    for i in s1:e1
                        particle = particles[species][i]

                        E = E + particle.w * ((particle.v[1] - v[1])^2
                                              + (particle.v[2] - v[2])^2
                                              + (particle.v[3] - v[3])^2)
                    end
                    E *= 0.5 * species_data[species].mass / (n * k_B)
                    T = (2.0/3.0) * E
                end
        
                phys_props.np[cell,species] = np
                phys_props.n[cell,species] = n * grid.cells[cell].inv_V
                phys_props.v[1,cell,species] = v[1]
                phys_props.v[2,cell,species] = v[2]
                phys_props.v[3,cell,species] = v[3]
                phys_props.T[cell,species] = T
            end
        end
    end
end

"""
    compute_props_sorted!(particles::Vector{ParticleVector{D}}, pia, species_data, phys_props, grid::AbstractGrid) where {G<:AbstractGrid,D}

Compute the physical properties of all species in all cells and store the result in a `PhysProps` instance,
assuming the particles are sorted.
This function does not compute the total moments, even if `phys_props.n_moments > 0`.
If `ndens_not_Np` is `true`, the number density will be computed based on the volumes of the grid cells;
otherwise, the number of physical particles in each cell will be computed. Currently this does not
compute the length of the particle array.

# Positional arguments
* `particles`: the `Vector` of `ParticleVector`s containing all the particles in a simulation
* `pia`: the `ParticleIndexerArray` instance
* `species_data`: the `Vector` of `SpeciesData`
* `phys_props`: the `PhysProps` instance in which the computed physical properties are stored
* `grid`: the physical grid
"""
@inline function compute_props_sorted!(particles::Vector{ParticleVector{D}}, pia, species_data, phys_props, grid::G) where {G<:AbstractGrid,D}
    compute_props_sorted!(particles, pia, species_data, phys_props, grid, 1:phys_props.n_cells)
end

"""
    compute_mixed_moment(particles, pia, cell, species, powers; sum_scaler=1.0, res_scaler=1.0)

Compute mixed velocity moment of particles in a cell: ``\\sum_i w_i v_{x,i}^{p_x} v_{y,i}^{p_y} v_{z,i}^{p_z}``.

# Positional arguments
* `particles`: the `Vector` of `ParticleVector`s containing all the particles in a simulation
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the cell in which the moment is being computed
* `species`: the species for which the mixed moment is being computed
* `powers`: a `Vector` of the powers to which the x-, y-, and z-components of the velocity are to be raised

# Keyword arguments
* `sum_scaler`: during the summation over the particles, scale each summand by this factor
    potentially reduce round-off issues
* `res_scaler`: scaling factor by which to multiply the result at the end
"""
function compute_mixed_moment(particles, pia, cell, species, powers; sum_scaler=1.0, res_scaler=1.0)
    # sum scaler is used inside the particle summation loops to potentially reduce round-off issues
    # res_scaler can be used as inverse of sum_scaler (e.g. to get the full moment)
    # or set to any other quantity to get scaling of result as well
    # e.g. sum_scaler=(1.0/ndens) (computed separately), res_scaler=1.0 would compute normalized moment
    result = 0.0


    s1 = pia.indexer[cell,species].start1
    e1 = pia.indexer[cell,species].end1

    @inbounds for i in s1:e1
        result += particles[species][i].w * sum_scaler * (particles[species][i].v[1]^powers[1]) *
                  (particles[species][i].v[2]^powers[2]) *
                  (particles[species][i].v[3]^powers[3])
    end

    if pia.indexer[cell,species].n_group2 > 0
        s2 = pia.indexer[cell,species].start2
        e2 = pia.indexer[cell,species].end2
        @inbounds for i in s2:e2
            result += particles[species][i].w * sum_scaler * (particles[species][i].v[1]^powers[1]) *
                        (particles[species][i].v[2]^powers[2]) *
                        (particles[species][i].v[3]^powers[3])
        end
    end

    return result * res_scaler
end

end