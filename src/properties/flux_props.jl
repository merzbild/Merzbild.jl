using LinearAlgebra

@muladd begin

"""
    FluxProps

Structure to store computed flux densities in a physical cell.
The kinetic energy flux density of a species is computed as
``\\frac{m}{2V} \\sum_i w_i \\mathbf{C}_i C_i^2``, where ``\\mathbf{C}_i`` is the peculiar velocity
of a particle, ``m`` is the molecular mass of the species, ``V`` is the cell volume.
The ``k,l`` component of the the momentum flux density tensor is computed as
``\\frac{m}{2V} \\sum_i w_i \\C_{i,k} \\C_{i,l}``.
The ``xx``, ``yy``, ``zz`` components thereof are stored in `diagonal_momentum_flux`.
The ``xy``, ``xz``, ``yz`` components thereof are stored in `off_diagonal_momentum_flux`.


# Fields
* `n_cells`: number of physical cells
* `n_species`: number of species
* `kinetic_energy_flux`: per-species kinetic energy flux in a cell (array of shape `(3, n_cells, n_species)`)
* `diagonal_momentum_flux`: diagonal components of the per-species momentum flux tensor
    in a cell (array of shape `(3, n_cells, n_species)`)
* `off_diagonal_momentum_flux`: off-diagonal components of the per-species momentum flux tensor
    in a cell (array of shape `(3, n_cells, n_species)`)
"""
mutable struct FluxProps
    n_cells::Int64
    n_species::Int64
    kinetic_energy_flux::Array{Float64,3}  # kinetic energy flux: component x cells x species
    diagonal_momentum_flux::Array{Float64,3}  # diagonal components of the momentum flux density tensor: component x cells x species
    off_diagonal_momentum_flux::Array{Float64,3}  # off-diagonal components of the momentum flux density tensor: component x cells x species
end

"""
    FluxProps(n_cells, n_species)

Construct a `FluxProps` instance given the number of cells and species.

# Positional arguments
* `n_cells`: number of cells
* `n_species`: number of species
"""
FluxProps(n_cells, n_species) = FluxProps(n_cells, n_species,
                                          zeros(3, n_cells, n_species),
                                          zeros(3, n_cells, n_species),
                                          zeros(3, n_cells, n_species))

"""
    FluxProps(pia)

Construct a `FluxProps` instance given a `ParticleIndexerArray` instance.

# Positional arguments
* `pia`: the `ParticleIndexerArray` instance
"""
FluxProps(pia) = FluxProps(size(pia.indexer)[1], size(pia.indexer)[2])

"""
    compute_flux_props!(particles, pia, species_data, phys_props::PhysProps, flux_props::FluxProps, grid::G) where {G<:AbstractGrid}

Compute the fluxes of all species in all cells and store the result in a `FluxProps` instance.
This uses the pre-computed species-wise mean velocities from a `PhysProps` instance, which needs to be computed
at the same timestep before calling this function.
Particles of a cell can be distributed across `group1` and `group2` of indices pointed to by the `ParticleIndexerArray`
instance, i.e. this function can be used to compute fluxes immediately after collisions (but before convection).
In case particles are sorted, `compute_flux_props_sorted!` will be more efficient.

# Positional arguments
* `particles`: the `Vector` of `ParticleVector`s containing all the particles in a simulation
* `pia`: the `ParticleIndexerArray` instance
* `species_data`: the `Vector` of `SpeciesData`
* `phys_props`: the `PhysProps` instance computed at the same timestep
* `flux_props`: the `FluxProps` instance in which the computed fluxes are stored
* `grid`: the physical grid
"""
function compute_flux_props!(particles, pia, species_data, phys_props::PhysProps, flux_props::FluxProps, grid::G) where {G<:AbstractGrid}
    for species in 1:flux_props.n_species
        @inbounds for cell in 1:flux_props.n_cells

            kefd = SVector{3,Float64}(0.0, 0.0, 0.0)  # kinetic_energy_flux
            dmfd = SVector{3,Float64}(0.0, 0.0, 0.0)  # diagonal_momentum_flux
            odmfd = SVector{3,Float64}(0.0, 0.0, 0.0)  # off_diagonal_momentum_flux
            c = SVector{3,Float64}(0.0, 0.0, 0.0)

            for i in pia.indexer[cell,species].start1:pia.indexer[cell,species].end1
                c = particles[species][i].v - SVector{3,Float64}(phys_props.v[1, cell, species],
                                                                 phys_props.v[2, cell, species],
                                                                 phys_props.v[3, cell, species])
                cxsq = c[1]^2
                cysq = c[2]^2
                czsq = c[3]^2
                
                csq = cxsq + cysq + czsq
                kefd = kefd + particles[species][i].w * c * csq
                dmfd = dmfd + particles[species][i].w * SVector{3,Float64}(cxsq, cysq, czsq)
                odmfd = odmfd + particles[species][i].w * SVector{3,Float64}(c[1]*c[2], c[1]*c[3], c[2]*c[3])
            end
        
            if pia.indexer[cell,species].n_group2 > 0
                for i in pia.indexer[cell,species].start2:pia.indexer[cell,species].end2
                    c = particles[species][i].v - SVector{3,Float64}(phys_props.v[1, cell, species],
                                                                    phys_props.v[2, cell, species],
                                                                    phys_props.v[3, cell, species])
                    cxsq = c[1]^2
                    cysq = c[2]^2
                    czsq = c[3]^2
                    csq = cxsq + cysq + czsq

                    kefd = kefd + particles[species][i].w * c * csq
                    dmfd = dmfd + particles[species][i].w * SVector{3,Float64}(cxsq, cysq, czsq)
                    odmfd = odmfd + particles[species][i].w * SVector{3,Float64}(c[1]*c[2], c[1]*c[3], c[2]*c[3])
                end
            end

            flux_props.kinetic_energy_flux[:,cell,species] = 0.5 * kefd * species_data[species].mass * grid.cells[cell].inv_V
            flux_props.diagonal_momentum_flux[:,cell,species] = dmfd * species_data[species].mass * grid.cells[cell].inv_V
            flux_props.off_diagonal_momentum_flux[:,cell,species] = odmfd * species_data[species].mass * grid.cells[cell].inv_V
        end
    end
end

"""
    clear_props!(flux_props::FluxProps)

Clear all data from `FluxProps`, for use when flux densities are averaged over timesteps
and averaging over a new set of timesteps needs to be started.

# Positional arguments
* `flux_props`: the `FluxProps` instance to be cleared
"""
function clear_props!(flux_props::FluxProps)
    flux_props.kinetic_energy_flux[:,:,:] .= 0.0
    flux_props.diagonal_momentum_flux[:,:,:] .= 0.0
    flux_props.off_diagonal_momentum_flux[:,:,:] .= 0.0
end

"""
    avg_props!(flux_props_avg::FluxProps, flux_props::FluxProps, n_avg_timesteps)

Used to time-average computed flux densities.
For each instantaneous value of a property computed and stored in `flux_props`,
it is divided by `n_avg_timesteps` and added to `flux_props_avg`.

# Positional arguments
* `flux_props_avg`: the `FluxProps` instance used to store the time-averaged properties
* `flux_props`: the `FluxProps` instance holding the current values of the properties
    to be used for the averaging at the current timestep
* `n_avg_timesteps`: the number of timesteps over which the averaging is performed
"""
function avg_props!(flux_props_avg::FluxProps, flux_props::FluxProps, n_avg_timesteps)
    inv_nt_avg = 1.0 / n_avg_timesteps

    for species in 1:flux_props.n_species
        @inbounds @simd for cell in 1:flux_props.n_cells
            flux_props_avg.kinetic_energy_flux[1,cell,species] = flux_props_avg.kinetic_energy_flux[1,cell,species] + 
                                                                 flux_props.kinetic_energy_flux[1,cell,species] * inv_nt_avg
            flux_props_avg.kinetic_energy_flux[2,cell,species] = flux_props_avg.kinetic_energy_flux[2,cell,species] + 
                                                                 flux_props.kinetic_energy_flux[2,cell,species] * inv_nt_avg
            flux_props_avg.kinetic_energy_flux[3,cell,species] = flux_props_avg.kinetic_energy_flux[3,cell,species] + 
                                                                 flux_props.kinetic_energy_flux[3,cell,species] * inv_nt_avg

            flux_props_avg.diagonal_momentum_flux[1,cell,species] = flux_props_avg.diagonal_momentum_flux[1,cell,species] + 
                                                                    flux_props.diagonal_momentum_flux[1,cell,species] * inv_nt_avg
            flux_props_avg.diagonal_momentum_flux[2,cell,species] = flux_props_avg.diagonal_momentum_flux[2,cell,species] + 
                                                                    flux_props.diagonal_momentum_flux[2,cell,species] * inv_nt_avg
            flux_props_avg.diagonal_momentum_flux[3,cell,species] = flux_props_avg.diagonal_momentum_flux[3,cell,species] + 
                                                                    flux_props.diagonal_momentum_flux[3,cell,species] * inv_nt_avg

            flux_props_avg.off_diagonal_momentum_flux[1,cell,species] = flux_props_avg.off_diagonal_momentum_flux[1,cell,species] + 
                                                                        flux_props.off_diagonal_momentum_flux[1,cell,species] * inv_nt_avg
            flux_props_avg.off_diagonal_momentum_flux[2,cell,species] = flux_props_avg.off_diagonal_momentum_flux[2,cell,species] + 
                                                                        flux_props.off_diagonal_momentum_flux[2,cell,species] * inv_nt_avg
            flux_props_avg.off_diagonal_momentum_flux[3,cell,species] = flux_props_avg.off_diagonal_momentum_flux[3,cell,species] + 
                                                                        flux_props.off_diagonal_momentum_flux[3,cell,species] * inv_nt_avg
        end
    end
end

"""
    compute_flux_props_sorted!(particles, pia, species_data, phys_props, flux_props, grid::G, cell_chunk) where {G<:AbstractGrid}

Compute the flux densities of all species in a
subset of cells and store the result in a `FluxProps` instance,
assuming the particles are sorted.
This uses the pre-computed species-wise mean velocities from a `PhysProps` instance, which needs to be computed
at the same timestep before calling this function for the same subset of cells.

# Positional arguments
* `particles`: the `Vector` of `ParticleVector`s containing all the particles in a simulation
* `pia`: the `ParticleIndexerArray` instance
* `species_data`: the `Vector` of `SpeciesData`
* `phys_props`: the `PhysProps` instance computed at the same timestep for the same subset of cells
* `flux_props`: the `FluxProps` instance in which the computed fluxes are stored
* `grid`: the physical grid
* `cell_chunk`: the list of cell indices or range of cell indices in which to compute the properties
"""
function compute_flux_props_sorted!(particles, pia, species_data, phys_props, flux_props, grid::G, cell_chunk) where {G<:AbstractGrid}
    for species in 1:phys_props.n_species
        for cell in cell_chunk

            kefd = SVector{3,Float64}(0.0, 0.0, 0.0)  # kinetic_energy_flux
            dmfd = SVector{3,Float64}(0.0, 0.0, 0.0)  # diagonal_momentum_flux
            odmfd = SVector{3,Float64}(0.0, 0.0, 0.0)  # off_diagonal_momentum_flux
            c = SVector{3,Float64}(0.0, 0.0, 0.0)

            for i in pia.indexer[cell,species].start1:pia.indexer[cell,species].end1
                c = particles[species][i].v - SVector{3,Float64}(phys_props.v[1, cell, species],
                                                                 phys_props.v[2, cell, species],
                                                                 phys_props.v[3, cell, species])
                cxsq = c[1]^2
                cysq = c[2]^2
                czsq = c[3]^2
                
                csq = cxsq + cysq + czsq
                kefd = kefd + particles[species][i].w * c * csq
                dmfd = dmfd + particles[species][i].w * SVector{3,Float64}(cxsq, cysq, czsq)
                odmfd = odmfd + particles[species][i].w * SVector{3,Float64}(c[1]*c[2], c[1]*c[3], c[2]*c[3])
            end

            flux_props.kinetic_energy_flux[:,cell,species] = 0.5 * kefd * species_data[species].mass * grid.cells[cell].inv_V
            flux_props.diagonal_momentum_flux[:,cell,species] = dmfd * species_data[species].mass * grid.cells[cell].inv_V
            flux_props.off_diagonal_momentum_flux[:,cell,species] = odmfd * species_data[species].mass * grid.cells[cell].inv_V
        end
    end
end

"""
    compute_flux_props_sorted!(particles, pia, species_data, phys_props, flux_props, grid::G) where {G<:AbstractGrid}

Compute the flux densities of all species in all
cells and store the result in a `FluxProps` instance,
assuming the particles are sorted.
This uses the pre-computed species-wise mean velocities from a `PhysProps` instance, which needs to be computed
at the same timestep before calling this function.

# Positional arguments
* `particles`: the `Vector` of `ParticleVector`s containing all the particles in a simulation
* `pia`: the `ParticleIndexerArray` instance
* `species_data`: the `Vector` of `SpeciesData`
* `phys_props`: the `PhysProps` instance computed at the same timestep
* `flux_props`: the `FluxProps` instance in which the computed fluxes are stored
* `grid`: the physical grid
"""
@inline function compute_flux_props_sorted!(particles, pia, species_data, phys_props, flux_props, grid::G) where {G<:AbstractGrid}
    compute_flux_props_sorted!(particles, pia, species_data, phys_props, flux_props, grid, 1:phys_props.n_cells)
end

end