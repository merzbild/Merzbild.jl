@muladd begin

"""
    accelerate_constant_field_x!(particles, pia, cell, species, species_data, E, Δt)

Accelerate particles with a constant electric field in the X direction; no
sorting of particles is required since the field is constant.

# Positional arguments
* `particles`: vector-like structure of particles to be accelerated
* `pia`: ParticleIndexerArray instance
* `cell`: index of the cell in which particles are being accelerated
* `species`: index of the species of the particles being accelerated
* `species_data`: a `Vector{Species}` instance with the species' data
* `E`: value of the electric field in V/m
* `Δt`: timestep for which the acceleration is performed
"""
function accelerate_constant_field_x!(particles, pia, cell, species, species_data, E, Δt)
    # for 1-D problems, assume E-field is in the X direction
    dv = species_data[species].charge_div_mass * E * Δt

    @inbounds s1 = pia.indexer[cell,species].start1
    @inbounds e1 = pia.indexer[cell,species].end1

    @inbounds for i in s1:e1
        p = particles[i]
        p.v = SVector{3, Float64}(p.v[1] + dv,
                                  p.v[2],
                                  p.v[3])
    end

    @inbounds s2 = pia.indexer[cell,species].start2
    if s2 > 0
        @inbounds e2 = pia.indexer[cell,species].end2
        @inbounds for i in s2:e2
            p = particles[i]
            p.v = SVector{3, Float64}(p.v[1] + dv,
                                    p.v[2],
                                    p.v[3])
        end
    end
end

"""
    accelerate_electric_field_x!(grid::Grid1DUniform, particles, pia, cell, species, species_data, field_props, Δt)

Accelerate the particles of a species in a single cell of a 1-D uniform grid with the self-consistent
electric field stored in the nodes, using first-order (cloud-in-cell) interpolation of the field to
the position of a particle. The interpolation uses the same weights as [`deposit_charge!`](@ref).

The particles are assumed to be sorted on the grid (so no particles are indexed by group2).

# Positional arguments
* `grid`: the `Grid1DUniform` grid
* `particles`: the `ParticleVector` of the particles being accelerated
* `pia`: the `ParticleIndexerArray` instance
* `cell`: index of the cell in which particles are being accelerated
* `species`: index of the species of the particles being accelerated
* `species_data`: a `Vector{Species}` instance with the species' data
* `field_props`: the `ElectrostaticFieldProps` instance holding the electric field
* `Δt`: timestep for which the acceleration is performed
"""
function accelerate_electric_field_x!(grid::Grid1DUniform, particles, pia, cell, species, species_data,
                                      field_props, Δt)
    # for 1-D problems, assume E-field is in the X direction
    @inbounds dv_factor = species_data[species].charge_div_mass * Δt

    inv_Δx = grid.inv_Δx
    shift = cell - 1

    @inbounds E_left = field_props.electric_field[cell]
    @inbounds ΔE = field_props.electric_field[cell+1] - E_left

    @inbounds s1 = pia.indexer[cell,species].start1
    @inbounds e1 = pia.indexer[cell,species].end1

    @inbounds for i in s1:e1
        p = particles[i]
        ξ = p.x[1] * inv_Δx - shift
        dv = dv_factor * (E_left + ξ * ΔE)
        p.v = SVector{3, Float64}(p.v[1] + dv,
                                  p.v[2],
                                  p.v[3])
    end
end

"""
    accelerate_electric_field_x!(grid::Grid1DUniform, particles, pia, species, species_data, field_props, Δt)

Accelerate the particles of a species in all cells of a 1-D uniform grid with the self-consistent
electric field stored in the nodes.

The particles are assumed to be sorted on the grid.

# Positional arguments
* `grid`: the `Grid1DUniform` grid
* `particles`: the `ParticleVector` of the particles being accelerated
* `pia`: the `ParticleIndexerArray` instance
* `species`: index of the species of the particles being accelerated
* `species_data`: a `Vector{Species}` instance with the species' data
* `field_props`: the `ElectrostaticFieldProps` instance holding the electric field
* `Δt`: timestep for which the acceleration is performed
"""
function accelerate_electric_field_x!(grid::Grid1DUniform, particles, pia, species, species_data,
                                      field_props, Δt)
    for cell in 1:grid.n_cells
        accelerate_electric_field_x!(grid, particles, pia, cell, species, species_data, field_props, Δt)
    end
end

end