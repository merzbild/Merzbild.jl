@muladd begin

"""
    accelerate_constant_field_x!(particles, pia, cell, species, species_data, E, Δt)

Accelerate particles with a constant electric field in the X direction

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

end