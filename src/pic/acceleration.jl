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

    for i in pia.indexer[cell,species].start1:pia.indexer[cell,species].end1
        particles[i].v = SVector{3, Float64}(particles[i].v[1] + dv,
                                             particles[i].v[2],
                                             particles[i].v[3])
    end

    if pia.indexer[cell,species].start2 > 0
        for i in pia.indexer[cell,species].start2:pia.indexer[cell,species].end2
            particles[i].v = SVector{3, Float64}(particles[i].v[1] + dv,
                                                 particles[i].v[2],
                                                 particles[i].v[3])
        end
    end
end

end