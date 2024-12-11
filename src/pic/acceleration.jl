"""
Accelerate particles with a constant electric field in the X direction
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