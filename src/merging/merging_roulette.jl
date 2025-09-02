"""
    merge_roulette!(rng, particles, pia, cell, species, target_np)

Perform roulette merging - delete random particles until target number of particles is reached, and re-weight remaining particles
to conserve number density.

# Positional arguments
* `rng`: the random number generator instance
* `particles`: the `ParticleVector` instance containing the particles to be merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the grid cell in which particles are being merged
* `species`: the index of the species being merged
* `target_np`: the target post-merge number of particles

# References
* J. Watrous, D.B. Seidel, C.H. Moore, W. McDoniel,
    Improvements to Particle Merge Algorithms for Sandia National Laboratories Plasma Physics Modeling
    Code, EMPIRE. [Presentation, 2023](https://www.osti.gov/servlets/purl/2431184).
"""
function merge_roulette!(rng, particles, pia, cell, species, target_np)
    current_count = pia.indexer[cell, species].n_local

    n_to_delete = current_count - target_np

    w_deleted = 0.0

    for _ in 1:n_to_delete
        i_delete = floor(Int64, rand(rng, Float64) * pia.indexer[cell, species].n_local)
        i_delete = map_cont_index(pia.indexer[cell, species], i_delete)

        w_deleted += particles[i_delete].w
        delete_particle!(pv, pia, cell, species, i_delete)
    end

    w_total = 0.0
    s1 = pia.indexer[cell,species].start1
    e1 = pia.indexer[cell,species].end1
    for i in s1:e1
        w_total += particles[i].w
    end

    if pia.indexer[cell, species].n_group2 > 0
        s2 = pia.indexer[cell,species].start2
        e2 = pia.indexer[cell,species].end2
        for i in s2:e2
            w_total += particles[i].w
        end
    end

    scale_factor = (w_total + w_deleted) / w_total

    s1 = pia.indexer[cell,species].start1
    e1 = pia.indexer[cell,species].end1
    for i in s1:e1
        particles[i].w *= scale_factor
    end

    if pia.indexer[cell, species].n_group2 > 0
        s2 = pia.indexer[cell,species].start2
        e2 = pia.indexer[cell,species].end2
        for i in s2:e2
            particles[i].w *= scale_factor
        end
    end
end