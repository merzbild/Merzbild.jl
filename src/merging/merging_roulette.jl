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

# Keyword arguments
* `conservative`: if `true`, the post-merge particles' velocities will be corrected to ensure conservation of
    momentum and energy

# References
* J. Watrous, D.B. Seidel, C.H. Moore, W. McDoniel,
    Improvements to Particle Merge Algorithms for Sandia National Laboratories Plasma Physics Modeling
    Code, EMPIRE. [Presentation, 2023](https://www.osti.gov/servlets/purl/2431184).
"""
function merge_roulette!(rng, particles, pia, cell, species, target_np; conservative=false)
    current_count = pia.indexer[cell, species].n_local

    n_to_delete = current_count - target_np

    # initial density, velocity, energy
    w_total0 = 0.0
    v0 = SVector{3,Float64}([0.0, 0.0, 0.0])
    E0 = 0.0
    if conservative
        w_total0 = 0.0

        @inbounds s1 = pia.indexer[cell,species].start1
        @inbounds e1 = pia.indexer[cell,species].end1
        @inbounds for i in s1:e1
            w_total0 += particles[i].w
            v0 += particles[i].w * particles[i].v
        end

        @inbounds if pia.indexer[cell, species].n_group2 > 0
            @inbounds s2 = pia.indexer[cell,species].start2
            @inbounds e2 = pia.indexer[cell,species].end2
            @inbounds for i in s2:e2
                w_total0 += particles[i].w
                v0 += particles[i].w * particles[i].v
            end
        end

        v0 /= w_total0

        # energy compute 
        @inbounds s1 = pia.indexer[cell,species].start1
        @inbounds e1 = pia.indexer[cell,species].end1
        @inbounds for i in s1:e1
            E0 += particles[i].w * ((particles[i].v[1] - v0[1])^2
                                    + (particles[i].v[2] - v0[2])^2
                                    + (particles[i].v[3] - v0[3])^2)
        end

        @inbounds if pia.indexer[cell, species].n_group2 > 0
            @inbounds s2 = pia.indexer[cell,species].start2
            @inbounds e2 = pia.indexer[cell,species].end2
            @inbounds for i in s2:e2
                E0 += particles[i].w * ((particles[i].v[1] - v0[1])^2
                                        + (particles[i].v[2] - v0[2])^2
                                        + (particles[i].v[3] - v0[3])^2)
            end
        end

        E0 /= w_total0
    end

    w_deleted = 0.0

    @inbounds for _ in 1:n_to_delete
        i_delete = floor(Int64, rand(rng, Float64) * pia.indexer[cell, species].n_local)
        i_delete = map_cont_index(pia.indexer[cell, species], i_delete)

        w_deleted += particles[i_delete].w
        delete_particle!(particles, pia, cell, species, i_delete)
    end

    scale_factor = 1.0

    # we most likely break continuity
    # and checking that we only deleted from group2 in the last domain cell
    # is kind of too much work for this merge
    @inbounds pia.contiguous[species] = false

    if conservative
        scale_factor = w_total0 / (w_total0 - w_deleted)
    else
        # first need to compute the remaining weight of the particles
        w_total_new = 0.0
        @inbounds s1 = pia.indexer[cell,species].start1
        @inbounds e1 = pia.indexer[cell,species].end1
        @inbounds for i in s1:e1
            w_total_new += particles[i].w
        end

        @inbounds if pia.indexer[cell, species].n_group2 > 0
            @inbounds s2 = pia.indexer[cell,species].start2
            @inbounds e2 = pia.indexer[cell,species].end2
            @inbounds for i in s2:e2
                w_total_new += particles[i].w
            end
        end
        scale_factor = (w_total_new + w_deleted) / w_total_new
    end

    # re-scale weights
    @inbounds s1 = pia.indexer[cell,species].start1
    @inbounds e1 = pia.indexer[cell,species].end1
    @inbounds for i in s1:e1
        particles[i].w *= scale_factor
    end

    @inbounds if pia.indexer[cell, species].n_group2 > 0
        @inbounds s2 = pia.indexer[cell,species].start2
        @inbounds e2 = pia.indexer[cell,species].end2
        @inbounds for i in s2:e2
            particles[i].w *= scale_factor
        end
    end

    # now we re-scale velocities if needed
    if conservative
        v_mean_new = SVector{3,Float64}([0.0, 0.0, 0.0])
        E_new = 0.0

        @inbounds s1 = pia.indexer[cell,species].start1
        @inbounds e1 = pia.indexer[cell,species].end1
        @inbounds for i in s1:e1
            v_mean_new += particles[i].w * particles[i].v
        end

        @inbounds if pia.indexer[cell, species].n_group2 > 0
            @inbounds s2 = pia.indexer[cell,species].start2
            @inbounds e2 = pia.indexer[cell,species].end2
            @inbounds for i in s2:e2
                v_mean_new += particles[i].w * particles[i].v
            end
        end
        v_mean_new /= w_total0

        @inbounds s1 = pia.indexer[cell,species].start1
        @inbounds e1 = pia.indexer[cell,species].end1
        @inbounds for i in s1:e1
            E_new += particles[i].w * ((particles[i].v[1] - v_mean_new[1])^2
                                       + (particles[i].v[2] - v_mean_new[2])^2
                                       + (particles[i].v[3] - v_mean_new[3])^2)
        end

        @inbounds if pia.indexer[cell, species].n_group2 > 0
            @inbounds s2 = pia.indexer[cell,species].start2
            @inbounds e2 = pia.indexer[cell,species].end2
            @inbounds for i in s2:e2
                E_new += particles[i].w * ((particles[i].v[1] - v_mean_new[1])^2
                                           + (particles[i].v[2] - v_mean_new[2])^2
                                           + (particles[i].v[3] - v_mean_new[3])^2)
            end
        end

        E_new /= w_total0

        E_scale = sqrt(E0 / E_new)

        @inbounds s1 = pia.indexer[cell,species].start1
        @inbounds e1 = pia.indexer[cell,species].end1
        @inbounds for i in s1:e1
            particles[i].v = E_scale*(particles[i].v - v_mean_new) + v0
        end

        @inbounds if pia.indexer[cell, species].n_group2 > 0
            @inbounds s2 = pia.indexer[cell,species].start2
            @inbounds e2 = pia.indexer[cell,species].end2
            @inbounds for i in s2:e2
                particles[i].v = E_scale*(particles[i].v - v_mean_new) + v0
            end
        end
    end
end