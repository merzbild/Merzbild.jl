mutable struct CollisionFactors
    n1::Int64
    n2::Int64
    sigma_g_w_max::Float64  # (σ(g) * g * w)_max
    n_coll::Int64  # number of collisions tested
    n_coll_performed::Int64  # number of collisions performed
    n_eq_w_coll_performed::Int64  # number of collisions of particles with equal weights (no splitting), for debugging/etc
end

CollisionFactors() = CollisionFactors(0, 0.0, 0.0, 0, 0, 0)

function create_collision_factors_array(n_species)
    coll_factor_array = Array{CollisionFactors, 3}(undef, (n_species, n_species, 1))
    for k in 1:n_species
        for i in 1:n_species
            coll_factor_array[i,k,1] = CollisionFactors(0, 0.0, 0.0, 0, 0, 0)
        end
    end
    return coll_factor_array
end

function create_collision_factors_array(n_species, n_cells)
    coll_factor_array = Array{CollisionFactors, 3}(undef, (n_species, n_species, n_cells))
    for k in 1:n_cells
        for j in 1:n_species
            for i in 1:n_species
                coll_factor_array[i,j,k] = CollisionFactors(0, 0.0, 0.0, 0, 0, 0)
            end
        end
    end
    return coll_factor_array
end

function compute_n_coll_single_species(rng, collision_factors, np, Δt, V)
    return 0.5 * Δt * np * (np - 1) * collision_factors.sigma_g_w_max / V +
        rand(rng, Float64)
end

function compute_n_coll_two_species(rng, collision_factors, np1, np2, Δt, V)
    return Δt * np1 * np2 * collision_factors.sigma_g_w_max / V + rand(rng, Float64)
end

function ntc!(rng, collision_factors, collision_data, interaction, particles, pia,
              cell, species, Δt, V)
    # single-species ntc
    # compute ncoll
    # loop over particles
    # update sigma_g_w_max
    # collide
    collision_factors.n1 = pia.indexer[cell, species].n_local
    collision_factors.n2 = pia.indexer[cell, species].n_local
    n_coll_float = compute_n_coll_single_species(rng, collision_factors, pia.indexer[cell, species].n_local, Δt, V)
    n_coll_int = floor(Int64, n_coll_float)
    # println(n_coll_float, ", ", n_coll_int)

    collision_factors.n_coll = n_coll_int
    collision_factors.n_coll_performed = 0
    collision_factors.n_eq_w_coll_performed = 0

    for _ in 1:n_coll_int
        # TODO: check if 0 or 1 !!!
        i = floor(Int64, rand(rng, Float64) * pia.indexer[cell, species].n_local)
        k = floor(Int64, rand(rng, Float64) * pia.indexer[cell, species].n_local)

        while (i == k)
            k = floor(Int64, rand(rng, Float64) * pia.indexer[cell, species].n_local)
        end

        # example: bounds from [1,4], [7,9]; n_total = 7
        # n_group1 = 4, n_group2 = 3
        # i = 0,1,2,3 - [1,4]
        # i = 4,5,6 - [7,9]
        i = map_cont_index(pia.indexer[cell, species], i)
        k = map_cont_index(pia.indexer[cell, species], k)
        
        compute_g!(collision_data, particles[i], particles[k])
        # println("NTC: ", particle_indexer.n_total, ", ", length(particles))
        if (collision_data.g > eps())

            sigma = sigma_vhs(interaction[species, species], collision_data.g)
            sigma_g_w_max = sigma * collision_data.g * max(particles[i].w, particles[k].w)

            # update (σ g w)_max if needed
            collision_factors.sigma_g_w_max = sigma_g_w_max > collision_factors.sigma_g_w_max ? sigma_g_w_max : collision_factors.sigma_g_w_max
            
            if (rand(rng, Float64) < sigma_g_w_max / collision_factors.sigma_g_w_max)
                collision_factors.n_coll_performed += 1
                compute_com!(collision_data, interaction[species, species], particles[i], particles[k])
                # do collision
                if (particles[i].w == particles[k].w)
                    collision_factors.n_eq_w_coll_performed += 1
                elseif (particles[i].w > particles[k].w)
                    # we split particle i, update velocity of i and k (split part remains unchanged)

                    # first need to grow particle array
                    if (length(particles) <= pia.n_total[species])
                        resize!(particles, length(particles)+DELTA_PARTICLES)
                    end

                    # first need to update the particle indexer struct
                    update_particle_indexer_new_particle(pia, cell, species)

                    Δw = particles[i].w - particles[k].w
                    particles[i].w = particles[k].w

                    particles[pia.n_total[species]] = Particle(Δw, particles[i].v, particles[i].x)
                else  # (particles[k].w > particles[i].w)
                    if (length(particles) <= pia.n_total[species])
                        resize!(particles, length(particles)+DELTA_PARTICLES)
                    end

                    update_particle_indexer_new_particle(pia, cell, species)

                    Δw = particles[k].w - particles[i].w
                    particles[k].w = particles[i].w

                    particles[pia.n_total[species]] = Particle(Δw, particles[k].v, particles[k].x)
                end
                scatter_vhs!(rng, collision_data, interaction[species, species], particles[i], particles[k])
            end
        end
    end
end

function ntc!(rng, collision_factors, collision_data, interaction,
              particles_1, particles_2, pia,
              cell, species1, species2, Δt, V)
    # compute ncoll
    # loop over particles
    # update sigma_g_w_max
    # collide
    collision_factors.n1 = pia.indexer[cell, species1].n_local
    collision_factors.n2 = pia.indexer[cell, species2].n_local
    n_coll_float = compute_n_coll_two_species(rng, collision_factors,
                                              pia.indexer[cell, species1].n_local, pia.indexer[cell, species2].n_local, Δt, V)
    n_coll_int = floor(Int64, n_coll_float)
    # println(n_coll_float, ", ", n_coll_int)

    collision_factors.n_coll = n_coll_int
    collision_factors.n_coll_performed = 0
    collision_factors.n_eq_w_coll_performed = 0

    for _ in 1:n_coll_int
        # TODO: check if 0 or 1 !!!
        i = floor(Int64, rand(rng, Float64) * pia.indexer[cell, species1].n_local)
        k = floor(Int64, rand(rng, Float64) * pia.indexer[cell, species2].n_local)

        # example: bounds from [1,4], [7,9]; n_total = 7
        # n_group1 = 4, n_group2 = 3
        # i = 0,1,2,3 - [1,4]
        # i = 4,5,6 - [7,9]
        i = map_cont_index(pia.indexer[cell, species1], i)
        k = map_cont_index(pia.indexer[cell, species2], k)
        
        compute_g!(collision_data, particles_1[i], particles_2[k])

        if (collision_data.g > eps())

            sigma = sigma_vhs(interaction[species1, species2], collision_data.g)
            sigma_g_w_max = sigma * collision_data.g * max(particles_1[i].w, particles_2[k].w)

            # update (σ g w)_max if needed
            collision_factors.sigma_g_w_max = sigma_g_w_max > collision_factors.sigma_g_w_max ? sigma_g_w_max : collision_factors.sigma_g_w_max

            if (rand(rng, Float64) < sigma_g_w_max / collision_factors.sigma_g_w_max)
                collision_factors.n_coll_performed += 1
                compute_com!(collision_data, interaction[species1, species2], particles_1[i], particles_2[k])
                # do collision
                if (particles_1[i].w == particles_2[k].w)
                    collision_factors.n_eq_w_coll_performed += 1
                elseif (particles_1[i].w > particles_2[k].w)
                    # we split particle i, update velocity of i and k (split part remains unchanged)
                    if (length(particles_1) <= pia.n_total[species_1])
                        resize!(particles_1, length(particles_1)+DELTA_PARTICLES)
                    end
                    # first need to update the particle indexer struct
                    update_particle_indexer_new_particle(pia, cell, species1)

                    Δw = particles_1[i].w - particles_2[k].w
                    particles_1[i].w = particles_2[k].w

                    particles_1[pia.n_total[species_1]] = Particle(Δw, particles_1[i].v, particles_1[i].x)
                else  # (particles[k].w > particles[i].w)
                    if (length(particles_2) <= pia.n_total[species_2])
                        resize!(particles_2, length(particles_2)+DELTA_PARTICLES)
                    end
                    update_particle_indexer_new_particle(pia, cell, species2)

                    Δw = particles_2[k].w - particles_1[i].w
                    particles_2[k].w = particles_1[i].w

                    particles_2[pia.n_total[species_2]] = Particle(Δw, particles_2[i].v, particles_2[i].x)
                end
                scatter_vhs!(rng, collision_data, interaction[species1, species2], particles_1[i], particles_2[k])
            end
        end
    end
end


function estimate_sigma_g_w_max_ntc_n_e!(rng, collision_factors, collision_data, interaction,
                                         n_e_interactions, n_e_cs, particles_n, particles_e,
                                         pia, cell, species_n, species_e, Δt, V; min_coll=5, n_loops=3)
    # compute ncoll
    # loop over particles
    # update sigma_g_w_max
    # collide
    collision_factors.n1 = pia.indexer[cell, species_n].n_local
    collision_factors.n2 = pia.indexer[cell, species_e].n_local

    for _ in 1:n_loops
        n_coll_float = compute_n_coll_two_species(rng, collision_factors,
                                                  pia.indexer[cell, species_n].n_local,
                                                  pia.indexer[cell, species_e].n_local, Δt, V) + min_coll
        n_coll_int = floor(Int64, n_coll_float)

        collision_factors.n_coll = n_coll_int
        collision_factors.n_coll_performed = 0
        collision_factors.n_eq_w_coll_performed = 0

        for _ in 1:n_coll_int
            # TODO: check if 0 or 1 !!!
            i = floor(Int64, rand(rng, Float64) * pia.indexer[cell, species_n].n_local)
            k = floor(Int64, rand(rng, Float64) * pia.indexer[cell, species_e].n_local)

            # example: bounds from [1,4], [7,9]; n_total = 7
            # n_group1 = 4, n_group2 = 3
            # i = 0,1,2,3 - [1,4]
            # i = 4,5,6 - [7,9]
            i = map_cont_index(pia.indexer[cell, species_n], i)
            k = map_cont_index(pia.indexer[cell, species_e], k)
            
            compute_g!(collision_data, particles_n[i], particles_e[k])


            if (collision_data.g > eps())
                collision_data.E_coll_eV = compute_cross_sections!(n_e_cs, interaction[species_n, species_e], collision_data.g, n_e_interactions, species_n)
                sigma_g_w_max = n_e_cs[species_n].cs_total * collision_data.g * max(particles_n[i].w, particles_e[k].w)

                # update (σ g w)_max if needed
                collision_factors.sigma_g_w_max = sigma_g_w_max > collision_factors.sigma_g_w_max ? sigma_g_w_max : collision_factors.sigma_g_w_max
            end
        end
    end
end

function ntc_n_e!(rng, collision_factors, collision_data, interaction,
                  n_e_interactions, n_e_cs, particles_n, particles_e, particles_ion,
                  pia, cell, species_n, species_e, species_ion, Δt, V)
    # no event splitting
    # compute ncoll
    # loop over particles
    # update sigma_g_w_max
    # collide
    collision_factors.n1 = pia.indexer[cell, species_n].n_local
    collision_factors.n2 = pia.indexer[cell, species_e].n_local
    n_coll_float = compute_n_coll_two_species(rng, collision_factors,
                                              pia.indexer[cell, species_n].n_local, pia.indexer[cell, species_e].n_local, Δt, V)
    n_coll_int = floor(Int64, n_coll_float)
    # println(n_coll_float, ", ", n_coll_int)

    collision_factors.n_coll = n_coll_int
    collision_factors.n_coll_performed = 0
    collision_factors.n_eq_w_coll_performed = 0

    for _ in 1:n_coll_int
        i = floor(Int64, rand(rng, Float64) * pia.indexer[cell, species_n].n_local)
        k = floor(Int64, rand(rng, Float64) * pia.indexer[cell, species_e].n_local)

        # example: bounds from [1,4], [7,9]; n_total = 7
        # n_group1 = 4, n_group2 = 3
        # i = 0,1,2,3 - [1,4]
        # i = 4,5,6 - [7,9]
        i = map_cont_index(pia.indexer[cell, species_n], i)
        k = map_cont_index(pia.indexer[cell, species_e], k)
        
        compute_g!(collision_data, particles_n[i], particles_e[k])

        if (collision_data.g > eps())

            collision_data.E_coll_electron_eV = compute_cross_sections!(n_e_cs, interaction[species_n, species_e], collision_data.g, n_e_interactions, species_n)
            sigma_g_w_max = get_cs_total(n_e_interactions, n_e_cs, species_n) * collision_data.g * max(particles_n[i].w, particles_e[k].w)

            # update (σ g w)_max if needed
            collision_factors.sigma_g_w_max = sigma_g_w_max > collision_factors.sigma_g_w_max ? sigma_g_w_max : collision_factors.sigma_g_w_max

            # if (collision_data.E_coll_electron_eV > 1500)
            #     println(collision_factors.sigma_g_w_max, ", ", i, ", ", k, ", ", particles_n[i].w, ", ", particles_e[k].w)
            # end

            if (rand(rng, Float64) < sigma_g_w_max / collision_factors.sigma_g_w_max)
                collision_factors.n_coll_performed += 1
                compute_com!(collision_data, interaction[species_n, species_e], particles_n[i], particles_e[k])

                # do collision

                if (particles_n[i].w > particles_e[k].w)
                    # we split particle i, update velocity of i and k (split part remains unchanged)
                    if (length(particles_n) <= pia.n_total[species_n])
                        resize!(particles_n, length(particles_n)+DELTA_PARTICLES)
                    end
                    # first need to update the particle indexer struct
                    update_particle_indexer_new_particle(pia, cell, species_n)

                    Δw = particles_n[i].w - particles_e[k].w
                    particles_n[i].w = particles_e[k].w

                    particles_n[pia.n_total[species_n]] = Particle(Δw, particles_n[i].v, particles_n[i].x)
                elseif (particles_n[i].w == particles_e[k].w)
                    collision_factors.n_eq_w_coll_performed += 1
                else  # (particles[k].w > particles[i].w)
                    if (length(particles_e) <= pia.n_total[species_e])
                        resize!(particles_e, length(particles_e)+DELTA_PARTICLES)
                    end
                    update_particle_indexer_new_particle(pia, cell, species_e)

                    Δw = particles_e[k].w - particles_n[i].w
                    particles_e[k].w = particles_n[i].w

                    particles_e[pia.n_total[species_e]] = Particle(Δw, particles_e[k].v, particles_e[k].x)
                end

                # now we collide the 2 equal-weight particles
                R = rand(rng, Float64)  # decide which process we do

                # println(n_e_cs[species_n].prob_vec, ", ", n_e_cs[species_n].cs_total, ", ", n_e_cs[species_n].cs_elastic)
                if (R < n_e_cs[species_n].prob_vec[1])  # TODO: fix indexing! 
                    # elastic collision
                    scatter_vhs!(rng, collision_data, interaction[species_n, species_e], particles_n[i], particles_e[k])
                else
                    # perform the ionization: 
                    if (length(particles_ion) <= pia.n_total[species_ion])
                        resize!(particles_ion, length(particles_ion)+DELTA_PARTICLES)
                    end
                    # create the ion particle
                    update_particle_indexer_new_particle(pia, cell, species_ion)
                    particles_ion[pia.n_total[species_ion]] = Particle(particles_n[i].w, particles_n[i].v, particles_n[i].x)

                    # add a second electron
                    if (length(particles_e) <= pia.n_total[species_e])
                        resize!(particles_e, length(particles_e)+DELTA_PARTICLES)
                    end
                    update_particle_indexer_new_particle(pia, cell, species_e)
                    particles_e[pia.n_total[species_e]] = Particle(particles_n[i].w, particles_e[k].v, particles_e[k].x)

                    # set neutral particle weight to 0
                    particles_n[i].w = 0.0

                    # compute energy split across the primare and secondary electrons
                    compute_g_new_ionization!(collision_data, interaction[species_n, species_e],
                                              get_ionization_threshold(n_e_interactions, species_n), get_electron_energy_split(n_e_interactions, species_n))

                    scatter_ionization_electrons!(rng, collision_data, particles_e, k, pia.n_total[species_e])
                end
            end
        end
    end
end

function ntc_n_e_es!(rng, collision_factors, collision_data, interaction,
    n_e_interactions, n_e_cs, particles_n, particles_e, particles_ion, 
    pia, cell, species_n, species_e, species_ion, Δt, V)
    # event splitting
    # compute ncoll
    # loop over particles
    # update sigma_g_w_max
    # collide
    collision_factors.n1 = pia.indexer[cell, species_n].n_local
    collision_factors.n2 = pia.indexer[cell, species_e].n_local
    n_coll_float = compute_n_coll_two_species(rng, collision_factors,
                                              pia.indexer[cell, species_n].n_local, pia.indexer[cell, species_e].n_local, Δt, V)
    n_coll_int = floor(Int64, n_coll_float)
    # println(n_coll_float, ", ", n_coll_int)

    collision_factors.n_coll = n_coll_int
    collision_factors.n_coll_performed = 0
    collision_factors.n_eq_w_coll_performed = 0

    for _ in 1:n_coll_int
        i = floor(Int64, rand(rng, Float64) * pia.indexer[cell, species_n].n_local)
        k = floor(Int64, rand(rng, Float64) * pia.indexer[cell, species_e].n_local)

        # example: bounds from [1,4], [7,9]; n_total = 7
        # n_group1 = 4, n_group2 = 3
        # i = 0,1,2,3 - [1,4]
        # i = 4,5,6 - [7,9]
        i = map_cont_index(pia.indexer[cell, species_n], i)
        k = map_cont_index(pia.indexer[cell, species_e], k)
        
        compute_g!(collision_data, particles_n[i], particles_e[k])

        if (collision_data.g > eps())

            collision_data.E_coll_electron_eV = compute_cross_sections!(n_e_cs, interaction[species_n, species_e], collision_data.g, n_e_interactions, species_n)
            sigma_g_w_max = get_cs_total(n_e_interactions, n_e_cs, species_n) * collision_data.g * max(particles_n[i].w, particles_e[k].w)

            # update (σ g w)_max if needed
            collision_factors.sigma_g_w_max = sigma_g_w_max > collision_factors.sigma_g_w_max ? sigma_g_w_max : collision_factors.sigma_g_w_max

            # if (collision_data.E_coll_electron_eV > 1500)
            #     println(collision_factors.sigma_g_w_max, ", ", i, ", ", k, ", ", particles_n[i].w, ", ", particles_e[k].w)
            # end

            if (rand(rng, Float64) < sigma_g_w_max / collision_factors.sigma_g_w_max)
                collision_factors.n_coll_performed += 1
                compute_com!(collision_data, interaction[species_n, species_e], particles_n[i], particles_e[k])

                # do collision

                if (particles_n[i].w > particles_e[k].w)
                    # we split particle i, update velocity of i and k (split part remains unchanged)
                    if (length(particles_n) <= pia.n_total[species_n])
                        resize!(particles_n, length(particles_n)+DELTA_PARTICLES)
                    end
                    # first need to update the particle indexer struct
                    update_particle_indexer_new_particle(pia, cell, species_n)

                    Δw = particles_n[i].w - particles_e[k].w
                    particles_n[i].w = particles_e[k].w

                    particles_n[pia.n_total[species_n]] = Particle(Δw, particles_n[i].v, particles_n[i].x)
                elseif (particles_n[i].w == particles_e[k].w)
                    collision_factors.n_eq_w_coll_performed += 1
                else  # (particles[k].w > particles[i].w)
                    if (length(particles_e) <= pia.n_total[species_e])
                        resize!(particles_e, length(particles_e)+DELTA_PARTICLES)
                    end
                    update_particle_indexer_new_particle(pia, cell, species_e)

                    Δw = particles_e[k].w - particles_n[i].w
                    particles_e[k].w = particles_n[i].w

                    particles_e[pia.n_total[species_e]] = Particle(Δw, particles_e[k].v, particles_e[k].x)
                end

                # now we collide the 2 equal-weight particles
                
                if (n_e_cs[species_n].prob_vec[2] == 0.0)
                    scatter_vhs!(rng, collision_data, interaction[species_n, species_e], particles_n[i], particles_e[k])
                else
                    w_ionized = particles_e[k].w * n_e_cs[species_n].prob_vec[2]
    
                    particles_n[i].w -= w_ionized
                    particles_e[k].w -= w_ionized

                    if (length(particles_ion) <= pia.n_total[species_ion])
                        resize!(particles_ion, length(particles_ion)+DELTA_PARTICLES)
                    end
                    # create the ion particle
                    update_particle_indexer_new_particle(pia, cell, species_ion)
                    particles_ion[pia.n_total[species_ion]] = Particle(w_ionized, particles_n[i].v, particles_n[i].x)

                    # add 2 electrons (split + secondary)
                    if (length(particles_e) < pia.n_total[species_e] + 2)
                        resize!(particles_e, length(particles_e)+DELTA_PARTICLES)
                    end
                    update_particle_indexer_new_particle(pia, cell, species_e)
                    particles_e[pia.n_total[species_e]] = Particle(w_ionized, particles_e[k].v, particles_e[k].x)
                    k1 = pia.n_total[species_e]
                    update_particle_indexer_new_particle(pia, cell, species_e)
                    particles_e[pia.n_total[species_e]] = Particle(w_ionized, particles_e[k].v, particles_e[k].x)
                    k2 = pia.n_total[species_e]

                    # elastic scattering
                    scatter_vhs!(rng, collision_data, interaction[species_n, species_e], particles_n[i], particles_e[k])

                    compute_g_new_ionization!(collision_data, interaction[species_n, species_e],
                                              get_ionization_threshold(n_e_interactions, species_n), get_electron_energy_split(n_e_interactions, species_n))

                    scatter_ionization_electrons!(rng, collision_data, particles_e, k1, k2)
                end
            end
        end
    end
end