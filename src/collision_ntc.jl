mutable struct CollisionFactors
    n1::Int64
    n2::Int64
    sigma_g_w_max::Float64  # (σ(g) * g * w)_max
    n_coll::Int64  # number of collisions tested
    n_coll_performed::Int64  # number of collisions performed
    n_eq_w_coll_performed::Int64  # number of collisions of particles with equal weights (no splitting), for debugging/etc
end

function create_collision_factors()
    return CollisionFactors(0, 0.0, 0.0, 0, 0, 0)
end

function create_collision_factors(n_species)
    coll_factor_array = Array{CollisionFactors, 2}(undef, (n_species, n_species))
    for k in 1:n_species
        for i in 1:n_species
            coll_factor_array[i,k] = CollisionFactors(0, 0.0, 0.0, 0, 0, 0)
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

function ntc!(species, cell, rng, collision_factors, pia, collision_data, interaction, particles,
    Δt, V)
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

            sigma = sigma_vhs(interaction, collision_data.g)
            sigma_g_w_max = sigma * collision_data.g * max(particles[i].w, particles[k].w)

            # update (σ g w)_max if needed
            collision_factors.sigma_g_w_max = sigma_g_w_max > collision_factors.sigma_g_w_max ? sigma_g_w_max : collision_factors.sigma_g_w_max
            
            if (rand(rng, Float64) < sigma_g_w_max / collision_factors.sigma_g_w_max)
                collision_factors.n_coll_performed += 1
                compute_com!(collision_data, interaction, particles[i], particles[k])

                # do collision
                if (particles[i].w == particles[k].w)
                    collision_factors.n_eq_w_coll_performed += 1
                    scatter_vhs(rng, collision_data, interaction, particles[i], particles[k])
                elseif (particles[i].w > particles[k].w)
                    # we split particle i, update velocity of i and k (split part remains unchanged)

                    # first need to grow particle array
                    if (length(particles) <= pia.n_total[species])
                        resize!(particles, length(particles)+DELTA_PARTICLES)
                    end

                    # first need to update the particle indexer struct
                    update_particle_indexer_new_particle(species, cell, pia)

                    Δw = particles[i].w - particles[k].w
                    particles[i].w = particles[k].w

                    particles[pia.n_total[species]] = Particle(Δw, particles[i].v, particles[i].x)
                else  # (particles[k].w > particles[i].w)
                    if (length(particles) <= pia.n_total[species])
                        resize!(particles, length(particles)+DELTA_PARTICLES)
                    end

                    update_particle_indexer_new_particle(species, cell, pia)

                    Δw = particles[k].w - particles[i].w
                    particles[k].w = particles[i].w

                    particles[pia.n_total[species]] = Particle(Δw, particles[k].v, particles[k].x)
                end
                scatter_vhs(rng, collision_data, interaction, particles[i], particles[k])
            end
        end
    end
end



function ntc!(species1, species2, cell, rng, collision_factors, pia,
    collision_data, interaction, particles_1, particles_2,
    Δt, V)
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

            sigma = sigma_vhs(interaction, collision_data.g)
            sigma_g_w_max = sigma * collision_data.g * max(particles_1[i].w, particles_2[k].w)

            # update (σ g w)_max if needed
            collision_factors.sigma_g_w_max = sigma_g_w_max > collision_factors.sigma_g_w_max ? sigma_g_w_max : collision_factors.sigma_g_w_max

            if (rand(rng, Float64) < sigma_g_w_max / collision_factors.sigma_g_w_max)
                collision_factors.n_coll_performed += 1
                compute_com!(collision_data, interaction, particles_1[i], particles_2[k])
                # do collision
                if (particles_1[i].w == particles_2[k].w)
                    collision_factors.n_eq_w_coll_performed += 1
                    scatter_vhs(rng, collision_data, interaction, particles_1[i], particles_2[k])
                elseif (particles_1[i].w > particles_2[k].w)
                    # we split particle i, update velocity of i and k (split part remains unchanged)

                    # first need to update the particle indexer struct
                    update_particle_indexer_new_particle(species1, cell, pia)

                    Δw = particles_1[i].w - particles_2[k].w
                    particles_1[i].w = particles_2[k].w

                    particles_1[n_total].w = Δw
                    particles_1[n_total].v = particles_1[i].v
                    particles_1[n_total].x = particles_1[i].x
                else  # (particles[k].w > particles[i].w)
                    update_particle_indexer_new_particle(species2, cell, pia)

                    Δw = particles_2[k].w - particles_1[i].w
                    particles_2[k].w = particles_1[i].w

                    particles_2[n_total].w = Δw
                    particles_2[n_total].v = particles_2[k].v
                    particles_2[n_total].x = particles_2[k].x
                end
                scatter_vhs(rng, collision_data, interaction, particles_1[i], particles_2[k])
            end
        end
    end
end