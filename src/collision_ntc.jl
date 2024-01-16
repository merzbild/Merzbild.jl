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

function compute_n_coll_single_species(rng, collision_factors, particle_indexer, Δt, V)
    return 0.5 * Δt * particle_indexer.n_total * (particle_indexer.n_total - 1) * collision_factors.sigma_g_w_max / V +
        rand(rng, Float64)
end

function compute_n_coll_two_species(rng, collision_factors, particle_indexer_1, particle_indexer_2, Δt, V)
    return Δt * particle_indexer_1.n_total * particle_indexer_2.n_total * collision_factors.sigma_g_w_max / V +
        rand(rng, Float64)
end

function ntc!(rng, collision_factors, particle_indexer, collision_data, interaction, particles,
    Δt, V)
    # compute ncoll
    # loop over particles
    # update sigma_g_w_max
    # collide
    collision_factors.n1 = particle_indexer.n_total
    collision_factors.n2 = particle_indexer.n_total
    n_coll_float = compute_n_coll_single_species(rng, collision_factors, particle_indexer, Δt, V)
    n_coll_int = floor(Int64, n_coll_float)
    # println(n_coll_float, ", ", n_coll_int)

    collision_factors.n_coll = n_coll_int
    collision_factors.n_coll_performed = 0
    collision_factors.n_eq_w_coll_performed = 0

    for _ in 1:n_coll_int
        # TODO: check if 0 or 1 !!!
        i = floor(Int64, rand(rng, Float64) * particle_indexer.n_total)
        k = floor(Int64, rand(rng, Float64) * particle_indexer.n_total)

        while (i == k)
            k = floor(Int64, rand(rng, Float64) * particle_indexer.n_total)
        end

        # example: bounds from [1,4], [7,9]; n_total = 7
        # n_group1 = 4, n_group2 = 3
        # i = 0,1,2,3 - [1,4]
        # i = 4,5,6 - [7,9]
        i = map_cont_index(particle_indexer, i)
        k = map_cont_index(particle_indexer, k)
        
        compute_g!(collision_data, particles[i], particles[k])

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

                    # first need to update the particle indexer struct
                    update_particle_indexer_new_particle(particle_indexer)

                    Δw = particles[i].w - particles[k].w
                    particles[i].w = particles[k].w

                    particles[n_total].w = Δw
                    particles[n_total].v = particles[i].v
                    particles[n_total].x = particles[i].x
                else  # (particles[k].w > particles[i].w)
                    update_particle_indexer_new_particle(particle_indexer)

                    Δw = particles[k].w - particles[i].w
                    particles[k].w = particles[i].w

                    particles[n_total].w = Δw
                    particles[n_total].v = particles[k].v
                    particles[n_total].x = particles[k].x
                end
                scatter_vhs(rng, collision_data, interaction, particles[i], particles[k])
            end
        end
    end
end



function ntc!(rng, collision_factors, particle_indexer_1, particle_indexer_2,
    collision_data, interaction, particles_1, particles_2,
    Δt, V)
    # compute ncoll
    # loop over particles
    # update sigma_g_w_max
    # collide
    collision_factors.n1 = particle_indexer_1.n_total
    collision_factors.n2 = particle_indexer_2.n_total
    n_coll_float = compute_n_coll_two_species(rng, collision_factors,
    particle_indexer_1, particle_indexer_2, Δt, V)
    n_coll_int = floor(Int64, n_coll_float)
    # println(n_coll_float, ", ", n_coll_int)

    collision_factors.n_coll = n_coll_int
    collision_factors.n_coll_performed = 0
    collision_factors.n_eq_w_coll_performed = 0

    for _ in 1:n_coll_int
        # TODO: check if 0 or 1 !!!
        i = floor(Int64, rand(rng, Float64) * particle_indexer_1.n_total)
        k = floor(Int64, rand(rng, Float64) * particle_indexer_2.n_total)

        # example: bounds from [1,4], [7,9]; n_total = 7
        # n_group1 = 4, n_group2 = 3
        # i = 0,1,2,3 - [1,4]
        # i = 4,5,6 - [7,9]
        i = map_cont_index(particle_indexer_1, i)
        k = map_cont_index(particle_indexer_2, k)
        
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
                    update_particle_indexer_new_particle(particle_indexer)

                    Δw = particles_1[i].w - particles_2[k].w
                    particles_1[i].w = particles_2[k].w

                    particles_1[n_total].w = Δw
                    particles_1[n_total].v = particles_1[i].v
                    particles_1[n_total].x = particles_1[i].x
                else  # (particles[k].w > particles[i].w)
                    update_particle_indexer_new_particle(particle_indexer)

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