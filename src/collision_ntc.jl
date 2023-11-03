mutable struct CollisionFactors
    n1::Int64
    n2::Int64
    sigma_g_w_max::Float64  # (σ(g) * g * w)_max
    n_coll::Int64
end

function compute_n_coll_single_species(rng, collision_factors, particle_indexer, Δt, V)
    return Δt * particle_indexer.n_total * (particle_indexer.n_total - 1) * collision_factors.sigma_g_w_max / V +
        rand(rng, Float64)
end

function ntc_single_species(rng, collision_factors, particle_indexer, collision_data, interaction, particles,
    Δt, V)
    # compute ncoll
    # loop over particles
    # update sigma_g_w_max
    # collide
    collision_factors.n_coll = compute_n_coll_single_species(rng, collision_factors, particle_indexer, Δt, V)
    n_coll_int = floor(Int64, collision_factors.n_coll)

    for _ in 1:n_coll_int
        # TODO: check if 0 or 1 !!!
        i = floor(Int64, rand(rng, Float64) * particle_indexer.n_total) + 1
        k = floor(Int64, rand(rng, Float64) * particle_indexer.n_total) + 1

        while (i == k)
            k = floor(Int64, rand(rng, Float64) * particle_indexer.n_total) + 1
        end

        # TODO: check bounds!
        i = i > particle_indexer.n_group1 ? (i - particle_indexer.n_group1) + particle_indexer.start2 : i + particle_indexer.start1
        # k = ...
    end
end