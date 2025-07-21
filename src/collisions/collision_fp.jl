@muladd begin

"""
Model single-species elastic collisions using a linear Fokker-Planck approximation 
"""
function fp!(rng, collision_data_fp, interaction, species_data, particles, pia, cell, species, Δt, V)
    @inbounds indexer = pia.indexer[cell, species]
    n_local = indexer.n_local
    n_begin = indexer.start1
    n_end = indexer.end1

    if n_local < 7
        return nothing
    end

    local_w = 0.0
    es_old = 0.0
    es_new = 0.0
    Krot = 0.0
    Kvib = 0.0
    collision_data_fp.vel_ave = SVector{3,Float64}(0.0, 0.0, 0.0)

    if length(collision_data_fp.xvel_rand) < n_local
        resize!(collision_data_fp.xvel_rand, n_local+DELTA_PARTICLES)
        resize!(collision_data_fp.yvel_rand, n_local+DELTA_PARTICLES)
        resize!(collision_data_fp.zvel_rand, n_local+DELTA_PARTICLES)
    end

    sample_normal_rands!(rng, collision_data_fp, n_local)
    scale_norm_rands!(collision_data_fp, n_local)

    for part_id in n_begin:n_end
        @inbounds collision_data_fp.vel_ave = collision_data_fp.vel_ave + particles[part_id].v * particles[part_id].w
        @inbounds local_w += particles[part_id].w
    end

    if indexer.n_group2 > 0
        for part_id in indexer.start2:indexer.end2
            @inbounds collision_data_fp.vel_ave = vel_ave + particles[part_id].v * particles[part_id].w
            @inbounds local_w += particles[part_id].w
        end
    end

    collision_data_fp.vel_ave = collision_data_fp.vel_ave / local_w

    for part_id in n_begin:n_end
        @inbounds particles[part_id].v = particles[part_id].v - collision_data_fp.vel_ave
        @inbounds es_old = es_old + 0.5 * (particles[part_id].v[1]^2
                                  + particles[part_id].v[2]^2
                                  + particles[part_id].v[3]^2) * particles[part_id].w
    end

    if indexer.n_group2 > 0
        for part_id in indexer.start2:indexer.end2
            @inbounds particles[part_id].v = particles[part_id].v - collision_data_fp.vel_ave
            @inbounds es_old = es_old + 0.5 * (particles[part_id].v[1]^2
                                      + particles[part_id].v[2]^2
                                      + particles[part_id].v[3]^2) * particles[part_id].w
        end
    end

    es_old = es_old / local_w

    τ = compute_relaxation_time(interaction, species_data, particles, pia, cell, species, V, es_old, local_w)
    A = exp(-Δt/τ)
    B = (1.0 - A)*τ

    C = sqrt(((2.0 / 3.0) * es_old + (Krot + Kvib) * τ) * (1.0 - exp(-2.0 * Δt/τ)))

    for part_id in n_begin:n_end
        i = part_id - n_begin + 1

        @inbounds particles[part_id].v = particles[part_id].v*A + C*SVector{3,Float64}(collision_data_fp.xvel_rand[i],
                                                                             collision_data_fp.yvel_rand[i],
                                                                             collision_data_fp.zvel_rand[i])
        @inbounds es_new = es_new + 0.5 * (particles[part_id].v[1]^2
                                  + particles[part_id].v[2]^2
                                  + particles[part_id].v[3]^2) * particles[part_id].w
    end

    if indexer.n_group2 > 0
        for part_id in indexer.start2:indexer.end2
            i = part_id - n_begin + 1

            @inbounds particles[part_id].v = particles[part_id].v*A + C*SVector{3,Float64}(collision_data_fp.xvel_rand[i],
                                                                                 collision_data_fp.yvel_rand[i],
                                                                                 collision_data_fp.zvel_rand[i])
            @inbounds es_new = es_new + 0.5 * (particles[part_id].v[1]^2
                                      + particles[part_id].v[2]^2
                                      + particles[part_id].v[3]^2) * particles[part_id].w
        end
    end

    es_new = es_new / local_w

    alpha = sqrt(es_old / es_new)

    for part_id in n_begin:n_end
        @inbounds particles[part_id].v = alpha*particles[part_id].v + collision_data_fp.vel_ave
    end

    if indexer.n_group2 > 0
        for part_id in indexer.start2:indexer.end2
            @inbounds particles[part_id].v = alpha*particles[part_id].v + collision_data_fp.vel_ave
        end
    end
end

@inline function compute_relaxation_time(interaction, species_data, particles, pia, cell, species, V, es_old, local_w)
    T = es_old * species_data[species].mass / ((3.0 / 2.0) * k_B)

    nrho = local_w / V
    p = nrho * k_B * T
    μ = interaction.vhs_muref * (T / interaction.vhs_Tref)^interaction.vhs_o

    return 2.0 * μ / p;
end

@inline function sample_normal_rands!(rng, collision_data_fp, n_local)
    for i in 1:n_local
        @inbounds collision_data_fp.xvel_rand[i] = randn(rng)
        @inbounds collision_data_fp.yvel_rand[i] = randn(rng)
        @inbounds collision_data_fp.zvel_rand[i] = randn(rng)
    end
end

@inline function scale_norm_rands!(collision_data_fp, n_local)
    collision_data_fp.mean = SVector{3,Float64}(0.0, 0.0, 0.0)
    collision_data_fp.stddev = SVector{3,Float64}(0.0, 0.0, 0.0)

    for i in 1:n_local
        @inbounds collision_data_fp.mean = collision_data_fp.mean + SVector{3,Float64}(collision_data_fp.xvel_rand[i],
                                                                             collision_data_fp.yvel_rand[i],
                                                                             collision_data_fp.zvel_rand[i])
    end

    collision_data_fp.mean = collision_data_fp.mean / n_local

    for i in 1:n_local
        @inbounds collision_data_fp.xvel_rand[i] -= collision_data_fp.mean[1]
        @inbounds collision_data_fp.yvel_rand[i] -= collision_data_fp.mean[2]
        @inbounds collision_data_fp.zvel_rand[i] -= collision_data_fp.mean[3]

        @inbounds collision_data_fp.stddev = collision_data_fp.stddev + SVector{3,Float64}(collision_data_fp.xvel_rand[i]^2,
                                                                                 collision_data_fp.yvel_rand[i]^2,
                                                                                 collision_data_fp.zvel_rand[i]^2)
    end

    collision_data_fp.stddev = SVector{3,Float64}(sqrt.(n_local ./ collision_data_fp.stddev))

    for i in 1:n_local
        @inbounds collision_data_fp.xvel_rand[i] *= collision_data_fp.stddev[1];
        @inbounds collision_data_fp.yvel_rand[i] *= collision_data_fp.stddev[2];
        @inbounds collision_data_fp.zvel_rand[i] *= collision_data_fp.stddev[3];
    end
end

end