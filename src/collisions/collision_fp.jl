@muladd begin

"""
Model single-species elastic collisions using a linear Fokker-Planck approximation 
"""
function fp!(rng, collision_data_fp, interaction, species_data, particles, pia, cell, species, Δt, V)
    indexer = pia.indexer[cell, species]
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

    rands = randn(rng, Float64, (n_local, 3))
    scale_norm_rands!(collision_data_fp, rands)

    for part_id in n_begin:n_end
        collision_data_fp.vel_ave = collision_data_fp.vel_ave + particles[part_id].v * particles[part_id].w
        local_w += particles[part_id].w
    end

    if indexer.n_group2 > 0
        for part_id in indexer.start2:indexer.end2
            collision_data_fp.vel_ave = vel_ave + particles[part_id].v * particles[part_id].w
            local_w += particles[part_id].w
        end
    end

    collision_data_fp.vel_ave = collision_data_fp.vel_ave / local_w

    for part_id in n_begin:n_end
        particles[part_id].v = particles[part_id].v - collision_data_fp.vel_ave
        es_old = es_old + 0.5 * (particles[part_id].v[1]^2
                                + particles[part_id].v[2]^2
                                + particles[part_id].v[3]^2) * particles[part_id].w
    end

    if indexer.n_group2 > 0
        for part_id in indexer.start2:indexer.end2
            particles[part_id].v = particles[part_id].v - collision_data_fp.vel_ave
            es_old = es_old + 0.5 * (particles[part_id].v[1]^2
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

        particles[part_id].v = particles[part_id].v*A + C*rands[i, :]
        es_new = es_new + 0.5 * (particles[part_id].v[1]^2
                                 + particles[part_id].v[2]^2
                                 + particles[part_id].v[3]^2) * particles[part_id].w
    end

    if indexer.n_group2 > 0
        for part_id in indexer.start2:indexer.end2
            i = part_id - n_begin + 1

            particles[part_id].v = particles[part_id].v*A + C*rands[i, :]
            es_new = es_new + 0.5 * (particles[part_id].v[1]^2
                                    + particles[part_id].v[2]^2
                                    + particles[part_id].v[3]^2) * particles[part_id].w
        end
    end

    es_new = es_new / local_w

    alpha = sqrt(es_old / es_new)

    for part_id in n_begin:n_end
        particles[part_id].v = alpha*particles[part_id].v + collision_data_fp.vel_ave
    end

    if indexer.n_group2 > 0
        for part_id in indexer.start2:indexer.end2
            particles[part_id].v = alpha*particles[part_id].v + collision_data_fp.vel_ave
        end
    end

    return nothing
end

@inline function compute_relaxation_time(interaction, species_data, particles, pia, cell, species, V, es_old, local_w)
    T = es_old * species_data[species].mass / ((3.0 / 2.0) * k_B)

    nrho = local_w / V
    p = nrho * k_B * T
    μ = interaction.vhs_muref * (T / interaction.vhs_Tref)^interaction.vhs_o

    return 2.0 * μ / p;
end

function scale_norm_rands!(collision_data_fp, norm_rands)
    n_length = size(norm_rands)[1]
    collision_data_fp.mean = SVector{3,Float64}(0.0, 0.0, 0.0)
    collision_data_fp.stddev = SVector{3,Float64}(0.0, 0.0, 0.0)

    for i in 1:n_length
        collision_data_fp.mean = collision_data_fp.mean + norm_rands[i, :]
    end
    collision_data_fp.mean = collision_data_fp.mean / n_length

    for i in 1:n_length
        norm_rands[i, 1] -= collision_data_fp.mean[1]
        norm_rands[i, 2] -= collision_data_fp.mean[2]
        norm_rands[i, 3] -= collision_data_fp.mean[3]

        collision_data_fp.stddev = collision_data_fp.stddev + SVector{3,Float64}(norm_rands[i, 1]^2, norm_rands[i, 2]^2, norm_rands[i, 3]^2)
    end

    collision_data_fp.stddev = 1.0 ./ sqrt.(collision_data_fp.stddev / n_length)

    for i in 1:n_length
        norm_rands[i, 1] *= collision_data_fp.stddev[1];
        norm_rands[i, 2] *= collision_data_fp.stddev[2];
        norm_rands[i, 3] *= collision_data_fp.stddev[3];
    end
end

end