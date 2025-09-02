@muladd begin

"""
    fp_linear!(rng, collision_data_fp, interaction, species_data, particles, pia, cell, species, Δt, V)

Model single-species elastic collisions using a linear Fokker-Planck approximation.

# Positional arguments
* `rng`: the random number generator
* `collision_data_fp`: `CollisionDataFP` instance used for storing collisional quantities
* `interaction`: 2-dimensional array of `Interaction` instances for all possible species pairs
* `species_data`: the vector of `SpeciesData`
* `particles`: `ParticleVector` of the particles being collided
* `pia`: the `ParticleIndexerArray`
* `cell`: the index of the cell in which collisions are performed
* `species`: the index of the species for which collisions are performed
* `Δt`: timestep
* `V`: cell volume

# References
* M.H. Gorji, M. Torrilhon, P. Jenny, Fokker–Planck model for computational studies of monatomic rarefied gas flows.
    [J. Fluid Mech., 2011](https://doi.org/10.1017/jfm.2011.188).
"""
function fp_linear!(rng, collision_data_fp, interaction, species_data, particles, pia, cell, species, Δt, V)
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

    τ = compute_relaxation_time(interaction, species_data, species, V, es_old, local_w)
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

"""
    compute_relaxation_time(interaction, species_data, species, V, es_old, local_w)

Compute the relaxation time for a single-species flow, as given by ``\\tau = 2\\mu / p``.

# Positional arguments
* `interaction`: the `Interaction` instance for the self-interaction of the species in question
* `species_data`: the vector of `SpeciesData`
* `species`: the index of the species for which collisions are performed
* `V`: cell volume
* `es_old`: the total kinetic energy of the particles in the cell, divided by the species' molecular mass
* `local_w`: the total computational weight of the particles in the cell

# Returns
The elastic collision relaxation time.
"""
@inline function compute_relaxation_time(interaction, species_data, species, V, es_old, local_w)
    T = es_old * species_data[species].mass / ((3.0 / 2.0) * k_B)

    nrho = local_w / V
    p = nrho * k_B * T
    μ = interaction.vhs_muref * (T / interaction.vhs_Tref)^interaction.vhs_o

    return 2.0 * μ / p;
end

"""
    sample_normal_rands!(rng, collision_data_fp, n_local)

Sample `3*n_local` normally distributed random numbers (with mean 0 and variance 1), one for each velocity component,
and write them to a `CollisionDataFP` instance.

# Positional arguments
* `rng`: the random number generator
* `collision_data_fp`: `CollisionDataFP` instance used for storing collisional quantities
* `n_local`: the number of particles to sample the numbers for
"""
@inline function sample_normal_rands!(rng, collision_data_fp, n_local)
    for i in 1:n_local
        @inbounds collision_data_fp.xvel_rand[i] = randn(rng)
        @inbounds collision_data_fp.yvel_rand[i] = randn(rng)
        @inbounds collision_data_fp.zvel_rand[i] = randn(rng)
    end
end

"""
    scale_norm_rands!(collision_data_fp, n_local)

Scale sampled normally distributed random numbers velocity component-wise, so that their
means for each component are exactly 0, and their variances are exactly 1.

# Positional arguments
* `collision_data_fp`: `CollisionDataFP` instance used for storing collisional quantities
* `n_local`: the number of particles for which the numbers were sampled for
"""
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