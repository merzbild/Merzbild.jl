# using NNLS
# using TimerOutputs

mutable struct NNLSMerging
    v0::SVector{3,Float64}
    vref::Float64  # used for scaling
    inv_vref::Float64  # used for scaling

    Ex::Float64  # for additional particles
    Ey::Float64
    Ez::Float64
    w_total::Float64
    minvx::Float64
    maxvx::Float64
    minvy::Float64
    maxvy::Float64
    minvz::Float64
    maxvz::Float64

    n_moments::Int32
    rhs_vector::Vector{Float64}
    residual::Vector{Float64}
    mim::Vector{Vector{Int32}}  # mult-index moments
    tot_order::Vector{Int32}
    work::NNLSWorkspace
end

function create_nnls_merging(multi_index_moments, init_np)
    base_moments = base_multi_index_moments()
    filtered_mim = [x for x in multi_index_moments if !(x in base_moments)]
    n_total_conserved = length(base_moments) + length(filtered_mim)
    append!(base_moments, filtered_mim)
    tot_order = [sum(mim) for mim in base_moments]

    return NNLSMerging(SVector{3,Float64}(0.0, 0.0, 0.0), 1.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                       n_total_conserved,
                       zeros(n_total_conserved), zeros(n_total_conserved),
                       base_moments, tot_order,
                       NNLSWorkspace(zeros(n_total_conserved, init_np), zeros(n_total_conserved)))
end

function create_nnls_merging_rate_conserving(multi_index_moments, init_np)
    base_moments = base_multi_index_moments()
    filtered_mim = [x for x in multi_index_moments if !(x in base_moments)]
    n_total_conserved = length(base_moments) + length(filtered_mim) + 2
    append!(base_moments, filtered_mim)
    tot_order = [sum(mim) for mim in base_moments]

    return NNLSMerging(SVector{3,Float64}(0.0, 0.0, 0.0), 1.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                       0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                       n_total_conserved - 2,  # we conserve N moments and 2 rates, so n_moments = n_total_conserved - 2
                       zeros(n_total_conserved), zeros(n_total_conserved),
                       base_moments, tot_order,
                       NNLSWorkspace(zeros(n_total_conserved, init_np), zeros(n_total_conserved)))
end

function vx_sign(octant)
    if octant % 2 == 1
        return -1
    else
        return 1
    end
end

function vy_sign(octant)
    if (octant == 3) || (octant == 4) || (octant == 7) || (octant == 8)
        return 1
    else
        return -1
    end
end

function vz_sign(octant)
    if octant >= 5
        return 1
    else
        return -1
    end
end

function check_speed_bound(speed_val, minval, maxval)
    return max(min(speed_val, maxval), minval)
end

function check_speed_bound(speed_val, minval, maxval, multiplier)
    ms = multiplier * speed_val
    if ms > maxval
        return maxval * multiplier
    elseif ms < minval
        return minval * multiplier
    else
        return ms
    end
end

function base_multi_index_moments()
    # mass/momentum/directional energy conservation
    return [[0, 0, 0],
            [1, 0, 0], [0, 1, 0], [0, 0, 1],
            [2, 0, 0], [0, 2, 0], [0, 0, 2]]
end

function compute_multi_index_moments(n)
    result = []
    for i in 0:n
        for j in 0:n-i
            k = n - i - j
            if k >= 0
                push!(result, [i, j, k])
            end
        end
    end
    return result
end

function compute_w_total_v0!(nnls_merging, particles, pia, cell, species)
    nnls_merging.v0 = SVector{3,Float64}(0.0, 0.0, 0.0)
    nnls_merging.w_total = 0.0

    nnls_merging.minvx = c_light
    nnls_merging.maxvx = -c_light

    nnls_merging.minvy = c_light
    nnls_merging.maxvy = -c_light

    nnls_merging.minvz = c_light
    nnls_merging.maxvz = -c_light

    for i in pia.indexer[cell,species].start1:pia.indexer[cell,species].end1
        nnls_merging.w_total += particles[i].w
        nnls_merging.v0 = nnls_merging.v0 + particles[i].w * particles[i].v

        if (particles[i].v[1] > nnls_merging.maxvx)
            nnls_merging.maxvx = particles[i].v[1]
        end
        if (particles[i].v[1] < nnls_merging.minvx)
            nnls_merging.minvx = particles[i].v[1]
        end

        if (particles[i].v[2] > nnls_merging.maxvy)
            nnls_merging.maxvy = particles[i].v[2]
        end
        if (particles[i].v[2] < nnls_merging.minvy)
            nnls_merging.minvy = particles[i].v[2]
        end

        if (particles[i].v[3] > nnls_merging.maxvz)
            nnls_merging.maxvz = particles[i].v[3]
        end
        if (particles[i].v[3] < nnls_merging.minvz)
            nnls_merging.minvz = particles[i].v[3]
        end
    end

    if pia.indexer[cell,species].start2 > 0
        for i in pia.indexer[cell,species].start2:pia.indexer[cell,species].end2
            nnls_merging.w_total += particles[i].w
            nnls_merging.v0 = nnls_merging.v0 + particles[i].w * particles[i].v

            if (particles[i].v[1] > nnls_merging.maxvx)
                nnls_merging.maxvx = particles[i].v[1]
            end
            if (particles[i].v[1] < nnls_merging.minvx)
                nnls_merging.minvx = particles[i].v[1]
            end

            if (particles[i].v[2] > nnls_merging.maxvy)
                nnls_merging.maxvy = particles[i].v[2]
            end
            if (particles[i].v[2] < nnls_merging.minvy)
                nnls_merging.minvy = particles[i].v[2]
            end

            if (particles[i].v[3] > nnls_merging.maxvz)
                nnls_merging.maxvz = particles[i].v[3]
            end
            if (particles[i].v[3] < nnls_merging.minvz)
                nnls_merging.minvz = particles[i].v[3]
            end
        end
    end

    nnls_merging.v0 = nnls_merging.v0 / nnls_merging.w_total
    nnls_merging.minvx -= nnls_merging.v0[1]
    nnls_merging.minvy -= nnls_merging.v0[2]
    nnls_merging.minvz -= nnls_merging.v0[3]

    nnls_merging.maxvx -= nnls_merging.v0[1]
    nnls_merging.maxvy -= nnls_merging.v0[2]
    nnls_merging.maxvz -= nnls_merging.v0[3]
end


function ccm(v, v0, mim)
    # computed unweighted centered moment
    return (v[1] - v0[1])^mim[1] * (v[2] - v0[2])^mim[2] * (v[3] - v0[3])^mim[3]
end

function ccm(vx, vy, vz, vx0, vy0, vz0, mim)
    # computed unweighted centered moment
    return (vx - vx0)^mim[1] * (vy - vy0)^mim[2] * (vz - vz0)^mim[3]
end

function compute_lhs_and_rhs!(nnls_merging, lhs_matrix,
                              particles, pia, cell, species)
    nnls_merging.Ex = 0.0
    nnls_merging.Ey = 0.0
    nnls_merging.Ez = 0.0
    nnls_merging.rhs_vector .= 0.0

    compute_w_total_v0!(nnls_merging, particles, pia, cell, species)

    col_index = 1
    for i in pia.indexer[cell,species].start1:pia.indexer[cell,species].end1
        # w_total += particles[i].w

        nnls_merging.Ex += particles[i].w * (particles[i].v[1] - nnls_merging.v0[1])^2
        nnls_merging.Ey += particles[i].w * (particles[i].v[2] - nnls_merging.v0[2])^2
        nnls_merging.Ez += particles[i].w * (particles[i].v[3] - nnls_merging.v0[3])^2

        for n_mom in 1:nnls_merging.n_moments
            tmp_ccm = ccm(particles[i].v, nnls_merging.v0, nnls_merging.mim[n_mom])
            nnls_merging.rhs_vector[n_mom] += particles[i].w * tmp_ccm
            lhs_matrix[n_mom, col_index] = tmp_ccm
        end
        col_index += 1
    end

    if pia.indexer[cell,species].start2 > 0
        for i in pia.indexer[cell,species].start2:pia.indexer[cell,species].end2

            nnls_merging.Ex += particles[i].w * (particles[i].v[1] - nnls_merging.v0[1])^2
            nnls_merging.Ey += particles[i].w * (particles[i].v[2] - nnls_merging.v0[2])^2
            nnls_merging.Ez += particles[i].w * (particles[i].v[3] - nnls_merging.v0[3])^2
            # w_total += particles[i].w
            for n_mom in 1:nnls_merging.n_moments
                tmp_ccm = ccm(particles[i].v, nnls_merging.v0, nnls_merging.mim[n_mom])
                nnls_merging.rhs_vector[n_mom] += particles[i].w * tmp_ccm
                lhs_matrix[n_mom, col_index] = tmp_ccm
            end
            col_index += 1
        end
    end

    nnls_merging.rhs_vector /= nnls_merging.w_total

    nnls_merging.Ex /= nnls_merging.w_total
    nnls_merging.Ey /= nnls_merging.w_total
    nnls_merging.Ez /= nnls_merging.w_total
    nnls_merging.Ex = sqrt(nnls_merging.Ex)
    nnls_merging.Ey = sqrt(nnls_merging.Ey)
    nnls_merging.Ez = sqrt(nnls_merging.Ez)

    return col_index
end


function compute_lhs_and_rhs_rate_preserving!(nnls_merging, lhs_matrix,
                                              interaction, electron_neutral_interactions, computed_cs,
                                              particles, pia, cell, species, neutral_species_index)
    nnls_merging.Ex = 0.0
    nnls_merging.Ey = 0.0
    nnls_merging.Ez = 0.0
    nnls_merging.rhs_vector .= 0.0

    compute_w_total_v0!(nnls_merging, particles, pia, cell, species)
    # v0 = norm(nnls_merging.v0)

    col_index = 1
    for i in pia.indexer[cell,species].start1:pia.indexer[cell,species].end1
        # w_total += particles[i].w

        nnls_merging.Ex += particles[i].w * (particles[i].v[1] - nnls_merging.v0[1])^2
        nnls_merging.Ey += particles[i].w * (particles[i].v[2] - nnls_merging.v0[2])^2
        nnls_merging.Ez += particles[i].w * (particles[i].v[3] - nnls_merging.v0[3])^2

        for n_mom in 1:nnls_merging.n_moments
            tmp_ccm = ccm(particles[i].v, nnls_merging.v0, nnls_merging.mim[n_mom])
            nnls_merging.rhs_vector[n_mom] += particles[i].w * tmp_ccm
            lhs_matrix[n_mom, col_index] = tmp_ccm
        end

        g = norm(particles[i].v)
        compute_cross_sections_only!(computed_cs, interaction, g, electron_neutral_interactions, neutral_species_index)
        lhs_matrix[nnls_merging.n_moments+1, col_index] = get_cs_elastic(electron_neutral_interactions, computed_cs, neutral_species_index) * g
        lhs_matrix[nnls_merging.n_moments+2, col_index] = get_cs_ionization(electron_neutral_interactions, computed_cs, neutral_species_index) * g
        nnls_merging.rhs_vector[nnls_merging.n_moments+1] += get_cs_elastic(electron_neutral_interactions, computed_cs, neutral_species_index) * g * particles[i].w
        nnls_merging.rhs_vector[nnls_merging.n_moments+2] += get_cs_ionization(electron_neutral_interactions, computed_cs, neutral_species_index) * g * particles[i].w

        col_index += 1
    end

    if pia.indexer[cell,species].start2 > 0
        for i in pia.indexer[cell,species].start2:pia.indexer[cell,species].end2

            nnls_merging.Ex += particles[i].w * (particles[i].v[1] - nnls_merging.v0[1])^2
            nnls_merging.Ey += particles[i].w * (particles[i].v[2] - nnls_merging.v0[2])^2
            nnls_merging.Ez += particles[i].w * (particles[i].v[3] - nnls_merging.v0[3])^2
            # w_total += particles[i].w
            for n_mom in 1:nnls_merging.n_moments
                tmp_ccm = ccm(particles[i].v, nnls_merging.v0, nnls_merging.mim[n_mom])
                nnls_merging.rhs_vector[n_mom] += particles[i].w * tmp_ccm
                lhs_matrix[n_mom, col_index] = tmp_ccm
            end

            g = norm(particles[i].v)
            compute_cross_sections_only!(computed_cs, interaction, g, electron_neutral_interactions, neutral_species_index)
            lhs_matrix[nnls_merging.n_moments+1, col_index] = get_cs_elastic(electron_neutral_interactions, computed_cs, neutral_species_index) * g
            lhs_matrix[nnls_merging.n_moments+2, col_index] = get_cs_ionization(electron_neutral_interactions, computed_cs, neutral_species_index) * g
            nnls_merging.rhs_vector[nnls_merging.n_moments+1] += get_cs_elastic(electron_neutral_interactions, computed_cs, neutral_species_index) * g * particles[i].w
            nnls_merging.rhs_vector[nnls_merging.n_moments+2] += get_cs_ionization(electron_neutral_interactions, computed_cs, neutral_species_index) * g * particles[i].w

            col_index += 1
        end
    end

    nnls_merging.rhs_vector /= nnls_merging.w_total

    nnls_merging.Ex /= nnls_merging.w_total
    nnls_merging.Ey /= nnls_merging.w_total
    nnls_merging.Ez /= nnls_merging.w_total
    nnls_merging.Ex = sqrt(nnls_merging.Ex)
    nnls_merging.Ey = sqrt(nnls_merging.Ey)
    nnls_merging.Ez = sqrt(nnls_merging.Ez)

    return col_index
end


function compute_lhs_particles_additional!(rng, col_index, nnls_merging, lhs_matrix,
                                           particles, pia, cell, species,
                                           n_rand_pairs)
    # a particle centered at local zero
    for n_mom in 1:nnls_merging.n_moments
        lhs_matrix[n_mom, col_index] = 0.0
    end
    col_index += 1

    # add 8 additional particles
    for i in 1:8
        vxs = vx_sign(i)
        vys = vy_sign(i)
        vzs = vz_sign(i)

        vx = check_speed_bound(vxs * nnls_merging.Ex, nnls_merging.minvx, nnls_merging.maxvx, 0.5)
        vy = check_speed_bound(vys * nnls_merging.Ey, nnls_merging.minvy, nnls_merging.maxvy, 0.5)
        vz = check_speed_bound(vzs * nnls_merging.Ez, nnls_merging.minvz, nnls_merging.maxvz, 0.5)
        for n_mom in 1:nnls_merging.n_moments
            lhs_matrix[n_mom, col_index] = ccm(vx, vy, vz,
                                               0.0, 0.0, 0.0,
                                               nnls_merging.mim[n_mom])
        end
        col_index += 1
    end

    # add 8 additional particles
    for i in 1:8
        vxs = vx_sign(i)
        vys = vy_sign(i)
        vzs = vz_sign(i)
        
        vx = check_speed_bound(vxs * nnls_merging.Ex, nnls_merging.minvx, nnls_merging.maxvx, 1.0)
        vy = check_speed_bound(vys * nnls_merging.Ey, nnls_merging.minvy, nnls_merging.maxvy, 1.0)
        vz = check_speed_bound(vzs * nnls_merging.Ez, nnls_merging.minvz, nnls_merging.maxvz, 1.0)
        for n_mom in 1:nnls_merging.n_moments
            lhs_matrix[n_mom, col_index] = ccm(vx, vy, vz,
                                               0.0, 0.0, 0.0,
                                               nnls_merging.mim[n_mom])
        end
        col_index += 1
    end

    for _ in 1:n_rand_pairs
        i = floor(Int64, rand(rng, Float64) * pia.indexer[cell, species].n_local)
        k = floor(Int64, rand(rng, Float64) * pia.indexer[cell, species].n_local)

        i = map_cont_index(pia.indexer[cell, species], i)
        k = map_cont_index(pia.indexer[cell, species], k)

        if (i != k)
            vx = 0.5 * (particles[i].v[1] + particles[k].v[1])
            vy = 0.5 * (particles[i].v[2] + particles[k].v[2])
            vz = 0.5 * (particles[i].v[3] + particles[k].v[3])
            for n_mom in 1:nnls_merging.n_moments
                lhs_matrix[n_mom, col_index] = ccm(vx, vy, vz,
                                                   nnls_merging.v0[1], nnls_merging.v0[2], nnls_merging.v0[3],
                                                   nnls_merging.mim[n_mom])
            end
        end
        col_index += 1
    end
end


function compute_lhs_particles_additional_rate_preserving!(rng, col_index, nnls_merging, lhs_matrix,
                                           interaction, electron_neutral_interactions, computed_cs, 
                                           particles, pia, cell, species, neutral_species_index,
                                           n_rand_pairs)

    v0 = norm(nnls_merging.v0)
    # a particle centered at local zero
    for n_mom in 1:nnls_merging.n_moments
        lhs_matrix[n_mom, col_index] = 0.0
    end

    compute_cross_sections_only!(computed_cs, interaction, v0, electron_neutral_interactions, neutral_species_index)
    lhs_matrix[nnls_merging.n_moments+1, col_index] = get_cs_elastic(electron_neutral_interactions, computed_cs, neutral_species_index) * v0
    lhs_matrix[nnls_merging.n_moments+2, col_index] = get_cs_ionization(electron_neutral_interactions, computed_cs, neutral_species_index) * v0

    col_index += 1

    # add 8 additional particles
    for i in 1:8
        vxs = vx_sign(i)
        vys = vy_sign(i)
        vzs = vz_sign(i)

        vx = check_speed_bound(vxs * nnls_merging.Ex, nnls_merging.minvx, nnls_merging.maxvx, 0.5)
        vy = check_speed_bound(vys * nnls_merging.Ey, nnls_merging.minvy, nnls_merging.maxvy, 0.5)
        vz = check_speed_bound(vzs * nnls_merging.Ez, nnls_merging.minvz, nnls_merging.maxvz, 0.5)
        for n_mom in 1:nnls_merging.n_moments
            lhs_matrix[n_mom, col_index] = ccm(vx, vy, vz,
                                               0.0, 0.0, 0.0,
                                               nnls_merging.mim[n_mom])
        end

        g = sqrt((nnls_merging.v0[1] + vx)^2 + (nnls_merging.v0[2] + vy)^2 + (nnls_merging.v0[3] + vz)^2)
        compute_cross_sections_only!(computed_cs, interaction, g, electron_neutral_interactions, neutral_species_index)
        lhs_matrix[nnls_merging.n_moments+1, col_index] = get_cs_elastic(electron_neutral_interactions, computed_cs, neutral_species_index) * g
        lhs_matrix[nnls_merging.n_moments+2, col_index] = get_cs_ionization(electron_neutral_interactions, computed_cs, neutral_species_index) * g

        col_index += 1
    end

    # add 8 additional particles
    for i in 1:8
        vxs = vx_sign(i)
        vys = vy_sign(i)
        vzs = vz_sign(i)

        vx = check_speed_bound(vxs * nnls_merging.Ex, nnls_merging.minvx, nnls_merging.maxvx, 1.0)
        vy = check_speed_bound(vys * nnls_merging.Ey, nnls_merging.minvy, nnls_merging.maxvy, 1.0)
        vz = check_speed_bound(vzs * nnls_merging.Ez, nnls_merging.minvz, nnls_merging.maxvz, 1.0)
        for n_mom in 1:nnls_merging.n_moments
            lhs_matrix[n_mom, col_index] = ccm(vx, vy, vz,
                                               0.0, 0.0, 0.0,
                                               nnls_merging.mim[n_mom])
        end

        g = sqrt((nnls_merging.v0[1] + vx)^2 + (nnls_merging.v0[2] + vy)^2 + (nnls_merging.v0[3] + vz)^2)
        compute_cross_sections_only!(computed_cs, interaction, g, electron_neutral_interactions, neutral_species_index)
        lhs_matrix[nnls_merging.n_moments+1, col_index] = get_cs_elastic(electron_neutral_interactions, computed_cs, neutral_species_index) * g
        lhs_matrix[nnls_merging.n_moments+2, col_index] = get_cs_ionization(electron_neutral_interactions, computed_cs, neutral_species_index) * g
        
        col_index += 1
    end

    for _ in 1:n_rand_pairs
        i = floor(Int64, rand(rng, Float64) * pia.indexer[cell, species].n_local)
        k = floor(Int64, rand(rng, Float64) * pia.indexer[cell, species].n_local)

        i = map_cont_index(pia.indexer[cell, species], i)
        k = map_cont_index(pia.indexer[cell, species], k)

        if (i != k)
            vx = 0.5 * (particles[i].v[1] + particles[k].v[1])
            vy = 0.5 * (particles[i].v[2] + particles[k].v[2])
            vz = 0.5 * (particles[i].v[3] + particles[k].v[3])
            for n_mom in 1:nnls_merging.n_moments
                lhs_matrix[n_mom, col_index] = ccm(vx, vy, vz,
                                                   nnls_merging.v0[1], nnls_merging.v0[2], nnls_merging.v0[3],
                                                   nnls_merging.mim[n_mom])
            end

            g = sqrt((nnls_merging.v0[1] + vx)^2 + (nnls_merging.v0[2] + vy)^2 + (nnls_merging.v0[3] + vz)^2)
            compute_cross_sections_only!(computed_cs, interaction, g, electron_neutral_interactions, neutral_species_index)
            lhs_matrix[nnls_merging.n_moments+1, col_index] = get_cs_elastic(electron_neutral_interactions, computed_cs, neutral_species_index) * g
            lhs_matrix[nnls_merging.n_moments+2, col_index] = get_cs_ionization(electron_neutral_interactions, computed_cs, neutral_species_index) * g
            
        end
        col_index += 1
    end
end


function scale_lhs_rhs!(nnls_merging, lhs_matrix)
    for n_mom in 1:nnls_merging.n_moments
        lhs_matrix[n_mom, :] .*= nnls_merging.inv_vref^nnls_merging.tot_order[n_mom]
        nnls_merging.rhs_vector[n_mom] *= nnls_merging.inv_vref^nnls_merging.tot_order[n_mom]
    end
end

function scale_lhs_rhs_rate_preserving!(nnls_merging, lhs_matrix, ref_cs_elastic, ref_cs_ion)
    for n_mom in 1:nnls_merging.n_moments
        lhs_matrix[n_mom, :] .*= nnls_merging.inv_vref^nnls_merging.tot_order[n_mom]
        nnls_merging.rhs_vector[n_mom] *= nnls_merging.inv_vref^nnls_merging.tot_order[n_mom]
    end
    lhs_matrix[nnls_merging.n_moments+1, :] .*= nnls_merging.inv_vref / ref_cs_elastic
    lhs_matrix[nnls_merging.n_moments+2, :] .*= nnls_merging.inv_vref / ref_cs_ion

    nnls_merging.rhs_vector[nnls_merging.n_moments+1] *= nnls_merging.inv_vref / ref_cs_elastic
    nnls_merging.rhs_vector[nnls_merging.n_moments+2] *= nnls_merging.inv_vref / ref_cs_ion
end


function compute_post_merge_particles_nnls!(lhs_ncols, lhs_matrix,
                                            nnls_merging, particles, pia, cell, species)

    curr_particle_index = 0

    # lhs matrix has row 1 of ones (mass conservation)
    # row 2 are the vx components (Vx conservation)
    # row 3 are the vy components (Vy conservation)
    # row 4 are the vz components (Vz conservation)
    for j in 1:lhs_ncols
        if nnls_merging.work.x[j] > 0.0
            i = map_cont_index(pia.indexer[cell,species], curr_particle_index)
            curr_particle_index += 1
            particles[i].w = nnls_merging.work.x[j] * nnls_merging.w_total
            particles[i].v = nnls_merging.v0 + SVector{3, Float64}(lhs_matrix[2, j],
                                                                   lhs_matrix[3, j],
                                                                   lhs_matrix[4, j]) * nnls_merging.vref
        end
    end

    update_particle_indexer_new_lower_count(pia, cell, species, curr_particle_index)
end

function merge_nnls_based!(rng, nnls_merging, particles, pia, cell, species, vref; n_rand_pairs=0, max_err=1e-11)

    # create LHS matrix
    lhs_ncols = pia.indexer[cell, species].n_local + 8 + 1 + 8 + n_rand_pairs
    lhs_matrix = zeros(nnls_merging.n_moments, lhs_ncols)
    nnls_merging.vref = vref
    nnls_merging.inv_vref = 1.0 / vref

    # create LHS matrix and fill RHS vector using existing particles
    col_index = compute_lhs_and_rhs!(nnls_merging, lhs_matrix, particles, pia, cell, species)
    # and add more columns
    compute_lhs_particles_additional!(rng, col_index, nnls_merging, lhs_matrix,
                                      particles, pia, cell, species, n_rand_pairs)
    scale_lhs_rhs!(nnls_merging, lhs_matrix)

    load!(nnls_merging.work, lhs_matrix, nnls_merging.rhs_vector)
    solve!(nnls_merging.work, 2 * size(lhs_matrix, 2))

    nnls_merging.residual .= lhs_matrix * nnls_merging.work.x

    error = 0.0
    for i in 1:nnls_merging.n_moments
        error += abs(nnls_merging.residual[i] -  nnls_merging.rhs_vector[i])
    end

    nonzero = sum(nnls_merging.work.x .> 0.0)

    if error > max_err
        # println("Large L1 Error: ", error)
        return -1
    end

    if nonzero == pia.indexer[cell, species].n_local
        # println("Non-sparse solution ", error)
        return -1
    end

    nnls_merging.work.x =  nnls_merging.work.x / sum(nnls_merging.work.x)
    
    compute_post_merge_particles_nnls!(lhs_ncols, lhs_matrix,
                                       nnls_merging, particles, pia, cell, species)
    return 1
end

function merge_nnls_based_rate_preserving!(rng, nnls_merging,
                                           interaction, electron_neutral_interactions, computed_cs,
                                           particles, pia, cell, species, neutral_species_index,
                                           vref, ref_cs_elatic, ref_cs_ion; n_rand_pairs=0, max_err=1e-11)

    # create LHS matrix
    lhs_ncols = pia.indexer[cell, species].n_local + 8 + 1 + 8 + n_rand_pairs
    lhs_matrix = zeros(nnls_merging.n_moments+2, lhs_ncols)
    nnls_merging.vref = vref
    nnls_merging.inv_vref = 1.0 / vref

    # create LHS matrix and fill RHS vector using existing particles
    col_index = compute_lhs_and_rhs_rate_preserving!(nnls_merging, lhs_matrix,
                                                     interaction[species,neutral_species_index], electron_neutral_interactions, computed_cs, 
                                                     particles, pia, cell, species, neutral_species_index)
    # and add more columns
    compute_lhs_particles_additional_rate_preserving!(rng, col_index, nnls_merging, lhs_matrix,
                                      interaction[species,neutral_species_index], electron_neutral_interactions, computed_cs,
                                      particles, pia, cell, species, neutral_species_index, n_rand_pairs)
    scale_lhs_rhs_rate_preserving!(nnls_merging, lhs_matrix, ref_cs_elatic, ref_cs_ion)

    load!(nnls_merging.work, lhs_matrix, nnls_merging.rhs_vector)
    solve!(nnls_merging.work, 2 * size(lhs_matrix, 2))

    nnls_merging.residual .= lhs_matrix * nnls_merging.work.x

    error = 0.0
    for i in 1:nnls_merging.n_moments
        error += abs(nnls_merging.residual[i] -  nnls_merging.rhs_vector[i])
    end

    nonzero = sum(nnls_merging.work.x .> 0.0)

    if error > max_err
        # println("Large L1 Error: ", error)
        return -1
    end

    if nonzero == pia.indexer[cell, species].n_local
        # println("Non-sparse solution ", error)
        return -1
    end

    nnls_merging.work.x =  nnls_merging.work.x / sum(nnls_merging.work.x)
    
    compute_post_merge_particles_nnls!(lhs_ncols, lhs_matrix,
                                       nnls_merging, particles, pia, cell, species)
    return 1
end


# TODO: ResizableArray versions
# TODO: spatial moments!

# TODO: preallocate particle weight vector / vx, vy, vz vector, octree splitting, write out particles