@muladd begin

"""
    NNLSMerge
Struct for keeping track of merging-related quantities for NNLS-based merging.

# Fields
* `v0`: vector of the mean velocity of the particles
* `vref`: reference velocity magnitude for scaling
* `inv_vref`: inverse of reference velocity magnitude for scaling
* `Ex`: standard deviation of x velocity of particles
* `Ey`: standard deviation of y velocity of particles
* `Ez`: standard deviation of z velocity of particles
* `w_total`: total computational of the particles
* `minvx`: minimum x velocity of the particles
* `maxvx`: maximum x velocity of the particles
* `minvy`: minimum y velocity of the particles
* `maxvy`: maximum y velocity of the particles
* `minvz`: minimum z velocity of the particles
* `maxvz`: maximum z velocity of the particles
* `n_moments`: number of velocity moments to preserve
* `rhs_vector`: vector of computed moments
* `residual`: residual of solution
* `mim`: vector of 3-tuples of multi-indices for the velocity moments to preserve
* `tot_order`: vector of total orders of the velocity moments to preserve
* `work`: `NNLSWorkspace` instance
"""
mutable struct NNLSMerge
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

    @doc """
        NNLSMerge(multi_index_moments, init_np; rate_preserving=false)

    Create NNLS-based merging. Mass, momentum, directional energy are always conserved:
    if not in the list `multi_index_moments` of moments to preserve, the corresponding
    moment multi-indices will be added automatically. These indices are `(0,0,0)` for mass,
    `(1,0,0)`, `(0,1,0)`, `(0,0,1)` for momentum, and `(2,0,0)`, `(0,2,0)`, `(0,0,2)`
    for directional energies.

    # Positional arguments
    * `multi_index_moments`: vector of mixed moments to preserve of the form `[(i1, j1, k1), (i2, j2, k2), ...]``
    * `init_np`: assumption on pre-merge number of particles to pre-allocate memory for
    
    Keyword arguments:
    * `rate_preserving`: used for rate-preserving merging of electrons, preserves approximate elastic collision and ionization rates
    """
    function NNLSMerge(multi_index_moments, init_np; rate_preserving=false)
        add_length = 0
        if rate_preserving
            add_length = 2
        end

        base_moments = base_multi_index_moments()
        filtered_mim = [x for x in multi_index_moments if !(x in base_moments)]
        n_total_conserved = length(base_moments) + length(filtered_mim) + add_length
        append!(base_moments, filtered_mim)
        tot_order = [sum(mim) for mim in base_moments]
    
        return new(SVector{3,Float64}(0.0, 0.0, 0.0), 1.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   n_total_conserved - add_length,
                   zeros(n_total_conserved), zeros(n_total_conserved),
                   base_moments, tot_order,
                   NNLSWorkspace(zeros(n_total_conserved, init_np), zeros(n_total_conserved)))
    end
end

"""
    vx_sign(octant)
    
Return sign of velocity vx of an octant in velocity space.

# Positional arguments
* `octant`: index of the octant in velocity space

# Returns
Sign of the x-velocity corresponding to the octant.
"""
function vx_sign(octant)
    if octant % 2 == 1
        return -1
    else
        return 1
    end
end

"""
    vy_sign(octant)
    
Return sign of velocity vy of an octant in velocity space.

# Positional arguments
* `octant`: index of the octant in velocity space

# Returns
Sign of the y-velocity corresponding to the octant.
"""
function vy_sign(octant)
    if (octant == 3) || (octant == 4) || (octant == 7) || (octant == 8)
        return 1
    else
        return -1
    end
end

"""
    vz_sign(octant)

Return sign of velocity `vz` of an octant in velocity space.

# Positional arguments
* `octant`: index of the octant in velocity space

# Returns
Sign of the z-velocity corresponding to the octant
"""
function vz_sign(octant)
    if octant >= 5
        return 1
    else
        return -1
    end
end

"""
    check_speed_bound(speed_val, minval, maxval)
    
Compare `speed_val` to the bounds and return either the value if its within
bounds or the appropriate closest bound.

# Positional arguments
* `speed_val`: value to check
* `minval`: lower bound
* `maxval`: upper bound

# Returns
Either original value or value of closest bound.
"""
function check_speed_bound(speed_val, minval, maxval)
    return max(min(speed_val, maxval), minval)
end

"""
    check_speed_bound(speed_val, minval, maxval, multiplier)

Compare `speed_val` time a multiplier to the bounds and return either the value times
the multiplier if its within
bounds or the appropriate closest bound times the multiplier.

# Positional arguments
* `speed_val`: value to check
* `minval`: lower bound
* `maxval`: upper bound
* `multiplier`: the multiplier

# Returns
Either original value times the multiplier or value of closest bound times the multiplier.
"""
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

"""
    base_multi_index_moments()

Base multi indices corresponding to conservation of mass, momentum, and energy components.

# Positional arguments
None

# Returns
A vector of multi-indices (3-tuples) which corresponding to mass, momentum, and directional energy.
"""
function base_multi_index_moments()
    # mass/momentum/directional energy conservation
    return [[0, 0, 0],
            [1, 0, 0], [0, 1, 0], [0, 0, 1],
            [2, 0, 0], [0, 2, 0], [0, 0, 2]]
end

"""
    compute_multi_index_moments(n)

Compute all mixed moment multi-indices of total order up to n, i.e. all 3-tuples `(i,j,k)`` such that
`i+j+k <= n`.

# Positional arguments
* `n`: maximum total order

# Returns
Vector of 3-tuples of moment multi-indices.
"""
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

"""
    compute_w_total_v0!(nnls_merging, particles, pia, cell, species)

Compute total computational weight of particles and mean velocity, as well as velocity bounds of the
set of particles in each velocity direction.

# Positional arguments
* `nnls_merging`: the `NNLSMerge` instance where the computed values will be stored
* `particles`: the `ParticleVector` instance containing the particles that are being merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the grid cell in which particles are being merged
* `species`: the index of the species being merged
"""
function compute_w_total_v0!(nnls_merging, particles, pia, cell, species)
    nnls_merging.v0 = SVector{3,Float64}(0.0, 0.0, 0.0)
    nnls_merging.w_total = 0.0

    nnls_merging.minvx = c_light
    nnls_merging.maxvx = -c_light

    nnls_merging.minvy = c_light
    nnls_merging.maxvy = -c_light

    nnls_merging.minvz = c_light
    nnls_merging.maxvz = -c_light

    @inbounds s1 = pia.indexer[cell,species].start1
    @inbounds e1 = pia.indexer[cell,species].end1
    @inbounds for i in s1:e1
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

    @inbounds if pia.indexer[cell,species].n_group2 > 0
        @inbounds s2 = pia.indexer[cell,species].start2
        @inbounds e2 = pia.indexer[cell,species].end2
        @inbounds for i in s2:e2
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
    @inbounds nnls_merging.minvx -= nnls_merging.v0[1]
    @inbounds nnls_merging.minvy -= nnls_merging.v0[2]
    @inbounds nnls_merging.minvz -= nnls_merging.v0[3]

    @inbounds nnls_merging.maxvx -= nnls_merging.v0[1]
    @inbounds nnls_merging.maxvy -= nnls_merging.v0[2]
    @inbounds nnls_merging.maxvz -= nnls_merging.v0[3]
end

"""
    ccm(v, v0, mim)

Compute unweighted central moment.

# Positional arguments
* `v`: the 3-dimensional velocity vector
* `v0`: the 3-dimensional mean velocity vector
* `mim`: the 3-dimensional multi-index

# Returns
Computed unweighted central moment.
"""
@inline function ccm(v, v0, mim)
    @inbounds return (v[1] - v0[1])^mim[1] * (v[2] - v0[2])^mim[2] * (v[3] - v0[3])^mim[3]
end

"""
    ccm(vx, vy, vz, vx0, vy0, vz0, mim)

Compute unweighted central moment.

# Positional arguments
* `vx`: the x component of the velocity
* `vy`: the y component of the velocity
* `vz`: the z component of the velocity
* `vx0`: the x component of the mean velocity
* `vy0`: the y component of the mean velocity
* `vz0`: the z component of the mean velocity
* `mim`: the 3-dimensional multi-index

# Returns
Computed unweighted central moment.
"""
@inline function ccm(vx, vy, vz, vx0, vy0, vz0, mim)
    @inbounds return (vx - vx0)^mim[1] * (vy - vy0)^mim[2] * (vz - vz0)^mim[3]
end

"""
    compute_lhs_and_rhs!(nnls_merging, lhs_matrix, particles, pia, cell, species)

Compute LHS matrix and RHS vector for NNLS merging. Returns the pre-merge number of particles.

# Positional arguments
* `nnls_merging`: the `NNLSMerge` instance where the RHS vector will be stored
* `lhs_matrix`: the matrix of size `n_total_conserved x n_particles`, where
    `n_total_conserved` is the number of conserved moments and `n_particles` is the pre-merge number of particles
* `particles`: the `ParticleVector` instance containing the particles that are being merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the grid cell in which particles are being merged
* `species`: the index of the species being merged

# Returns
The pre-merge number of particles.
"""
function compute_lhs_and_rhs!(nnls_merging, lhs_matrix,
                              particles, pia, cell, species)
    n_moms = nnls_merging.n_moments
    
    nnls_merging.Ex = 0.0
    nnls_merging.Ey = 0.0
    nnls_merging.Ez = 0.0
    nnls_merging.rhs_vector .= 0.0

    compute_w_total_v0!(nnls_merging, particles, pia, cell, species)

    col_index = 1
    @inbounds s1 = pia.indexer[cell,species].start1
    @inbounds e1 = pia.indexer[cell,species].end1

    @inbounds for i in s1:e1
        # w_total += particles[i].w

        nnls_merging.Ex = nnls_merging.Ex + particles[i].w * (particles[i].v[1] - nnls_merging.v0[1])^2
        nnls_merging.Ey = nnls_merging.Ey + particles[i].w * (particles[i].v[2] - nnls_merging.v0[2])^2
        nnls_merging.Ez = nnls_merging.Ez + particles[i].w * (particles[i].v[3] - nnls_merging.v0[3])^2

        for n_mom in 1:n_moms
            tmp_ccm = ccm(particles[i].v, nnls_merging.v0, nnls_merging.mim[n_mom])
            nnls_merging.rhs_vector[n_mom] = nnls_merging.rhs_vector[n_mom] + particles[i].w * tmp_ccm
            lhs_matrix[n_mom, col_index] = tmp_ccm
        end
        col_index += 1
    end

    @inbounds if pia.indexer[cell,species].n_group2 > 0

        @inbounds s2 = pia.indexer[cell,species].start2
        @inbounds e2 = pia.indexer[cell,species].end2

        @inbounds for i in s2:e2
            nnls_merging.Ex = nnls_merging.Ex + particles[i].w * (particles[i].v[1] - nnls_merging.v0[1])^2
            nnls_merging.Ey = nnls_merging.Ey + particles[i].w * (particles[i].v[2] - nnls_merging.v0[2])^2
            nnls_merging.Ez = nnls_merging.Ez + particles[i].w * (particles[i].v[3] - nnls_merging.v0[3])^2
            # w_total += particles[i].w
            for n_mom in 1:nnls_merging.n_moments
                tmp_ccm = ccm(particles[i].v, nnls_merging.v0, nnls_merging.mim[n_mom])
                nnls_merging.rhs_vector[n_mom] = nnls_merging.rhs_vector[n_mom] + particles[i].w * tmp_ccm
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

"""
    compute_lhs_and_rhs_rate_preserving!(nnls_merging, lhs_matrix, particles, pia, cell, species)

Compute LHS matrix and RHS vector for the rate-preserving NNLS merging (for electrons). Approximate
    elastic scattering and electron-impact ionization rates are conserved.

# Positional arguments
* `nnls_merging`: the `NNLSMerge` instance where the RHS vector will be stored
* `lhs_matrix`: the matrix of size `n_total_conserved x n_particles`, where
    `n_total_conserved` is the number of conserved moments and `n_particles` is the pre-merge number of particles
* `interaction`: the `Interaction` instance describing the electron-neutral interaction being considered
* `electron_neutral_interactions`:  the `ElectronNeutralInteractions` instance storing the tabulated cross-section
    data used to compute the rates
* `computed_cs`: the vector of `ComputedCrossSection` instances in which the computed values will be stored
* `particles`: the `ParticleVector` instance containing the particles that are being merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the grid cell in which particles are being merged
* `species`: the index of the species being merged
* `neutral_species_index`: the index of the neutral species which is the collision partner in the electron-neutral
    collisions for which approximate rates are being preserved.
"""
function compute_lhs_and_rhs_rate_preserving!(nnls_merging, lhs_matrix,
                                              interaction, electron_neutral_interactions, computed_cs,
                                              particles, pia, cell, species, neutral_species_index)
    n_moms = nnls_merging.n_moments
    
    nnls_merging.Ex = 0.0
    nnls_merging.Ey = 0.0
    nnls_merging.Ez = 0.0
    nnls_merging.rhs_vector .= 0.0

    compute_w_total_v0!(nnls_merging, particles, pia, cell, species)
    # v0 = norm(nnls_merging.v0)

    col_index = 1
    @inbounds s1 = pia.indexer[cell,species].start1
    @inbounds e1 = pia.indexer[cell,species].end1
    @inbounds for i in s1:e1
        # w_total += particles[i].w

        nnls_merging.Ex = nnls_merging.Ex + particles[i].w * (particles[i].v[1] - nnls_merging.v0[1])^2
        nnls_merging.Ey = nnls_merging.Ey + particles[i].w * (particles[i].v[2] - nnls_merging.v0[2])^2
        nnls_merging.Ez = nnls_merging.Ez + particles[i].w * (particles[i].v[3] - nnls_merging.v0[3])^2

        for n_mom in 1:n_moms
            tmp_ccm = ccm(particles[i].v, nnls_merging.v0, nnls_merging.mim[n_mom])
            nnls_merging.rhs_vector[n_mom] = nnls_merging.rhs_vector[n_mom] + nnls_merging.rhs_vector[n_mom] + particles[i].w * tmp_ccm
            lhs_matrix[n_mom, col_index] = tmp_ccm
        end

        g = norm(particles[i].v)
        compute_cross_sections_only!(computed_cs, interaction, g, electron_neutral_interactions, neutral_species_index)
        lhs_matrix[nnls_merging.n_moments+1, col_index] = get_cs_elastic(electron_neutral_interactions, computed_cs, neutral_species_index) * g
        lhs_matrix[nnls_merging.n_moments+2, col_index] = get_cs_ionization(electron_neutral_interactions, computed_cs, neutral_species_index) * g
        nnls_merging.rhs_vector[nnls_merging.n_moments+1] = nnls_merging.rhs_vector[nnls_merging.n_moments+1] + get_cs_elastic(electron_neutral_interactions,
                                                                                                                               computed_cs, neutral_species_index) * g * particles[i].w
        nnls_merging.rhs_vector[nnls_merging.n_moments+2] = nnls_merging.rhs_vector[nnls_merging.n_moments+2] + get_cs_ionization(electron_neutral_interactions,
                                                                                                                                  computed_cs, neutral_species_index) * g * particles[i].w

        col_index += 1
    end

    @inbounds if pia.indexer[cell,species].n_group2 > 0
        @inbounds s2 = pia.indexer[cell,species].start2
        @inbounds e2 = pia.indexer[cell,species].end2

        @inbounds for i in s2:e2
            nnls_merging.Ex = nnls_merging.Ex + particles[i].w * (particles[i].v[1] - nnls_merging.v0[1])^2
            nnls_merging.Ey = nnls_merging.Ey + particles[i].w * (particles[i].v[2] - nnls_merging.v0[2])^2
            nnls_merging.Ez = nnls_merging.Ez + particles[i].w * (particles[i].v[3] - nnls_merging.v0[3])^2
            # w_total += particles[i].w
            for n_mom in 1:nnls_merging.n_moments
                tmp_ccm = ccm(particles[i].v, nnls_merging.v0, nnls_merging.mim[n_mom])
                nnls_merging.rhs_vector[n_mom] = nnls_merging.rhs_vector[n_mom] + particles[i].w * tmp_ccm
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

"""
    compute_lhs_particles_additional!(rng, col_index, nnls_merging, lhs_matrix,
                                      particles, pia, cell, species,
                                      n_rand_pairs, centered_at_mean, v_multipliers)

Compute additional LHS columns for additional particles. One particle is added at the local zero velocity (
i.e. mean velocity of the whole system of particles), 8 particles are added (1 per velocity space octant)
using the standard deviation of the velocities of the whole particle system in each direction to determine
their velocities (+ bounded via use of the [`check_speed_bound`](@ref) function). A second set of 8 particles
is added using a similar scheme, but the velocities are multiplied by 0.5.
Finally, `n_rand_pairs` particle pairs are chosen and their mean velocities are used to compute these
additional particles.

# Positional arguments
* `rng`: the random number generator instance
* `col_index`: the index of the first column where to start writing the additionally computed values
* `nnls_merging`: the `NNLSMerge` instance where the RHS vector will be stored
* `lhs_matrix`: the matrix of size `n_total_conserved x n_particles`, where
    `n_total_conserved` is the number of conserved moments and `n_particles` is the pre-merge number of particles
* `particles`: the `ParticleVector` instance containing the particles that are being merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the grid cell in which particles are being merged
* `species`: the index of the species being merged
* `n_rand_pairs`: the number of random particle pairs to choose
* `centered_at_mean`: whether to add a particle centered at the mean velocity of the system of particles
* `v_multipliers`: 
"""
function compute_lhs_particles_additional!(rng, col_index, nnls_merging, lhs_matrix,
                                           particles, pia, cell, species,
                                           n_rand_pairs, centered_at_mean, v_multipliers)

    if centered_at_mean
        n_moms = nnls_merging.n_moments
        # a particle centered at local zero
        @inbounds for n_mom in 1:n_moms
            lhs_matrix[n_mom, col_index] = 0.0
        end
        col_index += 1
    end

    # add 8 additional particles
    for v_mult in v_multipliers
        for i in 1:8
            vxs = vx_sign(i)
            vys = vy_sign(i)
            vzs = vz_sign(i)

            vx = check_speed_bound(vxs * nnls_merging.Ex, nnls_merging.minvx, nnls_merging.maxvx, v_mult)
            vy = check_speed_bound(vys * nnls_merging.Ey, nnls_merging.minvy, nnls_merging.maxvy, v_mult)
            vz = check_speed_bound(vzs * nnls_merging.Ez, nnls_merging.minvz, nnls_merging.maxvz, v_mult)
            @inbounds for n_mom in 1:nnls_merging.n_moments
                lhs_matrix[n_mom, col_index] = ccm(vx, vy, vz,
                                                             0.0, 0.0, 0.0,
                                                             nnls_merging.mim[n_mom])
            end
            col_index += 1
        end
    end

    for _ in 1:n_rand_pairs
        @inbounds i = floor(Int64, rand(rng, Float64) * pia.indexer[cell, species].n_local)
        @inbounds k = floor(Int64, rand(rng, Float64) * pia.indexer[cell, species].n_local)

        @inbounds i = map_cont_index(pia.indexer[cell, species], i)
        @inbounds k = map_cont_index(pia.indexer[cell, species], k)

        if (i != k)
            @inbounds vx = 0.5 * (particles[i].v[1] + particles[k].v[1])
            @inbounds vy = 0.5 * (particles[i].v[2] + particles[k].v[2])
            @inbounds vz = 0.5 * (particles[i].v[3] + particles[k].v[3])
            @inbounds for n_mom in 1:n_moms
                lhs_matrix[n_mom, col_index] = ccm(vx, vy, vz,
                                                   nnls_merging.v0[1], nnls_merging.v0[2], nnls_merging.v0[3],
                                                   nnls_merging.mim[n_mom])
            end
        end
        col_index += 1
    end
end

"""
    compute_lhs_particles_additional_rate_preserving!(rng, col_index, nnls_merging, lhs_matrix,
                                           interaction, electron_neutral_interactions, computed_cs, 
                                           particles, pia, cell, species, neutral_species_index,
                                           n_rand_pairs, centered_at_mean, v_multipliers)

Compute additional LHS columns for additional particles for the rate-preserving NNLS merging (for electrons).
One particle is added at the local zero velocity (i.e. mean velocity of the whole system of particles),
8 particles are added (1 per velocity space octant)
using the standard deviation of the velocities of the whole particle system in each direction to determine
their velocities (+ bounded via use of the [`check_speed_bound`](@ref) function). A second set of 8 particles
is added using a similar scheme, but the velocities are multiplied by 0.5.
Finally, `n_rand_pairs` particle pairs are chosen and their mean velocities are used to compute these
additional particles.

# Positional arguments
* `rng`: the random number generator instance
* `col_index`: the index of the first column where to start writing the additionally computed values
* `nnls_merging`: the `NNLSMerge` instance where the RHS vector will be stored
* `lhs_matrix`: the matrix of size `n_total_conserved x n_particles`, where
    `n_total_conserved` is the number of conserved moments and `n_particles` is the pre-merge number of particles
* `interaction`: the `Interaction` instance describing the electron-neutral interaction being considered
* `electron_neutral_interactions`:  the `ElectronNeutralInteractions` instance storing the tabulated cross-section
    data used to compute the rates
* `computed_cs`: the vector of `ComputedCrossSection` instances in which the computed values will be stored
* `particles`: the `ParticleVector` instance containing the particles that are being merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the grid cell in which particles are being merged
* `species`: the index of the species being merged
* `neutral_species_index`: the index of the neutral species which is the collision partner in the electron-neutral
    collisions for which approximate rates are being preserved.
* `n_rand_pairs`: the number of random particle pairs to choose
* `centered_at_mean`: whether to add a particle centered at the mean velocity of the system of particles
* `v_multipliers`: 
"""
function compute_lhs_particles_additional_rate_preserving!(rng, col_index, nnls_merging, lhs_matrix,
                                           interaction, electron_neutral_interactions, computed_cs, 
                                           particles, pia, cell, species, neutral_species_index,
                                           n_rand_pairs, centered_at_mean, v_multipliers)

    if centered_at_mean      
        n_moms = nnls_merging.n_moments
        v0 = norm(nnls_merging.v0)
        # a particle centered at local zero
        @inbounds for n_mom in 1:n_moms
            lhs_matrix[n_mom, col_index] = 0.0
        end

        compute_cross_sections_only!(computed_cs, interaction, v0, electron_neutral_interactions, neutral_species_index)
        @inbounds lhs_matrix[nnls_merging.n_moments+1, col_index] = get_cs_elastic(electron_neutral_interactions, computed_cs, neutral_species_index) * v0
        @inbounds lhs_matrix[nnls_merging.n_moments+2, col_index] = get_cs_ionization(electron_neutral_interactions, computed_cs, neutral_species_index) * v0

        col_index += 1
    end

    for v_mult in v_multipliers
        for i in 1:8
            vxs = vx_sign(i)
            vys = vy_sign(i)
            vzs = vz_sign(i)

            vx = check_speed_bound(vxs * nnls_merging.Ex, nnls_merging.minvx, nnls_merging.maxvx, v_mult)
            vy = check_speed_bound(vys * nnls_merging.Ey, nnls_merging.minvy, nnls_merging.maxvy, v_mult)
            vz = check_speed_bound(vzs * nnls_merging.Ez, nnls_merging.minvz, nnls_merging.maxvz, v_mult)
            @inbounds for n_mom in 1:nnls_merging.n_moments
                lhs_matrix[n_mom, col_index] = ccm(vx, vy, vz,
                                                             0.0, 0.0, 0.0,
                                                             nnls_merging.mim[n_mom])
            end

            @inbounds g = sqrt((nnls_merging.v0[1] + vx)^2 + (nnls_merging.v0[2] + vy)^2 + (nnls_merging.v0[3] + vz)^2)
            compute_cross_sections_only!(computed_cs, interaction, g, electron_neutral_interactions, neutral_species_index)
            @inbounds lhs_matrix[nnls_merging.n_moments+1, col_index] = get_cs_elastic(electron_neutral_interactions, computed_cs, neutral_species_index) * g
            @inbounds lhs_matrix[nnls_merging.n_moments+2, col_index] = get_cs_ionization(electron_neutral_interactions, computed_cs, neutral_species_index) * g

            col_index += 1
        end
    end

    for _ in 1:n_rand_pairs
        @inbounds i = floor(Int64, rand(rng, Float64) * pia.indexer[cell, species].n_local)
        @inbounds k = floor(Int64, rand(rng, Float64) * pia.indexer[cell, species].n_local)

        @inbounds i = map_cont_index(pia.indexer[cell, species], i)
        @inbounds k = map_cont_index(pia.indexer[cell, species], k)

        if (i != k)
            @inbounds vx = 0.5 * (particles[i].v[1] + particles[k].v[1])
            @inbounds vy = 0.5 * (particles[i].v[2] + particles[k].v[2])
            @inbounds vz = 0.5 * (particles[i].v[3] + particles[k].v[3])
            @inbounds for n_mom in 1:n_moms
                lhs_matrix[n_mom, col_index] = ccm(vx, vy, vz,
                                                   nnls_merging.v0[1], nnls_merging.v0[2], nnls_merging.v0[3],
                                                   nnls_merging.mim[n_mom])
            end

            @inbounds g = sqrt((nnls_merging.v0[1] + vx)^2 + (nnls_merging.v0[2] + vy)^2 + (nnls_merging.v0[3] + vz)^2)
            compute_cross_sections_only!(computed_cs, interaction, g, electron_neutral_interactions, neutral_species_index)
            @inbounds lhs_matrix[nnls_merging.n_moments+1, col_index] = get_cs_elastic(electron_neutral_interactions, computed_cs, neutral_species_index) * g
            @inbounds lhs_matrix[nnls_merging.n_moments+2, col_index] = get_cs_ionization(electron_neutral_interactions, computed_cs, neutral_species_index) * g
            
        end
        col_index += 1
    end
end

"""
    scale_lhs_rhs!(nnls_merging, lhs_matrix)

Scale the LHS and RHS of the NNLS system using the reference velocity ``v_{ref}``. Each moment is scaled
by ``(1/v_{ref})^{n_{tot}}``, where ``n_{tot}`` is the total order of the moment (i.e. for
a moment with multi-index `(i,j,k)` the total order is `i+j+k`).

# Positional arguments
* `nnls_merging`: the `NNLSMerge` instance
* `lhs_matrix`: the matrix of the LHS
"""
function scale_lhs_rhs!(nnls_merging, lhs_matrix)
    n_moms = nnls_merging.n_moments
    @inbounds for n_mom in 1:n_moms
        ref_val = nnls_merging.inv_vref^nnls_merging.tot_order[n_mom]
        lhs_matrix[n_mom, :] .*= ref_val
        nnls_merging.rhs_vector[n_mom] *= ref_val
    end
end

"""
    scale_lhs_rhs_rate_preserving!(nnls_merging, lhs_matrix, ref_cs_elastic, ref_cs_ion)

Scale the LHS and RHS of the NNLS system for the rate-preserving electron merging
using the reference velocity ``v_{ref}`` and
reference elastic scattering and ionization cross-sections. Each moment is scaled
by ``(1/v_{ref})^{n_{tot}}``, where ``n_{tot}`` is the total order of the moment (i.e. for
a moment with multi-index `(i,j,k)` the total order is `i+j+k`).
The entries in the LHS and RHS corresponding to the rates are scaled by
``1/(v_{ref} * \\sigma_{p})``, where ``\\sigma_{p}`` is the cross-section of
process ``p``. 

# Positional arguments
* `nnls_merging`: the `NNLSMerge` instance
* `lhs_matrix`: the matrix of the LHS
* `ref_cs_elastic`: the reference elastic scattering cross-section
* `ref_cs_ion`: the reference ionization cross-section
"""
function scale_lhs_rhs_rate_preserving!(nnls_merging, lhs_matrix, ref_cs_elastic, ref_cs_ion)
    scale_lhs_rhs!(nnls_merging, lhs_matrix)
    @inbounds lhs_matrix[nnls_merging.n_moments+1, :] .*= nnls_merging.inv_vref / ref_cs_elastic
    @inbounds lhs_matrix[nnls_merging.n_moments+2, :] .*= nnls_merging.inv_vref / ref_cs_ion

    @inbounds nnls_merging.rhs_vector[nnls_merging.n_moments+1] *= nnls_merging.inv_vref / ref_cs_elastic
    @inbounds nnls_merging.rhs_vector[nnls_merging.n_moments+2] *= nnls_merging.inv_vref / ref_cs_ion
end

"""
    compute_post_merge_particles_nnls!(lhs_ncols, lhs_matrix,
                                            nnls_merging, particles, pia, cell, species)

Compute post-merge particles based on the solution of the NNLS problem. This will replace
the particles with the post-merge ones and delete any extraneous particles.

# Positional arguments
* `lhs_ncols`: the number of columns in the LHS matrix
* `lhs_matrix`: the LHS matrix of the NNLS system
* `nnls_merging`: the `NNLSMerge` instance
* `particles`: the `ParticleVector` instance containing the particles that are being merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the grid cell in which particles are being merged
* `species`: the index of the species being merged
"""
function compute_post_merge_particles_nnls!(lhs_ncols, lhs_matrix,
                                            nnls_merging, particles, pia, cell, species)

    curr_particle_index = 0

    # lhs matrix has row 1 of ones (mass conservation)
    # row 2 are the vx components (Vx conservation)
    # row 3 are the vy components (Vy conservation)
    # row 4 are the vz components (Vz conservation)
    @inbounds for j in 1:lhs_ncols
        if nnls_merging.work.x[j] > 0.0
            i = map_cont_index(pia.indexer[cell,species], curr_particle_index)
            curr_particle_index += 1
            particles[i].w = nnls_merging.work.x[j] * nnls_merging.w_total
            particles[i].v = nnls_merging.v0 + SVector{3, Float64}(lhs_matrix[2, j],
                                                                   lhs_matrix[3, j],
                                                                   lhs_matrix[4, j]) * nnls_merging.vref
        end
    end

    @inbounds old_count = pia.indexer[cell,species].n_local
    n_particles_to_delete = old_count - curr_particle_index

    # if we delete from particles in last cell AND we delete less particles than were in group 2
    # then continuity is not broken
    # !(A && B) == !A || !B
    @inbounds if !(cell == size(pia.indexer)[1]) || (n_particles_to_delete > pia.indexer[cell,species].n_group2)
        pia.contiguous[species] = false
    end

    for _ in 1:n_particles_to_delete
        delete_particle_end!(particles, pia, cell, species)
    end
end

"""
    merge_nnls_based!(rng, nnls_merging, particles, pia, cell, species, vref; n_rand_pairs=0, max_err=1e-11,
                      centered_at_mean=true, v_multipliers=[0.5, 1.0], iteration_mult=2)

Perform NNLS-based merging. To improve stability, additional fictitious particles can be created
by randomly choosing particle pairs and creating new particles with a velocity and position that is a mean
of the velocity and position of the randomly chosen pair; these particles are added as additional columns to the matrix.

# Positional arguments
* `rng`: the random number generator instance
* `nnls_merging`: the `NNLSMerge` instance
* `particles`: the `ParticleVector` instance containing the particles to be merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the grid cell in which particles are being merged
* `species`: the index of the species being merged
* `vref`: the reference velocity used to scale the velocities

# Keyword arguments
* `n_rand_pairs`: the number of additional random pairs to sample to create additional entries in the matrix
    to potentially improve the stability of the algorithm
* `max_err`: maximum allowed value of the residual of the NNLS system
* `centered_at_mean`: whether to add a particle centered at the mean velocity of the system of particles
* `v_multipliers`: 
* `iteration_mult`: the number by which the number of columns of the NNLS system matrix is multiplied, this gives the maximum
    number of iterations of the NNLS algorithm

# Returns
If the residual exceeds `max_err` or
the number of non-zero elements in the solution vector is equal to the original number of particles,
`-1` is returned to signify a failure of the merging algorithm.
"""
function merge_nnls_based!(rng, nnls_merging, particles, pia, cell, species, vref;
                           n_rand_pairs=0, max_err=1e-11, centered_at_mean=true, v_multipliers=[0.5, 1.0], iteration_mult=2)

    # create LHS matrix
    n_add = centered_at_mean ? 1 : 0
    n_add += 8 * length(v_multipliers)
    @inbounds lhs_ncols = pia.indexer[cell, species].n_local + n_add + n_rand_pairs
    lhs_matrix = zeros(nnls_merging.n_moments, lhs_ncols)
    nnls_merging.vref = vref
    nnls_merging.inv_vref = 1.0 / vref

    # create LHS matrix and fill RHS vector using existing particles
    col_index = compute_lhs_and_rhs!(nnls_merging, lhs_matrix, particles, pia, cell, species)
    # and add more columns
    compute_lhs_particles_additional!(rng, col_index, nnls_merging, lhs_matrix,
                                      particles, pia, cell, species, n_rand_pairs,
                                      centered_at_mean, v_multipliers)
    scale_lhs_rhs!(nnls_merging, lhs_matrix)

    load!(nnls_merging.work, lhs_matrix, nnls_merging.rhs_vector)
    solve!(nnls_merging.work, iteration_mult * size(lhs_matrix, 2))

    nonzero = 0
    @inbounds for i in 1:lhs_ncols
        if nnls_merging.work.x[i] > 0.0
            nonzero += 1
        end
    end

    if nnls_merging.work.rnorm > max_err
        return -1
    end

    @inbounds if nonzero == pia.indexer[cell, species].n_local
        return -1
    end

    nnls_merging.work.x =  nnls_merging.work.x / sum(nnls_merging.work.x)
    
    compute_post_merge_particles_nnls!(lhs_ncols, lhs_matrix,
                                       nnls_merging, particles, pia, cell, species)
    return 1
end

"""
    merge_nnls_based_rate_preserving!(rng, nnls_merging,
                                           interaction, electron_neutral_interactions, computed_cs,
                                           particles, pia, cell, species, neutral_species_index,
                                           vref, ref_cs_elatic, ref_cs_ion; n_rand_pairs=0, max_err=1e-11,
                                           centered_at_mean=true, v_multipliers=[0.5, 1.0], iteration_mult=2)

Perform NNLS-based merging of electrons that conserves approximate elastic scattering and electron-impact ionization rates.
To improve stability, additional fictitious particles can be created
by randomly choosing particle pairs and creating new particles with a velocity and position that is a mean
of the velocity and position of the randomly chosen pair; these particles are added as additional columns to the matrix.

# Positional arguments
* `rng`: the random number generator instance
* `nnls_merging`: the `NNLSMerge` instance
* `interaction`: the `Interaction` instance describing the electron-neutral interaction being considered
* `electron_neutral_interactions`:  the `ElectronNeutralInteractions` instance storing the tabulated cross-section
    data used to compute the rates
* `computed_cs`: the vector of `ComputedCrossSection` instances in which the computed values will be stored
* `particles`: the `ParticleVector` instance containing the particles to be merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the grid cell in which particles are being merged
* `species`: the index of the species being merged
* `neutral_species_index`: the index of the neutral species which is the collision partner in the electron-neutral
    collisions for which approximate rates are being preserved.
* `vref`: the reference velocity used to scale the velocities
* `ref_cs_elatic`: the reference elastic scattering cross-section used to scale the rates
* `ref_cs_ion`: the reference electron-impact ionization cross-section used to scale the rates

# Keyword arguments
* `n_rand_pairs`: the number of additional random pairs to sample to create additional entries in the matrix
    to potentially improve the stability of the algorithm
* `max_err`: maximum allowed value of the residual of the NNLS system
* `centered_at_mean`: whether to add a particle centered at the mean velocity of the system of particles
* `v_multipliers`: 
* `iteration_mult`: the number by which the number of columns of the NNLS system matrix is multiplied, this gives the maximum
    number of iterations of the NNLS algorithm

# Returns
If the residual exceeds `max_err` or
the number of non-zero elements in the solution vector is equal to the original number of particles,
`-1` is returned to signify a failure of the merging algorithm.
"""
function merge_nnls_based_rate_preserving!(rng, nnls_merging,
                                           interaction, electron_neutral_interactions, computed_cs,
                                           particles, pia, cell, species, neutral_species_index,
                                           vref, ref_cs_elatic, ref_cs_ion; n_rand_pairs=0, max_err=1e-11,
                                           centered_at_mean=true, v_multipliers=[0.5, 1.0], iteration_mult=2)

    # create LHS matrix
    n_add = centered_at_mean ? 1 : 0
    n_add += 8 * length(v_multipliers)
    @inbounds lhs_ncols = pia.indexer[cell, species].n_local + n_add + n_rand_pairs
    lhs_matrix = zeros(nnls_merging.n_moments+2, lhs_ncols)
    nnls_merging.vref = vref
    nnls_merging.inv_vref = 1.0 / vref

    # create LHS matrix and fill RHS vector using existing particles
    @inbounds col_index = compute_lhs_and_rhs_rate_preserving!(nnls_merging, lhs_matrix,
                                                     interaction[species,neutral_species_index],
                                                     electron_neutral_interactions, computed_cs, 
                                                     particles, pia, cell, species, neutral_species_index)
    # and add more columns
    @inbounds compute_lhs_particles_additional_rate_preserving!(rng, col_index, nnls_merging, lhs_matrix,
                                      interaction[species,neutral_species_index], electron_neutral_interactions, computed_cs,
                                      particles, pia, cell, species, neutral_species_index, n_rand_pairs,
                                      centered_at_mean, v_multipliers)
    scale_lhs_rhs_rate_preserving!(nnls_merging, lhs_matrix, ref_cs_elatic, ref_cs_ion)

    load!(nnls_merging.work, lhs_matrix, nnls_merging.rhs_vector)
    solve!(nnls_merging.work, iteration_mult * size(lhs_matrix, 2))

    nonzero = 0
    @inbounds for i in 1:lhs_ncols
        if nnls_merging.work.x[i] > 0.0
            nonzero += 1
        end
    end

    if nnls_merging.work.rnorm > max_err
        return -1
    end

    @inbounds if nonzero == pia.indexer[cell, species].n_local
        return -1
    end

    nnls_merging.work.x = nnls_merging.work.x / sum(nnls_merging.work.x)
    
    compute_post_merge_particles_nnls!(lhs_ncols, lhs_matrix,
                                       nnls_merging, particles, pia, cell, species)
    return 1
end

end