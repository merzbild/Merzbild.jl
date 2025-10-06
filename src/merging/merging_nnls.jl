@muladd begin

"""
    NNLSMerge
Struct for keeping track of merging-related quantities for NNLS-based merging.

# Fields
* `v0`: vector of the mean velocity of the particles
* `x0`: vector of the mean position of the particles
* `vref`: reference velocity magnitude for scaling
* `inv_vref`: inverse of reference velocity magnitude for scaling
* `Ex`: standard deviation of x velocity of particles
* `Ey`: standard deviation of y velocity of particles
* `Ez`: standard deviation of z velocity of particles
* `Epx`: standard deviation of x position of particles
* `Epy`: standard deviation of y position of particles
* `Epz`: standard deviation of z position of particles
* `w_total`: total computational of the particles
* `minvx`: minimum x velocity of the particles
* `maxvx`: maximum x velocity of the particles
* `minvy`: minimum y velocity of the particles
* `maxvy`: maximum y velocity of the particles
* `minvz`: minimum z velocity of the particles
* `maxvz`: maximum z velocity of the particles
* `scalevx`: value used to scale the x velocity of particles
* `scalevy`: value used to scale the y velocity of particles
* `scalevz`: value used to scale the z velocity of particles
* `scalex`: value used to scale the x position of particles
* `scaley`: value used to scale the y position of particles
* `scalez`: value used to scale the z position of particles
* `n_total_conserved`: total number of moments conserved
* `n_moments_vel`: number of velocity moments to preserve
* `rhs_vector`: vector of computed moments
* `residual`: residual of solution
* `mim`: vector of 3-tuples of multi-indices for the velocity moments to preserve
* `tot_order`: vector of total orders of the velocity moments to preserve
* `n_moments_pos`: number of spatial moments to preserve
* `mim_pos`: vector of 3-tuples of multi-indices for the spatial moments to preserve
* `tot_order_pos`: vector of total orders of the spatial moments to preserve
* `pos_i_x`: index of the spatial moment corresponding to preservation of the center of mass in the x direction
* `pos_i_y`: index of the spatial moment corresponding to preservation of the center of mass in the y direction
* `pos_i_z`: index of the spatial moment corresponding to preservation of the center of mass in the z direction
* `lhs_matrices`: Vector of matrices for the LHS of the NNLS system
* `lhs_matrix_ncols_start`: the number of columns in the first matrix in `lhs_matrices`
* `lhs_matrix_ncols_end`: the number of columns in the last matrix in `lhs_matrices`
* `column_norms`: vector of vectors of the column-wise inverse norms of the LHS matrices
* `work`: Vector of `NNLSWorkspace` instances
"""
mutable struct NNLSMerge
    v0::SVector{3,Float64}
    x0::SVector{3,Float64}
    vref::Float64  # used for scaling
    inv_vref::Float64  # used for scaling

    Ex::Float64  # for additional particles
    Ey::Float64
    Ez::Float64
    Epx::Float64 
    Epy::Float64
    Epz::Float64
    w_total::Float64
    minvx::Float64
    maxvx::Float64
    minvy::Float64
    maxvy::Float64
    minvz::Float64
    maxvz::Float64

    scalevx::Float64
    scalevy::Float64
    scalevz::Float64
    scalex::Float64
    scaley::Float64
    scalez::Float64

    n_total_conserved::Int32
    n_moments_vel::Int32
    rhs_vector::Vector{Float64}
    mim::Vector{Vector{Int32}}  # mult-index moments
    tot_order::Vector{Int32}
    n_moments_pos::Int32
    mim_pos::Vector{Vector{Int32}}  # mult-index moments
    tot_order_pos::Vector{Int32}
    pos_i_x::Int32
    pos_i_y::Int32
    pos_i_z::Int32

    lhs_matrices::Vector{Matrix{Float64}}
    lhs_matrix_ncols_start::Int32
    lhs_matrix_ncols_end::Int32
    column_norms::Vector{Vector{Float64}}

    # work::NNLSWorkspace #
    work::Vector{NNLSWorkspace}

    @doc """
        NNLSMerge(multi_index_moments, init_np; rate_preserving=false, multi_index_moments_pos=[], matrix_ncol_nprealloc=0)

    Create NNLS-based merging. Mass, momentum, directional energy are always conserved:
    if not in the list `multi_index_moments` of moments to preserve, the corresponding
    moment multi-indices will be added automatically. These indices are `(0,0,0)` for mass,
    `(1,0,0)`, `(0,1,0)`, `(0,0,1)` for momentum, and `(2,0,0)`, `(0,2,0)`, `(0,0,2)`
    for directional energies. Spatial moments can be preserved by setting `multi_index_moments_pos`.
    If `multi_index_moments_pos` is non-empty, it is currently left to the user to include any relevant 1st order moments
    (corresponding center of mass conservation) i.e. `(1, 0, 0)`, `(0, 1, 0)`, `(0, 0, 1)`; otherwise
    the corresponding spatial moments including higher-order ones will not be conserved.
    By default, a single `NNLSWorkspace` is pre-allocated for the system, and no LHS matrices are pre-allocated.
    The number of columns in the LHS matrix is the number of particles to be merged (+ any fictitious particles);
    so it cannot be fixed in advance. By setting `matrix_ncol_nprealloc` to a value larger than 0,
    one can pre-allocate a range of matrices and
    corresponding workspaces with a fixed number of columns spanning `[init_np, init_np+matrix_ncol_nprealloc]`.
    In this case, `matrix_ncol_nprealloc+1` `NNLSWorkspace` instances are pre-allocated, with the last one being used
    in case the number of columns is not in the range (and the LHS matrix is constructed on the fly).
    A vector `column_norms` is also then pre-allocated, which is a vector of vectors, each containing the inverses of the
    column-wise norms of the LHS matrices, used for their scaling.

    # Positional arguments
    * `multi_index_moments`: vector of mixed moments to preserve of the form `[(i1, j1, k1), (i2, j2, k2), ...]``
    * `init_np`: assumption on pre-merge number of particles to pre-allocate memory for
    
    Keyword arguments:
    * `rate_preserving`: used for rate-preserving merging of electrons, preserves approximate elastic collision and ionization rates
    * `multi_index_moments_pos`: list of spatial moments to preserve
    * `matrix_ncol_nprealloc`: number of LHS matrices and NNLS workspaces with a fixed number of columns to pre-allocate
    """
    function NNLSMerge(multi_index_moments, init_np; rate_preserving=false, multi_index_moments_pos=[], matrix_ncol_nprealloc=0)
        add_length = 0
        if rate_preserving
            add_length = 2
        end

        base_moments = base_multi_index_moments()
        filtered_mim = [x for x in multi_index_moments if !(x in base_moments)]
        n_total_conserved = length(base_moments) + length(filtered_mim) + add_length + length(multi_index_moments_pos)
        append!(base_moments, filtered_mim)
        tot_order = [sum(mim) for mim in base_moments]
        tot_order_pos = [sum(mim) for mim in multi_index_moments_pos]

        pos_i_x = -1
        pos_i_y = -1
        pos_i_z = -1
        
        for i in 1:length(multi_index_moments_pos)
            if (multi_index_moments_pos[i][1] == 1) && (multi_index_moments_pos[i][2] == 0) && (multi_index_moments_pos[i][3] == 0)
                pos_i_x = i
            elseif (multi_index_moments_pos[i][1] == 0) && (multi_index_moments_pos[i][2] == 1) && (multi_index_moments_pos[i][3] == 0)
                pos_i_y = i
            elseif (multi_index_moments_pos[i][1] == 0) && (multi_index_moments_pos[i][2] == 0) && (multi_index_moments_pos[i][3] == 1)
                pos_i_z = i
            end 
        end

        matrices_preallocated = Vector{Matrix{Float64}}([])
        column_norms = Vector{Vector{Float64}}([])
        if matrix_ncol_nprealloc > 0
            for i in init_np:init_np+matrix_ncol_nprealloc
                push!(matrices_preallocated, zeros(n_total_conserved, i))
                push!(column_norms, ones(i))
            end
        end

        nnls_ws_preallocated = Vector{NNLSWorkspace}([])
        if matrix_ncol_nprealloc > 0
            for i in init_np:init_np+matrix_ncol_nprealloc+1
                push!(nnls_ws_preallocated, NNLSWorkspace(zeros(n_total_conserved, i), zeros(n_total_conserved)))
            end
        else
            push!(nnls_ws_preallocated, NNLSWorkspace(zeros(n_total_conserved, init_np), zeros(n_total_conserved)))
        end
        
        return new(SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0),
                   1.0, 1.0, # vref, inv_vref
                   0.0, 0.0, 0.0,   # Ex Ey Ez
                   0.0, 0.0, 0.0,   # Epx Epy Epz
                   0.0, # w_total
                   0.0, 0.0,  # min max vx
                   0.0, 0.0,  # min max vy
                   0.0, 0.0,  # min max vz
                   0.0, 0.0, 0.0,  # scale v
                   0.0, 0.0, 0.0,  # scale x
                   n_total_conserved,
                   length(base_moments),
                   zeros(n_total_conserved),
                   base_moments, tot_order,
                   length(multi_index_moments_pos), multi_index_moments_pos, tot_order_pos,
                   pos_i_x, pos_i_y, pos_i_z,
                   matrices_preallocated, init_np, init_np+matrix_ncol_nprealloc,
                   column_norms,
                   nnls_ws_preallocated)
                #    matrices_preallocated, init_np, init_np+matrix_ncol_nprealloc,
                #    NNLSWorkspace(zeros(n_total_conserved, init_np), zeros(n_total_conserved)))
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
        nnls_merging.x0 = nnls_merging.x0 + particles[i].w * particles[i].x

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
            nnls_merging.x0 = nnls_merging.x0 + particles[i].w * particles[i].x

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
    nnls_merging.x0 = nnls_merging.x0 / nnls_merging.w_total
    @inbounds nnls_merging.minvx -= nnls_merging.v0[1]
    @inbounds nnls_merging.minvy -= nnls_merging.v0[2]
    @inbounds nnls_merging.minvz -= nnls_merging.v0[3]

    @inbounds nnls_merging.maxvx -= nnls_merging.v0[1]
    @inbounds nnls_merging.maxvy -= nnls_merging.v0[2]
    @inbounds nnls_merging.maxvz -= nnls_merging.v0[3]
end

"""
    ccm(v, v0, mim)

Compute unweighted central velocity or spatial moment.

# Positional arguments
* `v`: the 3-dimensional velocity / position vector
* `v0`: the 3-dimensional mean velocity / position vector
* `mim`: the 3-dimensional multi-index

# Returns
Computed unweighted central moment.
"""
@inline function ccm(v, v0, mim)
    @inbounds return (v[1] - v0[1])^mim[1] * (v[2] - v0[2])^mim[2] * (v[3] - v0[3])^mim[3]
end

"""
    ccm(vx, vy, vz, vx0, vy0, vz0, mim)

Compute unweighted central velocity or spatial moment.

# Positional arguments
* `vx`: the x component of the velocity / position
* `vy`: the y component of the velocity / position
* `vz`: the z component of the velocity / position
* `vx0`: the x component of the mean velocity / position
* `vy0`: the y component of the mean velocity / position
* `vz0`: the z component of the mean velocity / position
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
    n_moms = nnls_merging.n_moments_vel
    
    nnls_merging.Ex = 0.0
    nnls_merging.Ey = 0.0
    nnls_merging.Ez = 0.0
    nnls_merging.Epx = 0.0
    nnls_merging.Epy = 0.0
    nnls_merging.Epz = 0.0
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

        nnls_merging.Epx = nnls_merging.Epx + particles[i].w * (particles[i].x[1] - nnls_merging.x0[1])^2
        nnls_merging.Epy = nnls_merging.Epy + particles[i].w * (particles[i].x[2] - nnls_merging.x0[2])^2
        nnls_merging.Epz = nnls_merging.Epz + particles[i].w * (particles[i].x[3] - nnls_merging.x0[3])^2

        for n_mom in 1:n_moms
            tmp_ccm = ccm(particles[i].v, nnls_merging.v0, nnls_merging.mim[n_mom])
            nnls_merging.rhs_vector[n_mom] = nnls_merging.rhs_vector[n_mom] + particles[i].w * tmp_ccm
            lhs_matrix[n_mom, col_index] = tmp_ccm
        end
        for n_mom in 1:nnls_merging.n_moments_pos
            tmp_ccm = ccm(particles[i].x, nnls_merging.x0, nnls_merging.mim_pos[n_mom])
            nnls_merging.rhs_vector[n_moms+n_mom] = nnls_merging.rhs_vector[n_moms+n_mom] + particles[i].w * tmp_ccm
            lhs_matrix[n_moms+n_mom, col_index] = tmp_ccm
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

            nnls_merging.Epx = nnls_merging.Epx + particles[i].w * (particles[i].x[1] - nnls_merging.x0[1])^2
            nnls_merging.Epy = nnls_merging.Epy + particles[i].w * (particles[i].x[2] - nnls_merging.x0[2])^2
            nnls_merging.Epz = nnls_merging.Epz + particles[i].w * (particles[i].x[3] - nnls_merging.x0[3])^2
            # w_total += particles[i].w
            for n_mom in 1:nnls_merging.n_moments_vel
                tmp_ccm = ccm(particles[i].v, nnls_merging.v0, nnls_merging.mim[n_mom])
                nnls_merging.rhs_vector[n_mom] = nnls_merging.rhs_vector[n_mom] + particles[i].w * tmp_ccm
                lhs_matrix[n_mom, col_index] = tmp_ccm
            end
            for n_mom in 1:nnls_merging.n_moments_pos
                tmp_ccm = ccm(particles[i].x, nnls_merging.x0, nnls_merging.mim_pos[n_mom])
                nnls_merging.rhs_vector[n_moms+n_mom] = nnls_merging.rhs_vector[n_moms+n_mom] + particles[i].w * tmp_ccm
                lhs_matrix[n_moms+n_mom, col_index] = tmp_ccm
            end
            col_index += 1
        end
    end

    for i in 1:nnls_merging.n_total_conserved
        nnls_merging.rhs_vector[i] /= nnls_merging.w_total
    end

    nnls_merging.Ex /= nnls_merging.w_total
    nnls_merging.Ey /= nnls_merging.w_total
    nnls_merging.Ez /= nnls_merging.w_total
    nnls_merging.Ex = sqrt(nnls_merging.Ex)
    nnls_merging.Ey = sqrt(nnls_merging.Ey)
    nnls_merging.Ez = sqrt(nnls_merging.Ez)

    nnls_merging.Epx /= nnls_merging.w_total
    nnls_merging.Epy /= nnls_merging.w_total
    nnls_merging.Epz /= nnls_merging.w_total
    nnls_merging.Epx = sqrt(nnls_merging.Epx)
    nnls_merging.Epy = sqrt(nnls_merging.Epy)
    nnls_merging.Epz = sqrt(nnls_merging.Epz)

    return col_index
end

"""
    compute_lhs_and_rhs_rate_preserving!(nnls_merging, lhs_matrix, particles, pia, cell, species, extend)

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
* `extend`: enum of `CSExtend` type that sets how out-of-range energy values are treated when computing cross-sections
"""
function compute_lhs_and_rhs_rate_preserving!(nnls_merging, lhs_matrix,
                                              interaction, electron_neutral_interactions, computed_cs,
                                              particles, pia, cell, species, neutral_species_index, extend)
    n_moms = nnls_merging.n_moments_vel
    
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
            nnls_merging.rhs_vector[n_mom] = nnls_merging.rhs_vector[n_mom] + particles[i].w * tmp_ccm
            lhs_matrix[n_mom, col_index] = tmp_ccm
        end

        g = norm(particles[i].v)
        compute_cross_sections_only!(computed_cs, interaction, g, electron_neutral_interactions, neutral_species_index, extend)
        cse_g = get_cs_elastic(electron_neutral_interactions, computed_cs, neutral_species_index) * g
        csi_g = get_cs_ionization(electron_neutral_interactions, computed_cs, neutral_species_index) * g
        lhs_matrix[nnls_merging.n_moments_vel+1, col_index] = cse_g
        lhs_matrix[nnls_merging.n_moments_vel+2, col_index] = csi_g
        nnls_merging.rhs_vector[nnls_merging.n_moments_vel+1] = nnls_merging.rhs_vector[nnls_merging.n_moments_vel+1] + cse_g * particles[i].w
        nnls_merging.rhs_vector[nnls_merging.n_moments_vel+2] = nnls_merging.rhs_vector[nnls_merging.n_moments_vel+2] + csi_g * particles[i].w

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
            for n_mom in 1:nnls_merging.n_moments_vel
                tmp_ccm = ccm(particles[i].v, nnls_merging.v0, nnls_merging.mim[n_mom])
                nnls_merging.rhs_vector[n_mom] = nnls_merging.rhs_vector[n_mom] + particles[i].w * tmp_ccm
                lhs_matrix[n_mom, col_index] = tmp_ccm
            end

            g = norm(particles[i].v)
            compute_cross_sections_only!(computed_cs, interaction, g, electron_neutral_interactions, neutral_species_index, extend)
            cse_g = get_cs_elastic(electron_neutral_interactions, computed_cs, neutral_species_index) * g
            csi_g = get_cs_ionization(electron_neutral_interactions, computed_cs, neutral_species_index) * g
            lhs_matrix[nnls_merging.n_moments_vel+1, col_index] = cse_g
            lhs_matrix[nnls_merging.n_moments_vel+2, col_index] = csi_g
            nnls_merging.rhs_vector[nnls_merging.n_moments_vel+1] = nnls_merging.rhs_vector[nnls_merging.n_moments_vel+1] + cse_g * particles[i].w
            nnls_merging.rhs_vector[nnls_merging.n_moments_vel+2] = nnls_merging.rhs_vector[nnls_merging.n_moments_vel+2] + csi_g * particles[i].w

            col_index += 1
        end
    end

    for i in 1:nnls_merging.n_total_conserved
        nnls_merging.rhs_vector[i] /= nnls_merging.w_total
    end

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
* `v_multipliers`: the multipliers for the velocity variances of the particles to add aditional particles to the system
"""
function compute_lhs_particles_additional!(rng, col_index, nnls_merging, lhs_matrix,
                                           particles, pia, cell, species,
                                           n_rand_pairs, centered_at_mean, v_multipliers)
    n_moms = nnls_merging.n_moments_vel

    if centered_at_mean
        
        # a particle centered at local zero
        @inbounds for n_mom in 1:n_moms
            lhs_matrix[n_mom, col_index] = 0.0
        end
        @inbounds for n_mom in 1:nnls_merging.n_moments_pos
            lhs_matrix[n_moms + n_mom, col_index] = 0.0
        end
        col_index += 1
    end

    # add 8 additional particles with velocity moments and position in center
    for v_mult in v_multipliers
        for i in 1:8
            vxs = vx_sign(i)
            vys = vy_sign(i)
            vzs = vz_sign(i)

            vx = check_speed_bound(vxs * nnls_merging.Ex, nnls_merging.minvx, nnls_merging.maxvx, v_mult)
            vy = check_speed_bound(vys * nnls_merging.Ey, nnls_merging.minvy, nnls_merging.maxvy, v_mult)
            vz = check_speed_bound(vzs * nnls_merging.Ez, nnls_merging.minvz, nnls_merging.maxvz, v_mult)
            @inbounds for n_mom in 1:nnls_merging.n_moments_vel
                lhs_matrix[n_mom, col_index] = ccm(vx, vy, vz,
                                                   0.0, 0.0, 0.0,
                                                   nnls_merging.mim[n_mom])
            end
            @inbounds for n_mom in 1:nnls_merging.n_moments_pos
                lhs_matrix[n_moms + n_mom, col_index] = 0.0
            end
            col_index += 1
        end
    end

    @inbounds for _ in 1:n_rand_pairs
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
            @inbounds x = 0.5 * (particles[i].x[1] + particles[k].x[1])
            @inbounds y = 0.5 * (particles[i].x[2] + particles[k].x[2])
            @inbounds z = 0.5 * (particles[i].x[3] + particles[k].x[3])
            @inbounds for n_mom in 1:nnls_merging.n_moments_pos
                lhs_matrix[n_moms + n_mom, col_index] = ccm(x, y, z,
                                                   nnls_merging.x0[1], nnls_merging.x0[2], nnls_merging.x0[3],
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
                                           n_rand_pairs, centered_at_mean, v_multipliers, extend)

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
* `v_multipliers`: the multipliers for the velocity variances of the particles to add aditional particles to the system
* `extend`: enum of `CSExtend` type that sets how out-of-range energy values are treated when computing cross-sections
"""
function compute_lhs_particles_additional_rate_preserving!(rng, col_index, nnls_merging, lhs_matrix,
                                           interaction, electron_neutral_interactions, computed_cs, 
                                           particles, pia, cell, species, neutral_species_index,
                                           n_rand_pairs, centered_at_mean, v_multipliers, extend)
    n_moms = nnls_merging.n_moments_vel

    if centered_at_mean      
        v0 = norm(nnls_merging.v0)
        # a particle centered at local zero
        @inbounds for n_mom in 1:n_moms
            lhs_matrix[n_mom, col_index] = 0.0
        end

        compute_cross_sections_only!(computed_cs, interaction, v0, electron_neutral_interactions, neutral_species_index, extend)
        @inbounds lhs_matrix[nnls_merging.n_moments_vel+1, col_index] = get_cs_elastic(electron_neutral_interactions, computed_cs, neutral_species_index) * v0
        @inbounds lhs_matrix[nnls_merging.n_moments_vel+2, col_index] = get_cs_ionization(electron_neutral_interactions, computed_cs, neutral_species_index) * v0

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
            @inbounds for n_mom in 1:nnls_merging.n_moments_vel
                lhs_matrix[n_mom, col_index] = ccm(vx, vy, vz,
                                                             0.0, 0.0, 0.0,
                                                             nnls_merging.mim[n_mom])
            end

            @inbounds g = sqrt((nnls_merging.v0[1] + vx)^2 + (nnls_merging.v0[2] + vy)^2 + (nnls_merging.v0[3] + vz)^2)
            compute_cross_sections_only!(computed_cs, interaction, g, electron_neutral_interactions, neutral_species_index, extend)
            @inbounds lhs_matrix[nnls_merging.n_moments_vel+1, col_index] = get_cs_elastic(electron_neutral_interactions, computed_cs, neutral_species_index) * g
            @inbounds lhs_matrix[nnls_merging.n_moments_vel+2, col_index] = get_cs_ionization(electron_neutral_interactions, computed_cs, neutral_species_index) * g

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
            compute_cross_sections_only!(computed_cs, interaction, g, electron_neutral_interactions, neutral_species_index, extend)
            @inbounds lhs_matrix[nnls_merging.n_moments_vel+1, col_index] = get_cs_elastic(electron_neutral_interactions, computed_cs, neutral_species_index) * g
            @inbounds lhs_matrix[nnls_merging.n_moments_vel+2, col_index] = get_cs_ionization(electron_neutral_interactions, computed_cs, neutral_species_index) * g
            
        end
        col_index += 1
    end
end

"""
    scale_lhs_rhs_vref!(nnls_merging, lhs_matrix, lhs_ncols)

Scale the LHS and RHS of the NNLS system using the reference velocity ``v_{ref}``. Each moment is scaled
by ``(1/v_{ref})^{n_{tot}}``, where ``n_{tot}`` is the total order of the moment (i.e. for
a moment with multi-index `(i,j,k)` the total order is `i+j+k`).

# Positional arguments
* `nnls_merging`: the `NNLSMerge` instance
* `lhs_matrix`: the matrix of the LHS
* `lhs_ncols`: number of columns in the LHS matrix
"""
function scale_lhs_rhs_vref!(nnls_merging, lhs_matrix, lhs_ncols)
    n_moms = nnls_merging.n_moments_vel
    @inbounds for n_mom in 1:n_moms
        ref_val = nnls_merging.inv_vref^nnls_merging.tot_order[n_mom]
        for col in 1:lhs_ncols
            lhs_matrix[n_mom, col] *= ref_val
        end
        nnls_merging.rhs_vector[n_mom] *= ref_val
    end
    nnls_merging.scalevx = nnls_merging.vref
    nnls_merging.scalevy = nnls_merging.vref
    nnls_merging.scalevz = nnls_merging.vref
end

"""
    scale_lhs_rhs_variance!(nnls_merging, lhs_matrix, lhs_ncols)

Scale the LHS and RHS of the NNLS system using the computed variances of the particles velocities
in the ``x``, ``y``, ``z`` directions. If any of the variances is smaller than 1e-6, then
``v_{ref}`` is used.
Each moment with multi-index `(i,j,k)` is scaled
by ``(1/Ex)^{i}(1/Ey)^{j}(1/Ez)^{k}``.

# Positional arguments
* `nnls_merging`: the `NNLSMerge` instance
* `lhs_matrix`: the matrix of the LHS
* `lhs_ncols`: number of columns in the LHS matrix
"""
function scale_lhs_rhs_variance!(nnls_merging, lhs_matrix, lhs_ncols)
    n_moms = nnls_merging.n_moments_vel
    inv_ex = nnls_merging.Ex > 1e-6 ? 1.0 / nnls_merging.Ex : nnls_merging.inv_vref
    inv_ey = nnls_merging.Ey > 1e-6 ? 1.0 / nnls_merging.Ey : nnls_merging.inv_vref
    inv_ez = nnls_merging.Ez > 1e-6 ? 1.0 / nnls_merging.Ez : nnls_merging.inv_vref
    @inbounds for n_mom in 1:n_moms
        ref_val = (inv_ex^nnls_merging.mim[n_mom][1]) * (inv_ey^nnls_merging.mim[n_mom][2]) * (inv_ez^nnls_merging.mim[n_mom][3])
        for col in 1:lhs_ncols
            lhs_matrix[n_mom, col] *= ref_val
        end
        nnls_merging.rhs_vector[n_mom] *= ref_val
    end
    nnls_merging.scalevx = 1.0 / inv_ex
    nnls_merging.scalevy = 1.0 / inv_ey
    nnls_merging.scalevz = 1.0 / inv_ez
end


"""
    scale_lhs_rhs_spatial_variance!(nnls_merging, lhs_matrix, lhs_ncols)

Scale the spatial moments in the LHS and RHS of the NNLS system using the computed variances of the particles positions
in the ``x``, ``y``, ``z`` directions. If any of the variances is smaller than 1e-6, then
no scaling is done.
Each moment with multi-index `(i,j,k)` is scaled
by ``(1/Epx)^{i}(1/Epy)^{j}(1/Epz)^{k}``.

# Positional arguments
* `nnls_merging`: the `NNLSMerge` instance
* `lhs_matrix`: the matrix of the LHS
* `lhs_ncols`: number of columns in the LHS matrix
"""
function scale_lhs_rhs_spatial_variance!(nnls_merging, lhs_matrix, lhs_ncols)
    n_moms = nnls_merging.n_moments_vel
    inv_epx = nnls_merging.Epx > 1e-6 ? 1.0 / nnls_merging.Epx : 1.0
    inv_epy = nnls_merging.Epy > 1e-6 ? 1.0 / nnls_merging.Epy : 1.0
    inv_epz = nnls_merging.Epz > 1e-6 ? 1.0 / nnls_merging.Epz : 1.0
    @inbounds for n_mom in 1:nnls_merging.n_moments_pos
        ref_val = (inv_epx^nnls_merging.mim_pos[n_mom][1]) * (inv_epy^nnls_merging.mim_pos[n_mom][2]) * (inv_epz^nnls_merging.mim_pos[n_mom][3])
        for col in 1:lhs_ncols
            lhs_matrix[n_moms + n_mom, col] *= ref_val
        end
        nnls_merging.rhs_vector[n_moms + n_mom] *= ref_val
    end
    nnls_merging.scalex = 1.0 / inv_epx
    nnls_merging.scaley = 1.0 / inv_epy
    nnls_merging.scalez = 1.0 / inv_epz
end


"""
    scale_lhs_rhs!(nnls_merging, lhs_matrix, scaling, lhs_ncols)

Scale the LHS and RHS of the NNLS system using either the reference velocity or the computed variances of the velocity
in each direction.

# Positional arguments
* `nnls_merging`: the `NNLSMerge` instance
* `lhs_matrix`: the matrix of the LHS
* `scaling`: how to scale entries (`:vref` or `:variance`)
* `lhs_ncols`: number of columns in the LHS matrix
"""
function scale_lhs_rhs!(nnls_merging, lhs_matrix, scaling, lhs_ncols)
    if scaling == :variance
        scale_lhs_rhs_variance!(nnls_merging, lhs_matrix, lhs_ncols)
        scale_lhs_rhs_spatial_variance!(nnls_merging, lhs_matrix, lhs_ncols)
    elseif scaling == :vref
        scale_lhs_rhs_vref!(nnls_merging, lhs_matrix, lhs_ncols)
        scale_lhs_rhs_spatial_variance!(nnls_merging, lhs_matrix, lhs_ncols)
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
* `scaling`: how to scale entries (`:vref` or `:variance`)
* `lhs_ncols`: number of columns in the LHS matrix
"""
function scale_lhs_rhs_rate_preserving!(nnls_merging, lhs_matrix, ref_cs_elastic, ref_cs_ion, scaling, lhs_ncols)
    scale_lhs_rhs!(nnls_merging, lhs_matrix, scaling, lhs_ncols)

    scaler_el = nnls_merging.inv_vref / ref_cs_elastic
    scaler_ion = nnls_merging.inv_vref / ref_cs_ion

    @inbounds for col in 1:lhs_ncols
        lhs_matrix[nnls_merging.n_moments_vel+1, col] *= scaler_el
        lhs_matrix[nnls_merging.n_moments_vel+2, col] *= scaler_ion
    end

    @inbounds nnls_merging.rhs_vector[nnls_merging.n_moments_vel+1] *= scaler_el
    @inbounds nnls_merging.rhs_vector[nnls_merging.n_moments_vel+2] *= scaler_ion
end

"""
    compute_post_merge_particles_nnls!(lhs_ncols, lhs_matrix,
                                            nnls_merging, particles, pia, cell, species)

Compute post-merge particles based on the solution of the NNLS problem. This will replace
the particles with the post-merge ones and delete any extraneous particles.

# Positional arguments
* `nnls_merging`: the `NNLSMerge` instance
* `x`: the solution vector of the NNLS system
* `particles`: the `ParticleVector` instance containing the particles that are being merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the grid cell in which particles are being merged
* `species`: the index of the species being merged
* `lhs_ncols`: the number of columns in the LHS matrix
* `lhs_matrix`: the LHS matrix of the NNLS system
* `max_err`: maximum allowed value of the residual of the NNLS system
* `w_threshold`: the relative (w.r.t the total computational weight of the particles being merge) value of the computational weight below which particles are discarded
* `work_index`: the index of the `NNLSWorkspace` used to solve the NNLS system
* `column_norms`: the vector of the column-wise norms of the LHS matrix
"""
function compute_post_merge_particles_nnls!(nnls_merging::NNLSMerge, x::Vector{Float64}, particles,
                                            pia, cell, species, lhs_ncols, lhs_matrix,
                                            max_err, w_threshold, work_index, column_norms)
    @inbounds if nnls_merging.work[work_index].rnorm > max_err
        return -1
    end

    # @inbounds x = nnls_merging.work[work_index].x

    nonzero = 0
    non_discarded_weight = 0.0
    @inbounds for i in 1:lhs_ncols
        x[i] = x[i] * column_norms[i]
        # w =  # * column_norms[i]
        if x[i] > w_threshold
            nonzero += 1
            non_discarded_weight += x[i]
        end
    end
    # println("ndw $al $(typeof(x)) $(typeof(non_discarded_weight)) $(typeof(column_norms)) $(w_threshold)")

    @inbounds if nonzero >= pia.indexer[cell, species].n_local
        return -1
    end

    # sum(non-discared elements in x) should be 1.0
    # nnls_merging.work.x = nnls_merging.work.x / sum(nnls_merging.work.x)
    # and we also scale back to the original weight
    # al = @allocated scale_factor = nnls_merging.w_total / non_discarded_weight
    nnls_merging.w_total = nnls_merging.w_total / non_discarded_weight
    
    curr_particle_index = 0

    # lhs matrix has row 1 of ones (mass conservation)
    # row 2 are the vx components (Vx conservation)
    # row 3 are the vy components (Vy conservation)
    # row 4 are the vz components (Vz conservation)
    @inbounds for j in 1:lhs_ncols
        if x[j] > w_threshold
            i = map_cont_index(pia.indexer[cell,species], curr_particle_index)
            curr_particle_index += 1
            particles[i].w = x[j] * nnls_merging.w_total# * column_norms[j]
            particles[i].v = SVector{3, Float64}(nnls_merging.v0[1] + lhs_matrix[2, j] * nnls_merging.scalevx,
                                                 nnls_merging.v0[2] + lhs_matrix[3, j] * nnls_merging.scalevy,
                                                 nnls_merging.v0[3] + lhs_matrix[4, j] * nnls_merging.scalevz)
            
            px = nnls_merging.pos_i_x > 0.0 ? lhs_matrix[nnls_merging.n_moments_vel + nnls_merging.pos_i_x, j] * nnls_merging.scalex : 0.0
            py = nnls_merging.pos_i_y > 0.0 ? lhs_matrix[nnls_merging.n_moments_vel + nnls_merging.pos_i_y, j] * nnls_merging.scaley : 0.0
            pz = nnls_merging.pos_i_z > 0.0 ? lhs_matrix[nnls_merging.n_moments_vel + nnls_merging.pos_i_z, j] * nnls_merging.scalez : 0.0

            particles[i].x = SVector{3,Float64}(nnls_merging.x0[1] + px, nnls_merging.x0[2] + py, nnls_merging.x0[3] + pz)
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

    return 1
end

"""
    merge_nnls_based!(rng, nnls_merging, particles, pia, cell, species, vref; n_rand_pairs=0, max_err=1e-11,
                      centered_at_mean=true, v_multipliers=[0.5, 1.0], iteration_mult=2)

Perform NNLS-based merging.
The NNLS system is scaled to improve numerical stability, the scaling algorithm is set by the `scaling` parameter.
Even if scaling is done using the computed variances, `vref` might be used in case those variances are small.

To improve stability, additional fictitious particles can be created
by randomly choosing particle pairs and creating new particles with a velocity and position that is a mean
of the velocity and position of the randomly chosen pair; these particles are added as additional columns to the matrix.
An additional fictitious particle can be created at the local mean velocity via the `centered_at_mean` parameter.
Additionally, particles can be created with velocities that are a multiple of the computed velocity variance of the system
of pre-merge particles in each direction.
For a given variance multiplier (an entry in the `v_multipliers` parameter, which is a vector of multipliers),
8 particles are added, one in each octant. Each of the particles velocity components is given either by
the value of the multiplier times the velocity variance of that component (times +1 or -1 depending on the octant),
or, if this value is outside of the bounding box, the closest bounding value of the velocity component in that direction
is used instead and also multiplied by the variance multiplier.

# Positional arguments
* `rng`: the random number generator instance
* `nnls_merging`: the `NNLSMerge` instance
* `particles`: the `ParticleVector` instance containing the particles to be merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the grid cell in which particles are being merged
* `species`: the index of the species being merged

# Keyword arguments
* `vref`: the reference velocity used to scale the velocities
* `scaling`: how to scale entries in the LHS and RHS of the NNLS system - either based on
    the reference velocity `vref` (`scaling=:vref`)
    or on the computed variances in each direction (`scaling=:variance`)
* `n_rand_pairs`: the number of additional random pairs to sample to create additional entries in the matrix
    to potentially improve the stability of the algorithm
* `max_err`: maximum allowed value of the residual of the NNLS system
* `centered_at_mean`: whether to add a particle centered at the mean velocity of the system of particles
* `v_multipliers`: the multipliers for the velocity variances of the particles to add aditional particles to the system
* `iteration_mult`: the number by which the number of columns of the NNLS system matrix is multiplied, this gives the maximum
    number of iterations of the NNLS algorithm
* `w_threshold`: any particles with a relative weight smaller than this value will be discarded (and the weight of the remaining particles re-scaled)

# Returns
If the residual exceeds `max_err` or
the number of non-zero (or smaller than `w_threshold`) elements in the solution vector is equal to the original number of particles,
`-1` is returned to signify a failure of the merging algorithm.

# References
* G. Oblapenko, A Non-Negative Least Squares-based Approach for Moment-Preserving Particle Merging.
    [arXiv preprint, 2025](https://doi.org/10.48550/arXiv.2412.12354).
"""
function merge_nnls_based!(rng, nnls_merging, particles, pia, cell, species;
                           vref=1.0, scaling=:variance,
                           n_rand_pairs=0, max_err=1e-11, centered_at_mean=true, v_multipliers=[0.5, 1.0], iteration_mult=2, w_threshold=0.0)
    # create LHS matrix
    n_add = centered_at_mean ? 1 : 0
    n_add = n_add + 8 * length(v_multipliers)
    @inbounds lhs_ncols = pia.indexer[cell, species].n_local + n_add + n_rand_pairs

    if (lhs_ncols > nnls_merging.lhs_matrix_ncols_end) || (lhs_ncols < nnls_merging.lhs_matrix_ncols_start)
        indexer = nnls_merging.lhs_matrix_ncols_end - nnls_merging.lhs_matrix_ncols_start + 2

        # we either have pre-allocation [1,2,...[Workspace]]
        # [1,2,Workspace]
        # [Workspace] - we don't pre-allocate if range is exactly 1 particle
        if indexer == 2
            indexer = 1
        end
         
        nnls_merging.work[indexer].QA = zeros(nnls_merging.n_total_conserved, lhs_ncols)
        lhs_matrix = zeros(nnls_merging.n_total_conserved, lhs_ncols)
        resize!(nnls_merging.work[indexer].x, lhs_ncols)
        resize!(nnls_merging.work[indexer].w, lhs_ncols)
        resize!(nnls_merging.work[indexer].idx, lhs_ncols)
        column_norms = ones(lhs_ncols)
    else
        indexer = lhs_ncols - nnls_merging.lhs_matrix_ncols_start + 1
        lhs_matrix = nnls_merging.lhs_matrices[lhs_ncols - nnls_merging.lhs_matrix_ncols_start + 1]
        column_norms = nnls_merging.column_norms[indexer]
    end
    nnls_merging.vref = vref
    nnls_merging.inv_vref = 1.0 / vref

    # create LHS matrix and fill RHS vector using existing particles
    col_index = compute_lhs_and_rhs!(nnls_merging, nnls_merging.work[indexer].QA, particles, pia, cell, species)
    # and add more columns
    compute_lhs_particles_additional!(rng, col_index, nnls_merging, nnls_merging.work[indexer].QA,
                                      particles, pia, cell, species, n_rand_pairs,
                                      centered_at_mean, v_multipliers)


    scale_lhs_rhs!(nnls_merging, nnls_merging.work[indexer].QA, scaling, lhs_ncols)

    lhs_matrix .= nnls_merging.work[indexer].QA
    
    scale_columns!(nnls_merging.work[indexer].QA, column_norms)
    
    nnls_merging.work[indexer].Qb .= nnls_merging.rhs_vector
    solve!(nnls_merging.work[indexer], iteration_mult * size(nnls_merging.work[indexer].QA, 2), 1e-14)
    return compute_post_merge_particles_nnls!(nnls_merging, nnls_merging.work[indexer].x, particles, pia, cell, species, lhs_ncols, lhs_matrix,
                                                             max_err, w_threshold, indexer, column_norms)
end

"""
    merge_nnls_based_rate_preserving!(rng, nnls_merging,
                                           interaction, electron_neutral_interactions, computed_cs,
                                           particles, pia, cell, species, neutral_species_index,
                                           ref_cs_elatic, ref_cs_ion; scaling=:variance,
                                           vref=1.0, n_rand_pairs=0, max_err=1e-11,
                                           centered_at_mean=true, v_multipliers=[0.5, 1.0], iteration_mult=2,
                                           extend::CSExtend=CSExtendConstant)

Perform NNLS-based merging of electrons that conserves approximate elastic scattering and electron-impact ionization rates.
The NNLS system is scaled to improve numerical stability, the scaling algorithm is set by the `scaling` parameter.
Even if scaling is done using the computed variances, `vref` might be used in case those variances are small.

To improve stability, additional fictitious particles can be created
by randomly choosing particle pairs and creating new particles with a velocity and position that is a mean
of the velocity and position of the randomly chosen pair; these particles are added as additional columns to the matrix.
An additional fictitious particle can be created at the local mean velocity via the `centered_at_mean` parameter.
Additionally, particles can be created with velocities that are a multiple of the computed velocity variance of the system
of pre-merge particles in each direction.
For a given variance multiplier (an entry in the `v_multipliers` parameter, which is a vector of multipliers),
8 particles are added, one in each octant. Each of the particles velocity components is given either by
the value of the multiplier times the velocity variance of that component (times +1 or -1 depending on the octant),
or, if this value is outside of the bounding box, the closest bounding value of the velocity component in that direction
is used instead and also multiplied by the variance multiplier.
The reference velocity is also used in conjunction with the reference cross-sections to scale the parts of
the NNLS matrix and RHS corresponding to conservation of electron-neutral collision rates.

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
* `ref_cs_elatic`: the reference elastic scattering cross-section used to scale the rates
* `ref_cs_ion`: the reference electron-impact ionization cross-section used to scale the rates

# Keyword arguments
* `vref`: the reference velocity used to scale the velocities and the electron-neutral collision rates
* `scaling`: how to scale entries in the LHS and RHS of the NNLS system - either based on
    the reference velocity `vref` (`scaling=:vref`)
    or on the computed variances in each direction (`scaling=:variance`)
* `n_rand_pairs`: the number of additional random pairs to sample to create additional entries in the matrix
    to potentially improve the stability of the algorithm
* `max_err`: maximum allowed value of the residual of the NNLS system
* `centered_at_mean`: whether to add a particle centered at the mean velocity of the system of particles
* `v_multipliers`: the multipliers for the velocity variances of the particles to add aditional particles to the system
* `iteration_mult`: the number by which the number of columns of the NNLS system matrix is multiplied, this gives the maximum
    number of iterations of the NNLS algorithm
* `w_threshold`: any particles with a relative weight smaller than this value will be discarded (and the weight of the remaining particles re-scaled)
* `extend`: enum of `CSExtend` type that sets how out-of-range energy values are treated when computing cross-sections

# Returns
If the residual exceeds `max_err` or
the number of non-zero (or smaller than `w_threshold`) elements in the solution vector is equal to the original number of particles,
`-1` is returned to signify a failure of the merging algorithm.

# References
* G. Oblapenko, A Non-Negative Least Squares-based Approach for Moment-Preserving Particle Merging.
    [arXiv preprint, 2025](https://doi.org/10.48550/arXiv.2412.12354).
"""
function merge_nnls_based_rate_preserving!(rng, nnls_merging,
                                           interaction, electron_neutral_interactions, computed_cs,
                                           particles, pia, cell, species, neutral_species_index,
                                           ref_cs_elatic, ref_cs_ion; vref=1.0, scaling=:variance,
                                           n_rand_pairs=0, max_err=1e-11,
                                           centered_at_mean=true, v_multipliers=[0.5, 1.0], iteration_mult=2, w_threshold=0.0,
                                           extend::CSExtend=CSExtendConstant)

    # create LHS matrix
    n_add = centered_at_mean ? 1 : 0
    n_add = n_add + 8 * length(v_multipliers)
    @inbounds lhs_ncols = pia.indexer[cell, species].n_local + n_add + n_rand_pairs

    if (lhs_ncols > nnls_merging.lhs_matrix_ncols_end) || (lhs_ncols < nnls_merging.lhs_matrix_ncols_start)
        indexer = nnls_merging.lhs_matrix_ncols_end - nnls_merging.lhs_matrix_ncols_start + 2

        # we either have pre-allocation [1,2,...[Workspace]]
        # [1,2,Workspace]
        # [Workspace] - we don't pre-allocate if range is exactly 1 particle
        if indexer == 2
            indexer = 1
        end

        nnls_merging.work[indexer].QA = zeros(nnls_merging.n_total_conserved, lhs_ncols)
        lhs_matrix = zeros(nnls_merging.n_total_conserved, lhs_ncols)
        resize!(nnls_merging.work[indexer].x, lhs_ncols)
        resize!(nnls_merging.work[indexer].w, lhs_ncols)
        resize!(nnls_merging.work[indexer].idx, lhs_ncols)
        column_norms = zeros(lhs_ncols)
    else
        indexer = lhs_ncols - nnls_merging.lhs_matrix_ncols_start + 1
        lhs_matrix = nnls_merging.lhs_matrices[lhs_ncols - nnls_merging.lhs_matrix_ncols_start + 1]
        column_norms = nnls_merging.column_norms[indexer]
    end
    nnls_merging.vref = vref
    nnls_merging.inv_vref = 1.0 / vref

    # create LHS matrix and fill RHS vector using existing particles
    @inbounds col_index = compute_lhs_and_rhs_rate_preserving!(nnls_merging, lhs_matrix,
                                                     interaction[species,neutral_species_index],
                                                     electron_neutral_interactions, computed_cs, 
                                                     particles, pia, cell, species, neutral_species_index, extend)
    # and add more columns
    @inbounds compute_lhs_particles_additional_rate_preserving!(rng, col_index, nnls_merging, lhs_matrix,
                                      interaction[species,neutral_species_index], electron_neutral_interactions, computed_cs,
                                      particles, pia, cell, species, neutral_species_index, n_rand_pairs,
                                      centered_at_mean, v_multipliers, extend)
    scale_lhs_rhs_rate_preserving!(nnls_merging, lhs_matrix, ref_cs_elatic, ref_cs_ion, scaling, lhs_ncols)
    scale_columns!(lhs_matrix, column_norms)

    lhs_matrix .= nnls_merging.work[indexer].QA

    nnls_merging.work[indexer].Qb .= nnls_merging.rhs_vector
    solve!(nnls_merging.work[indexer], iteration_mult * size(nnls_merging.work[indexer].QA, 2), 1e-14)
    return compute_post_merge_particles_nnls!(nnls_merging, nnls_merging.work[indexer].x, particles, pia, cell, species, lhs_ncols, lhs_matrix,
                                              max_err, w_threshold, indexer, column_norms)
end
end