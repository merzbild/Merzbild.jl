@muladd begin

"""
    NNLSMerge{D}
Struct for keeping track of merging-related quantities for NNLS-based merging of particles with D-dimensional position vectors.

# Fields
* `v0`: vector of the mean velocity of the particles
* `x0`: D-dimensional vector of the mean position of the particles
* `vref`: reference velocity magnitude for scaling
* `inv_vref`: inverse of reference velocity magnitude for scaling
* `Ev`: vector of standard deviation of velocities of particles
* `Ex`: vector of standard deviation of positions of particles
* `w_total`: total computational of the particles
* `scalev`: vector of values used to scale the velocity of particles
* `scalex`: vector of values used to scale the positions of particles
* `n_total_conserved`: total number of moments conserved
* `n_moments_vel`: number of velocity moments to preserve
* `rhs_vector`: vector of computed moments
* `row_scale`: scratch vector of the row-wise scaling factors applied to the LHS matrix
* `mim`: vector of 3-tuples of multi-indices for the velocity moments to preserve
* `tot_order`: vector of total orders of the velocity moments to preserve
* `vel_powers`: scratch table of powers of the centered velocity components of a particle
* `n_moments_pos`: number of spatial moments to preserve
* `mim_pos`: vector of 3-tuples of multi-indices for the spatial moments to preserve
* `tot_order_pos`: vector of total orders of the spatial moments to preserve
* `pos_powers`: scratch table of powers of the centered position components of a particle
* `pos_i_x`: index of the spatial moment corresponding to preservation of the center of mass in the x direction
* `pos_i_y`: index of the spatial moment corresponding to preservation of the center of mass in the y direction
* `pos_i_z`: index of the spatial moment corresponding to preservation of the center of mass in the z direction
* `lhs_matrix_ncols_start`: the number of columns in the first pre-allocated matrix
* `lhs_matrix_ncols_end`: the number of columns in the last pre-allocated matrix
* `column_norms`: vector of vectors of the column-wise inverse norms of the LHS matrices
* `vel_pos_matrices`: vector of matrices of size `(3+D)xNp` that store the velocities and positions of the particles
* `column_norms_scratch`: column-wise inverse norms used when the number of columns is outside the pre-allocated range
* `vel_pos_matrix_scratch`: velocities and positions used when the number of columns is outside the pre-allocated range
* `work`: Vector of `NNLSWorkspace` instances
"""
mutable struct NNLSMerge{D}
    v0::SVector{3,Float64}
    x0::SVector{D,Float64}
    vref::Float64  # used for scaling
    inv_vref::Float64  # used for scaling

    Ev::SVector{3,Float64}
    Ex::SVector{D,Float64}
    w_total::Float64
    scalev::SVector{3,Float64}
    scalex::SVector{D,Float64}

    n_total_conserved::Int64
    n_moments_vel::Int64
    rhs_vector::Vector{Float64}
    row_scale::Vector{Float64}
    mim::Vector{SVector{3,Int64}}  # mult-index moments
    tot_order::Vector{Int64}
    vel_powers::Matrix{Float64}
    n_moments_pos::Int64
    mim_pos::Vector{SVector{3,Int64}}  # mult-index moments
    tot_order_pos::Vector{Int64}
    pos_powers::Matrix{Float64}
    pos_i_x::Int64
    pos_i_y::Int64
    pos_i_z::Int64

    lhs_matrix_ncols_start::Int64
    lhs_matrix_ncols_end::Int64
    column_norms::Vector{Vector{Float64}}
    vel_pos_matrices::Vector{Matrix{Float64}}
    column_norms_scratch::Vector{Float64}
    vel_pos_matrix_scratch::Matrix{Float64}

    work::Vector{NNLSWorkspace{Float64, Int}}

    @doc """
        NNLSMerge{D}(multi_index_moments, init_np; rate_preserving=false, multi_index_moments_pos=[], matrix_ncol_nprealloc=0) where D

    Create NNLS-based merging. Mass, momentum, directional energy are always conserved:
    if not in the list `multi_index_moments` of moments to preserve, the corresponding
    moment multi-indices will be added automatically. These indices are `(0,0,0)` for mass,
    `(1,0,0)`, `(0,1,0)`, `(0,0,1)` for momentum, and `(2,0,0)`, `(0,2,0)`, `(0,0,2)`
    for directional energies. Spatial moments can be preserved by setting `multi_index_moments_pos`.
    If `multi_index_moments_pos` is non-empty, it is currently left to the user to include any relevant 1st order moments
    (corresponding center of mass conservation) i.e. `(1, 0, 0)`, `(0, 1, 0)`, `(0, 0, 1)`; otherwise
    the corresponding spatial moments including higher-order ones will not be conserved.
    By default, a single `NNLSWorkspace` is pre-allocated for the system.
    The number of columns in the LHS matrix is the number of particles to be merged (+ any fictitious particles);
    so it cannot be fixed in advance. By setting `matrix_ncol_nprealloc` to a value larger than 0,
    one can pre-allocate a range of workspaces with a fixed number of columns
    spanning `[init_np, init_np+matrix_ncol_nprealloc]`.
    In this case, `matrix_ncol_nprealloc+2` `NNLSWorkspace` instances are pre-allocated, with the last one being used
    in case the number of columns is not in the range (its LHS matrix is then re-allocated on the fly, but only
    when the number of columns differs from the previous call).
    Vectors `column_norms` and `vel_pos_matrices` are also then pre-allocated, holding respectively the inverses of the
    column-wise norms of the LHS matrices used for their scaling, and the particle velocities and positions.
    Outside the pre-allocated range the `column_norms_scratch` / `vel_pos_matrix_scratch` buffers are used instead;
    these grow as needed and are re-used across calls.

    # Positional arguments
    * `multi_index_moments`: vector of mixed moments to preserve of the form `[(i1, j1, k1), (i2, j2, k2), ...]``
    * `init_np`: assumption on pre-merge number of particles to pre-allocate memory for
    
    # Keyword arguments:
    * `rate_preserving`: used for rate-preserving merging of electrons, preserves approximate elastic collision and ionization rates
    * `multi_index_moments_pos`: list of spatial moments to preserve
    * `matrix_ncol_nprealloc`: number of NNLS workspaces with a fixed number of columns to pre-allocate

    # Throws
    * `ArgumentError`: if any of moment powers is negative
    """
    function NNLSMerge{D}(multi_index_moments, init_np; rate_preserving=false, multi_index_moments_pos=[], matrix_ncol_nprealloc=0) where D
        add_length = 0
        if rate_preserving
            add_length = 2
        end

        base_moments = base_multi_index_moments()
        filtered_mim = [x for x in multi_index_moments if !(x in base_moments)]

        mimpos = []

        if D == 3
            mimpos = copy(multi_index_moments_pos)
        elseif D == 2
            mimpos = [x for x in multi_index_moments_pos if (x[1] != 0) || (x[2] != 0)]  # filter out z-moments
        elseif D == 1
            mimpos = [x for x in multi_index_moments_pos if (x[1] != 0)]  # filter out y and z-moments
        end

        n_total_conserved = length(base_moments) + length(filtered_mim) + add_length + length(mimpos)
        append!(base_moments, filtered_mim)
        tot_order = [sum(mim) for mim in base_moments]
        tot_order_pos = [sum(mim) for mim in mimpos]

        pos_i_x = -1
        pos_i_y = -1
        pos_i_z = -1
        
        for i in 1:length(mimpos)
            if (mimpos[i][1] == 1) && (mimpos[i][2] == 0) && (mimpos[i][3] == 0)
                pos_i_x = i
            elseif (mimpos[i][1] == 0) && (mimpos[i][2] == 1) && (mimpos[i][3] == 0)
                pos_i_y = i
            elseif (mimpos[i][1] == 0) && (mimpos[i][2] == 0) && (mimpos[i][3] == 1)
                pos_i_z = i
            end
        end

        # moments are evaluated off tables of powers of the centered velocities/positions,
        # which only makes sense for non-negative multi-indices
        for m in Iterators.flatten((base_moments, mimpos))
            if minimum(m) < 0
                throw(ArgumentError("multi-index moment components must be non-negative, got $m"))
            end
        end

        mim_svec = SVector{3,Int64}[SVector{3,Int64}(m[1], m[2], m[3]) for m in base_moments]
        mim_pos_svec = SVector{3,Int64}[SVector{3,Int64}(m[1], m[2], m[3]) for m in mimpos]

        max_pow_vel = maximum(maximum(m) for m in mim_svec)
        max_pow_pos = length(mim_pos_svec) > 0 ? maximum(maximum(m) for m in mim_pos_svec) : 0

        column_norms = Vector{Vector{Float64}}([])
        nnls_ws_preallocated = Vector{NNLSWorkspace{Float64, Int}}([])
        vel_pos_matrices = Vector{Matrix{Float64}}([])
        if matrix_ncol_nprealloc > 0
            for i in init_np:init_np+matrix_ncol_nprealloc
                push!(column_norms, ones(i))
                push!(nnls_ws_preallocated, NNLSWorkspace(zeros(n_total_conserved, i), zeros(n_total_conserved)))
                push!(vel_pos_matrices, zeros(3+D, i))
            end

            # one extra workspace
            scratch_ncols = init_np + matrix_ncol_nprealloc + 1
            push!(nnls_ws_preallocated, NNLSWorkspace(zeros(n_total_conserved, scratch_ncols), zeros(n_total_conserved)))
        else
            scratch_ncols = init_np
            push!(nnls_ws_preallocated, NNLSWorkspace(zeros(n_total_conserved, scratch_ncols), zeros(n_total_conserved)))
        end

        return new{D}(SVector{3,Float64}(0.0, 0.0, 0.0), zero(SVector{D,Float64}),
                      1.0, 1.0, # vref, inv_vref
                      SVector{3,Float64}(0.0, 0.0, 0.0),   # std(v)
                      zero(SVector{D,Float64}),   # std(x)
                      0.0, # w_total
                      SVector{3,Float64}(0.0, 0.0, 0.0),   # scale v
                      zero(SVector{D,Float64}),   # scale x
                      n_total_conserved,
                      length(base_moments),
                      zeros(n_total_conserved),  # rhs_vector
                      zeros(n_total_conserved),  # row_scale
                      mim_svec, tot_order, zeros(3, max_pow_vel+1),
                      length(mim_pos_svec), mim_pos_svec, tot_order_pos, zeros(D, max_pow_pos+1),
                      pos_i_x, pos_i_y, pos_i_z,
                      init_np, init_np+matrix_ncol_nprealloc,
                      column_norms,
                      vel_pos_matrices,
                      ones(scratch_ncols),
                      zeros(3+D, scratch_ncols),
                      nnls_ws_preallocated)
    end

     @doc """
        NNLSMerge(multi_index_moments, init_np; rate_preserving=false, multi_index_moments_pos=[], matrix_ncol_nprealloc=0)

    Create NNLS-based merging for particles with a 3-dimensional position vector. Mass, momentum, directional energy are always conserved:
    if not in the list `multi_index_moments` of moments to preserve, the corresponding
    moment multi-indices will be added automatically. These indices are `(0,0,0)` for mass,
    `(1,0,0)`, `(0,1,0)`, `(0,0,1)` for momentum, and `(2,0,0)`, `(0,2,0)`, `(0,0,2)`
    for directional energies. Spatial moments can be preserved by setting `multi_index_moments_pos`.
    If `multi_index_moments_pos` is non-empty, it is currently left to the user to include any relevant 1st order moments
    (corresponding center of mass conservation) i.e. `(1, 0, 0)`, `(0, 1, 0)`, `(0, 0, 1)`; otherwise
    the corresponding spatial moments including higher-order ones will not be conserved.
    By default, a single `NNLSWorkspace` is pre-allocated for the system.
    The number of columns in the LHS matrix is the number of particles to be merged (+ any fictitious particles);
    so it cannot be fixed in advance. By setting `matrix_ncol_nprealloc` to a value larger than 0,
    one can pre-allocate a range of workspaces with a fixed number of columns
    spanning `[init_np, init_np+matrix_ncol_nprealloc]`.
    In this case, `matrix_ncol_nprealloc+2` `NNLSWorkspace` instances are pre-allocated, with the last one being used
    in case the number of columns is not in the range (its LHS matrix is then re-allocated on the fly, but only
    when the number of columns differs from the previous call).
    Vectors `column_norms` and `vel_pos_matrices` are also then pre-allocated, holding respectively the inverses of the
    column-wise norms of the LHS matrices used for their scaling, and the particle velocities and positions.
    Outside the pre-allocated range the `column_norms_scratch` / `vel_pos_matrix_scratch` buffers are used instead;
    these grow as needed and are re-used across calls.

    # Positional arguments
    * `multi_index_moments`: vector of mixed moments to preserve of the form `[(i1, j1, k1), (i2, j2, k2), ...]``
    * `init_np`: assumption on pre-merge number of particles to pre-allocate memory for
    
    # Keyword arguments:
    * `rate_preserving`: used for rate-preserving merging of electrons, preserves approximate elastic collision and ionization rates
    * `multi_index_moments_pos`: list of spatial moments to preserve
    * `matrix_ncol_nprealloc`: number of NNLS workspaces with a fixed number of columns to pre-allocate

    # Throws
    * `ArgumentError`: if any of moment powers is negative
    """
    function NNLSMerge(multi_index_moments, init_np; rate_preserving=false, multi_index_moments_pos=[], matrix_ncol_nprealloc=0)
        return NNLSMerge{3}(multi_index_moments, init_np; rate_preserving=rate_preserving, multi_index_moments_pos=multi_index_moments_pos, matrix_ncol_nprealloc=matrix_ncol_nprealloc)
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
    compute_w_total_v0!(nnls_merging, particles::ParticleVector{D}, pia, cell, species) where D

Compute total computational weight of particles and mean velocity, as well as velocity bounds of the
set of particles in each velocity direction.

# Positional arguments
* `nnls_merging`: the `NNLSMerge` instance where the computed values will be stored
* `particles`: the `ParticleVector` instance containing the particles that are being merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the grid cell in which particles are being merged
* `species`: the index of the species being merged
"""
function compute_w_total_v0!(nnls_merging, particles::ParticleVector{D}, pia, cell, species) where D
    nnls_merging.v0 = zero(SVector{3,Float64})
    nnls_merging.x0 = zero(SVector{D,Float64})
    nnls_merging.w_total = 0.0

    @inbounds s1 = pia.indexer[cell,species].start1
    @inbounds e1 = pia.indexer[cell,species].end1
    @inbounds for i in s1:e1
        p_i = particles[i]
        w = p_i.w
        v = p_i.v
        x = p_i.x

        nnls_merging.w_total += w
        nnls_merging.v0 = nnls_merging.v0 + w * v
        nnls_merging.x0 = nnls_merging.x0 + w * x
    end

    @inbounds if pia.indexer[cell,species].n_group2 > 0
        @inbounds s2 = pia.indexer[cell,species].start2
        @inbounds e2 = pia.indexer[cell,species].end2
        @inbounds for i in s2:e2
            p_i = particles[i]
            w = p_i.w
            v = p_i.v
            x = p_i.x

            nnls_merging.w_total += w
            nnls_merging.v0 = nnls_merging.v0 + w * v
            nnls_merging.x0 = nnls_merging.x0 + w * x
        end
    end

    nnls_merging.v0 = nnls_merging.v0 / nnls_merging.w_total
    nnls_merging.x0 = nnls_merging.x0 / nnls_merging.w_total
    return nothing
end

"""
    fill_powers!(powers, v, v0)

Fill a table of the powers of the components of the centered velocity / position vector `v - v0`,
so that `powers[j, e+1]` holds `(v[j] - v0[j])^e`. Only the first `size(powers, 1)` components
are considered, and powers up to `size(powers, 2) - 1` are computed.
Evaluating the moments off such a table avoids the runtime-exponent `^` calls that dominate
the cost of building the LHS matrix.

# Positional arguments
* `powers`: the table of powers to fill
* `v`: the velocity / position vector
* `v0`: the mean velocity / position vector
"""
@inline function fill_powers!(powers, v, v0)
    n_dim = size(powers, 1)
    n_pow = size(powers, 2)

    @inbounds for j in 1:n_dim
        powers[j, 1] = 1.0
    end

    @inbounds for e in 2:n_pow
        for j in 1:n_dim
            powers[j, e] = powers[j, e-1] * (v[j] - v0[j])
        end
    end
end

"""
    ccm_vel(vel_powers, mim)

Compute an unweighted central velocity moment from a table of powers filled by [`fill_powers!`](@ref).

# Positional arguments
* `vel_powers`: the table of powers of the centered velocity components
* `mim`: the 3-dimensional multi-index

# Returns
Computed unweighted central moment.
"""
@inline function ccm_vel(vel_powers, mim)
    @inbounds return vel_powers[1, mim[1]+1] * vel_powers[2, mim[2]+1] * vel_powers[3, mim[3]+1]
end

"""
    ccm_pos(pos_powers, mim, D)

Compute an unweighted central spatial moment from a table of powers filled by [`fill_powers!`](@ref).
Only the first `D` components of the multi-index are used, matching the dimensionality of the
particle position vectors.

# Positional arguments
* `pos_powers`: the table of powers of the centered position components
* `mim`: the 3-dimensional multi-index
* `D`: the dimensionality of the position vectors

# Returns
Computed unweighted central moment.
"""
@inline function ccm_pos(pos_powers, mim, D)
    res = 1.0
    @inbounds for j in 1:D
        res *= pos_powers[j, mim[j]+1]
    end
    return res
end


"""
    compute_lhs_and_rhs!(nnls_merging::NNLSMerge{D}, lhs_matrix, vel_pos_matrix, particles::ParticleVector{D}, pia, cell, species) where D

Compute LHS matrix and RHS vector for NNLS merging. Returns the pre-merge number of particles.

# Positional arguments
* `nnls_merging`: the `NNLSMerge` instance where the RHS vector will be stored
* `lhs_matrix`: the matrix of size `n_total_conserved x n_particles`, where
    `n_total_conserved` is the number of conserved moments and `n_particles` is the pre-merge number of particles
    + any fictitious particles
* `vel_pos_matrix`: the matrix of size `6 x n_particles`, where `n_particles` is the pre-merge number of particles
    + any fictitious particles, where their velocities and positions will be stored
* `particles`: the `ParticleVector` instance containing the particles that are being merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the grid cell in which particles are being merged
* `species`: the index of the species being merged

# Returns
The pre-merge number of particles.
"""
function compute_lhs_and_rhs!(nnls_merging::NNLSMerge{D}, lhs_matrix, vel_pos_matrix,
                              particles::ParticleVector{D}, pia, cell, species) where D
    n_moms = nnls_merging.n_moments_vel
    n_moms_pos = nnls_merging.n_moments_pos

    fill!(nnls_merging.rhs_vector, 0.0)

    compute_w_total_v0!(nnls_merging, particles, pia, cell, species)

    rhs_vector = nnls_merging.rhs_vector
    mim = nnls_merging.mim
    mim_pos = nnls_merging.mim_pos
    vel_powers = nnls_merging.vel_powers
    pos_powers = nnls_merging.pos_powers
    v0 = nnls_merging.v0
    x0 = nnls_merging.x0

    Ev = zero(SVector{3,Float64})
    Ex = zero(SVector{D,Float64})

    col_index = 1
    @inbounds s1 = pia.indexer[cell,species].start1
    @inbounds e1 = pia.indexer[cell,species].end1

    @inbounds for i in s1:e1
        p_i = particles[i]
        w = p_i.w
        v = p_i.v
        x = p_i.x

        Ev = Ev + w * (v - v0).^2
        Ex = Ex + w * (x - x0).^2

        vel_pos_matrix[1, col_index] = v[1]
        vel_pos_matrix[2, col_index] = v[2]
        vel_pos_matrix[3, col_index] = v[3]

        for j in 1:D
            vel_pos_matrix[3+j, col_index] = x[j]
        end

        fill_powers!(vel_powers, v, v0)
        for n_mom in 1:n_moms
            tmp_ccm = ccm_vel(vel_powers, mim[n_mom])
            rhs_vector[n_mom] = rhs_vector[n_mom] + w * tmp_ccm
            lhs_matrix[n_mom, col_index] = tmp_ccm
        end

        if n_moms_pos > 0
            fill_powers!(pos_powers, x, x0)
            for n_mom in 1:n_moms_pos
                tmp_ccm = ccm_pos(pos_powers, mim_pos[n_mom], D)
                rhs_vector[n_moms+n_mom] = rhs_vector[n_moms+n_mom] + w * tmp_ccm
                lhs_matrix[n_moms+n_mom, col_index] = tmp_ccm
            end
        end
        col_index += 1
    end

    @inbounds if pia.indexer[cell,species].n_group2 > 0

        @inbounds s2 = pia.indexer[cell,species].start2
        @inbounds e2 = pia.indexer[cell,species].end2

        @inbounds for i in s2:e2
            p_i = particles[i]
            w = p_i.w
            v = p_i.v
            x = p_i.x

            Ev = Ev + w * (v - v0).^2
            Ex = Ex + w * (x - x0).^2

            vel_pos_matrix[1, col_index] = v[1]
            vel_pos_matrix[2, col_index] = v[2]
            vel_pos_matrix[3, col_index] = v[3]

            for j in 1:D
                vel_pos_matrix[3+j, col_index] = x[j]
            end

            fill_powers!(vel_powers, v, v0)
            for n_mom in 1:n_moms
                tmp_ccm = ccm_vel(vel_powers, mim[n_mom])
                rhs_vector[n_mom] = rhs_vector[n_mom] + w * tmp_ccm
                lhs_matrix[n_mom, col_index] = tmp_ccm
            end

            if n_moms_pos > 0
                fill_powers!(pos_powers, x, x0)
                for n_mom in 1:n_moms_pos
                    tmp_ccm = ccm_pos(pos_powers, mim_pos[n_mom], D)
                    rhs_vector[n_moms+n_mom] = rhs_vector[n_moms+n_mom] + w * tmp_ccm
                    lhs_matrix[n_moms+n_mom, col_index] = tmp_ccm
                end
            end
            col_index += 1
        end
    end

    n_total_conserved = nnls_merging.n_total_conserved
    w_tot = nnls_merging.w_total

    @inbounds @simd for i in 1:n_total_conserved
        rhs_vector[i] /= w_tot
    end

    nnls_merging.Ev = sqrt.(Ev / w_tot)
    nnls_merging.Ex = sqrt.(Ex / w_tot)

    return col_index
end

"""
    compute_lhs_and_rhs_rate_preserving!(nnls_merging::NNLSMerge{D}, lhs_matrix, vel_pos_matrix,
                                         interaction, electron_neutral_interactions, computed_cs,
                                         particles::ParticleVector{D}, pia, cell, species, neutral_species_index, extend) where D

Compute LHS matrix and RHS vector for the rate-preserving NNLS merging (for electrons). Approximate
    elastic scattering and electron-impact ionization rates are conserved.

# Positional arguments
* `nnls_merging`: the `NNLSMerge` instance where the RHS vector will be stored
* `lhs_matrix`: the matrix of size `n_total_conserved x n_particles`, where
    `n_total_conserved` is the number of conserved moments and `n_particles` is the pre-merge number of particles
    + any fictitious particles
* `vel_pos_matrix`: the matrix of size `6 x n_particles`, where `n_particles` is the pre-merge number of particles
    + any fictitious particles, where their velocities and positions will be stored
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
function compute_lhs_and_rhs_rate_preserving!(nnls_merging::NNLSMerge{D}, lhs_matrix, vel_pos_matrix,
                                              interaction, electron_neutral_interactions, computed_cs,
                                              particles::ParticleVector{D}, pia, cell, species, neutral_species_index, extend) where D
    n_moms = nnls_merging.n_moments_vel

    fill!(nnls_merging.rhs_vector, 0.0)

    compute_w_total_v0!(nnls_merging, particles, pia, cell, species)

    rhs_vector = nnls_merging.rhs_vector
    mim = nnls_merging.mim
    vel_powers = nnls_merging.vel_powers
    v0 = nnls_merging.v0
    x0 = nnls_merging.x0

    Ev = zero(SVector{3,Float64})
    Ex = zero(SVector{D,Float64})

    col_index = 1
    @inbounds s1 = pia.indexer[cell,species].start1
    @inbounds e1 = pia.indexer[cell,species].end1
    @inbounds for i in s1:e1
        # w_total += particles[i].w
        p_i = particles[i]
        w = p_i.w
        v = p_i.v
        x = p_i.x

        Ev = Ev + w * (v - v0).^2
        Ex = Ex + w * (x - x0).^2

        vel_pos_matrix[1, col_index] = v[1]
        vel_pos_matrix[2, col_index] = v[2]
        vel_pos_matrix[3, col_index] = v[3]

        for j in 1:D
            vel_pos_matrix[3+j, col_index] = x[j]
        end

        fill_powers!(vel_powers, v, v0)
        for n_mom in 1:n_moms
            tmp_ccm = ccm_vel(vel_powers, mim[n_mom])
            rhs_vector[n_mom] = rhs_vector[n_mom] + w * tmp_ccm
            lhs_matrix[n_mom, col_index] = tmp_ccm
        end

        g = norm(v)
        compute_cross_sections_only!(computed_cs, interaction, g, electron_neutral_interactions, neutral_species_index, extend)
        cse_g = get_cs_elastic(electron_neutral_interactions, computed_cs, neutral_species_index) * g
        csi_g = get_cs_ionization(electron_neutral_interactions, computed_cs, neutral_species_index) * g
        lhs_matrix[n_moms+1, col_index] = cse_g
        lhs_matrix[n_moms+2, col_index] = csi_g
        rhs_vector[n_moms+1] = rhs_vector[n_moms+1] + cse_g * w
        rhs_vector[n_moms+2] = rhs_vector[n_moms+2] + csi_g * w

        col_index += 1
    end

    @inbounds if pia.indexer[cell,species].n_group2 > 0
        @inbounds s2 = pia.indexer[cell,species].start2
        @inbounds e2 = pia.indexer[cell,species].end2

        @inbounds for i in s2:e2
            p_i = particles[i]
            w = p_i.w
            v = p_i.v
            x = p_i.x

            Ev = Ev + w * (v - v0).^2
            Ex = Ex + w * (x - x0).^2

            vel_pos_matrix[1, col_index] = v[1]
            vel_pos_matrix[2, col_index] = v[2]
            vel_pos_matrix[3, col_index] = v[3]

            for j in 1:D
                vel_pos_matrix[3+j, col_index] = x[j]
            end

            # w_total += particles[i].w
            fill_powers!(vel_powers, v, v0)
            for n_mom in 1:n_moms
                tmp_ccm = ccm_vel(vel_powers, mim[n_mom])
                rhs_vector[n_mom] = rhs_vector[n_mom] + w * tmp_ccm
                lhs_matrix[n_mom, col_index] = tmp_ccm
            end

            g = norm(v)
            compute_cross_sections_only!(computed_cs, interaction, g, electron_neutral_interactions, neutral_species_index, extend)
            cse_g = get_cs_elastic(electron_neutral_interactions, computed_cs, neutral_species_index) * g
            csi_g = get_cs_ionization(electron_neutral_interactions, computed_cs, neutral_species_index) * g
            lhs_matrix[n_moms+1, col_index] = cse_g
            lhs_matrix[n_moms+2, col_index] = csi_g
            rhs_vector[n_moms+1] = rhs_vector[n_moms+1] + cse_g * w
            rhs_vector[n_moms+2] = rhs_vector[n_moms+2] + csi_g * w

            col_index += 1
        end
    end

    n_total_conserved = nnls_merging.n_total_conserved
    w_tot = nnls_merging.w_total

    @inbounds @simd for i in 1:n_total_conserved
        rhs_vector[i] /= w_tot
    end

    nnls_merging.Ev = sqrt.(Ev / w_tot)
    nnls_merging.Ex = sqrt.(Ex / w_tot)

    return col_index
end

"""
    compute_lhs_and_rhs_rate_preserving!(nnls_merging::NNLSMerge{D}, lhs_matrix, vel_pos_matrix,
                                         interaction, electron_neutral_interactions, computed_cs,
                                         particles::ParticleVector{D}, particles_neutral::ParticleVector{D}, pia, cell, species, neutral_species_index, extend) where D

Compute LHS matrix and RHS vector for the rate-preserving NNLS merging (for electrons). Exact
    elastic scattering and electron-impact ionization rates are conserved.

# Positional arguments
* `nnls_merging`: the `NNLSMerge` instance where the RHS vector will be stored
* `lhs_matrix`: the matrix of size `n_total_conserved x n_particles`, where
    `n_total_conserved` is the number of conserved moments and `n_particles` is the pre-merge number of particles
    + any fictitious particles
* `vel_pos_matrix`: the matrix of size `6 x n_particles`, where `n_particles` is the pre-merge number of particles
    + any fictitious particles, where their velocities and positions will be stored
* `interaction`: the `Interaction` instance describing the electron-neutral interaction being considered
* `electron_neutral_interactions`:  the `ElectronNeutralInteractions` instance storing the tabulated cross-section
    data used to compute the rates
* `computed_cs`: the vector of `ComputedCrossSection` instances in which the computed values will be stored
* `particles`: the `ParticleVector` instance containing the particles that are being merged
* `particles_neutrals`: the `ParticleVector` instance containing the neutral collision partner particles
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the grid cell in which particles are being merged
* `species`: the index of the species being merged
* `neutral_species_index`: the index of the neutral species which is the collision partner in the electron-neutral
    collisions for which approximate rates are being preserved.
* `extend`: enum of `CSExtend` type that sets how out-of-range energy values are treated when computing cross-sections
"""
function compute_lhs_and_rhs_rate_preserving!(nnls_merging::NNLSMerge{D}, lhs_matrix, vel_pos_matrix,
                                              interaction, electron_neutral_interactions, computed_cs,
                                              particles::ParticleVector{D}, particles_neutral::ParticleVector{D}, pia, cell, species, neutral_species_index, extend) where D
    n_moms = nnls_merging.n_moments_vel

    fill!(nnls_merging.rhs_vector, 0.0)

    compute_w_total_v0!(nnls_merging, particles, pia, cell, species)

    rhs_vector = nnls_merging.rhs_vector
    mim = nnls_merging.mim
    vel_powers = nnls_merging.vel_powers
    v0 = nnls_merging.v0
    x0 = nnls_merging.x0

    Ev = zero(SVector{3,Float64})
    Ex = zero(SVector{D,Float64})

    col_index = 1
    @inbounds s1 = pia.indexer[cell,species].start1
    @inbounds e1 = pia.indexer[cell,species].end1
    @inbounds k_s1 = pia.indexer[cell,neutral_species_index].start1
    @inbounds k_e1 = pia.indexer[cell,neutral_species_index].end1
    @inbounds k_s2 = pia.indexer[cell,neutral_species_index].start2
    @inbounds k_e2 = pia.indexer[cell,neutral_species_index].end2
    @inbounds for i in s1:e1
        p_i = particles[i]
        w = p_i.w
        v = p_i.v
        x = p_i.x

        Ev = Ev + w * (v - v0).^2
        Ex = Ex + w * (x - x0).^2

        vel_pos_matrix[1, col_index] = v[1]
        vel_pos_matrix[2, col_index] = v[2]
        vel_pos_matrix[3, col_index] = v[3]

        for j in 1:D
            vel_pos_matrix[3+j, col_index] = x[j]
        end

        fill_powers!(vel_powers, v, v0)
        for n_mom in 1:n_moms
            tmp_ccm = ccm_vel(vel_powers, mim[n_mom])
            rhs_vector[n_mom] = rhs_vector[n_mom] + w * tmp_ccm
            lhs_matrix[n_mom, col_index] = tmp_ccm
        end

        cse_g = 0.0
        csi_g = 0.0
        w_k = 0.0
        @inbounds for k in k_s1:k_e1
            p_k = particles_neutral[k]
            p_k_w = p_k.w

            g = norm(v - p_k.v)
            compute_cross_sections_only!(computed_cs, interaction, g, electron_neutral_interactions, neutral_species_index, extend)
            cse_g = cse_g + get_cs_elastic(electron_neutral_interactions, computed_cs, neutral_species_index) * g * p_k_w
            csi_g = csi_g + get_cs_ionization(electron_neutral_interactions, computed_cs, neutral_species_index) * g * p_k_w
            w_k = w_k + p_k_w
        end

        @inbounds for k in k_s2:k_e2
            p_k = particles_neutral[k]
            p_k_w = p_k.w

            g = norm(v - p_k.v)
            compute_cross_sections_only!(computed_cs, interaction, g, electron_neutral_interactions, neutral_species_index, extend)
            cse_g = cse_g + get_cs_elastic(electron_neutral_interactions, computed_cs, neutral_species_index) * g * p_k_w
            csi_g = csi_g + get_cs_ionization(electron_neutral_interactions, computed_cs, neutral_species_index) * g * p_k_w
            w_k = w_k + p_k_w
        end

        if w_k > 0
            cse_g = cse_g / w_k
            csi_g = csi_g / w_k
        end

        lhs_matrix[n_moms+1, col_index] = cse_g
        lhs_matrix[n_moms+2, col_index] = csi_g
        rhs_vector[n_moms+1] = rhs_vector[n_moms+1] + cse_g * w
        rhs_vector[n_moms+2] = rhs_vector[n_moms+2] + csi_g * w

        col_index += 1
    end

    @inbounds if pia.indexer[cell,species].n_group2 > 0
        @inbounds s2 = pia.indexer[cell,species].start2
        @inbounds e2 = pia.indexer[cell,species].end2

        @inbounds for i in s2:e2
            p_i = particles[i]
            w = p_i.w
            v = p_i.v
            x = p_i.x

            Ev = Ev + w * (v - v0).^2
            Ex = Ex + w * (x - x0).^2

            vel_pos_matrix[1, col_index] = v[1]
            vel_pos_matrix[2, col_index] = v[2]
            vel_pos_matrix[3, col_index] = v[3]

            for j in 1:D
                vel_pos_matrix[3+j, col_index] = x[j]
            end

            fill_powers!(vel_powers, v, v0)
            for n_mom in 1:n_moms
                tmp_ccm = ccm_vel(vel_powers, mim[n_mom])
                rhs_vector[n_mom] = rhs_vector[n_mom] + w * tmp_ccm
                lhs_matrix[n_mom, col_index] = tmp_ccm
            end

            cse_g = 0.0
            csi_g = 0.0
            w_k = 0.0
            @inbounds for k in k_s1:k_e1
                p_k = particles_neutral[k]
                p_k_w = p_k.w

                g = norm(v - p_k.v)
                compute_cross_sections_only!(computed_cs, interaction, g, electron_neutral_interactions, neutral_species_index, extend)
                cse_g = cse_g + get_cs_elastic(electron_neutral_interactions, computed_cs, neutral_species_index) * g * p_k_w
                csi_g = csi_g + get_cs_ionization(electron_neutral_interactions, computed_cs, neutral_species_index) * g * p_k_w
                w_k = w_k + p_k_w
            end
            @inbounds for k in k_s2:k_e2
                p_k = particles_neutral[k]
                p_k_w = p_k.w

                g = norm(v - p_k.v)
                compute_cross_sections_only!(computed_cs, interaction, g, electron_neutral_interactions, neutral_species_index, extend)
                cse_g = cse_g + get_cs_elastic(electron_neutral_interactions, computed_cs, neutral_species_index) * g * p_k_w
                csi_g = csi_g + get_cs_ionization(electron_neutral_interactions, computed_cs, neutral_species_index) * g * p_k_w
                w_k = w_k + p_k_w
            end

            if w_k > 0
                cse_g = cse_g / w_k
                csi_g = csi_g / w_k
            end

            lhs_matrix[n_moms+1, col_index] = cse_g
            lhs_matrix[n_moms+2, col_index] = csi_g
            rhs_vector[n_moms+1] = rhs_vector[n_moms+1] + cse_g * w
            rhs_vector[n_moms+2] = rhs_vector[n_moms+2] + csi_g * w

            col_index += 1
        end
    end

    n_total_conserved = nnls_merging.n_total_conserved
    w_tot = nnls_merging.w_total

    @inbounds @simd for i in 1:n_total_conserved
        rhs_vector[i] /= w_tot
    end

    nnls_merging.Ev = sqrt.(Ev / w_tot)
    nnls_merging.Ex = sqrt.(Ex / w_tot)

    return col_index
end

"""
    scale_lhs_rhs_vref!(nnls_merging::NNLSMerge{D}, lhs_matrix, lhs_ncols) where D

Scale the LHS and RHS of the NNLS system using the reference velocity ``v_{ref}``. Each moment is scaled
by ``(1/v_{ref})^{n_{tot}}``, where ``n_{tot}`` is the total order of the moment (i.e. for
a moment with multi-index `(i,j,k)` the total order is `i+j+k`).

# Positional arguments
* `nnls_merging`: the `NNLSMerge` instance
* `lhs_matrix`: the matrix of the LHS
* `lhs_ncols`: number of columns in the LHS matrix
"""
function scale_lhs_rhs_vref!(nnls_merging::NNLSMerge{D}, lhs_matrix, lhs_ncols) where D
    n_moms = nnls_merging.n_moments_vel
    row_scale = nnls_merging.row_scale
    @inbounds for n_mom in 1:n_moms
        ref_val = nnls_merging.inv_vref^nnls_merging.tot_order[n_mom]
        row_scale[n_mom] = ref_val
        nnls_merging.rhs_vector[n_mom] *= ref_val
    end
    @inbounds for col in 1:lhs_ncols
        @simd for n_mom in 1:n_moms
            lhs_matrix[n_mom, col] *= row_scale[n_mom]
        end
    end
    nnls_merging.scalev = SVector{3,Float64}(nnls_merging.vref, nnls_merging.vref, nnls_merging.vref)
    return nothing
end

"""
    scale_lhs_rhs_variance!(nnls_merging::NNLSMerge{D}, lhs_matrix, lhs_ncols) where D

Scale the LHS and RHS of the NNLS system using the computed variances of the particles velocities
in the ``x``, ``y``, ``z`` directions. If any of the variances is smaller than 1e-6, then
``v_{ref}`` is used.
Each moment with multi-index `(i,j,k)` is scaled
by ``(1/Ev[1])^{i}(1/Ev[2])^{j}(1/Ev[3])^{k}``.

# Positional arguments
* `nnls_merging`: the `NNLSMerge` instance
* `lhs_matrix`: the matrix of the LHS
* `lhs_ncols`: number of columns in the LHS matrix
"""
function scale_lhs_rhs_variance!(nnls_merging::NNLSMerge{D}, lhs_matrix, lhs_ncols) where D
    n_moms = nnls_merging.n_moments_vel
    row_scale = nnls_merging.row_scale
    inv_ev = ifelse.(nnls_merging.Ev .> 1e-6, 1.0 ./ nnls_merging.Ev, nnls_merging.inv_vref)
    @inbounds for n_mom in 1:n_moms
        ref_val = (inv_ev[1]^nnls_merging.mim[n_mom][1]) * (inv_ev[2]^nnls_merging.mim[n_mom][2]) * (inv_ev[3]^nnls_merging.mim[n_mom][3])
        row_scale[n_mom] = ref_val
        nnls_merging.rhs_vector[n_mom] *= ref_val
    end
    @inbounds for col in 1:lhs_ncols
        @simd for n_mom in 1:n_moms
            lhs_matrix[n_mom, col] *= row_scale[n_mom]
        end
    end
    nnls_merging.scalev = 1.0 ./ inv_ev
    return nothing
end

"""
    scale_lhs_rhs_spatial_variance!(nnls_merging::NNLSMerge{D}, lhs_matrix, lhs_ncols) where D

Scale the spatial moments in the LHS and RHS of the NNLS system using the computed variances of the particles positions
in the ``x``, ``y``, ``z`` directions. If any of the variances is smaller than 1e-6, then
no scaling is done.
Each spatial moment with multi-index `(i,j,k)` is scaled
by ``(1/Ex[1])^{i}(1/Ex[2])^{j}(1/Ex[3])^{k}``, for non-3-dimensional particle vectors only the first D components are considered and scaled.

# Positional arguments
* `nnls_merging`: the `NNLSMerge` instance
* `lhs_matrix`: the matrix of the LHS
* `lhs_ncols`: number of columns in the LHS matrix
"""
function scale_lhs_rhs_spatial_variance!(nnls_merging::NNLSMerge{D}, lhs_matrix, lhs_ncols) where D
    n_moms = nnls_merging.n_moments_vel
    n_moms_pos = nnls_merging.n_moments_pos
    row_scale = nnls_merging.row_scale

    inv_ex = ifelse.(nnls_merging.Ex .> 1e-6, 1.0 ./ nnls_merging.Ex, 1.0)
    @inbounds for n_mom in 1:n_moms_pos
        ref_val = 1.0
        for j in 1:D
            ref_val *= inv_ex[j]^nnls_merging.mim_pos[n_mom][j]
        end
        row_scale[n_moms + n_mom] = ref_val
        nnls_merging.rhs_vector[n_moms + n_mom] *= ref_val
    end
    @inbounds for col in 1:lhs_ncols
        @simd for n_mom in 1:n_moms_pos
            lhs_matrix[n_moms + n_mom, col] *= row_scale[n_moms + n_mom]
        end
    end
    nnls_merging.scalex = 1.0 ./ inv_ex
    return nothing
end

"""
    scale_lhs_rhs!(nnls_merging::NNLSMerge{D}, lhs_matrix, scaling, lhs_ncols) where D

Scale the LHS and RHS of the NNLS system using either the reference velocity or the computed variances of the velocity
in each direction.

# Positional arguments
* `nnls_merging`: the `NNLSMerge` instance
* `lhs_matrix`: the matrix of the LHS
* `scaling`: how to scale entries (`:vref` or `:variance`)
* `lhs_ncols`: number of columns in the LHS matrix
"""
function scale_lhs_rhs!(nnls_merging::NNLSMerge{D}, lhs_matrix, scaling, lhs_ncols) where D
    if scaling == :variance
        scale_lhs_rhs_variance!(nnls_merging, lhs_matrix, lhs_ncols)
        scale_lhs_rhs_spatial_variance!(nnls_merging, lhs_matrix, lhs_ncols)
    elseif scaling == :vref
        scale_lhs_rhs_vref!(nnls_merging, lhs_matrix, lhs_ncols)
        scale_lhs_rhs_spatial_variance!(nnls_merging, lhs_matrix, lhs_ncols)
    end
    return nothing
end

"""
    scale_lhs_rhs_rate_preserving!(nnls_merging, lhs_matrix, ref_k_elastic, ref_k_ion, scaling, lhs_ncols) where D

Scale the LHS and RHS of the NNLS system for the rate-preserving electron merging
using the reference velocity ``v_{ref}`` and
reference elastic scattering and ionization cross-sections. Each moment is scaled
by ``(1/v_{ref})^{n_{tot}}``, where ``n_{tot}`` is the total order of the moment (i.e. for
a moment with multi-index `(i,j,k)` the total order is `i+j+k`).
The entries in the LHS and RHS corresponding to the rates are scaled by
`ref_k_elastic`, `ref_k_ion`. 

# Positional arguments
* `nnls_merging`: the `NNLSMerge` instance
* `lhs_matrix`: the matrix of the LHS
* `ref_k_elastic`: the reference elastic collision rate coefficient
* `ref_k_ion`: the reference ionization rate coefficient
* `scaling`: how to scale entries (`:vref` or `:variance`)
* `lhs_ncols`: number of columns in the LHS matrix
"""
function scale_lhs_rhs_rate_preserving!(nnls_merging::NNLSMerge{D}, lhs_matrix, ref_k_elastic, ref_k_ion, scaling, lhs_ncols) where D
    scale_lhs_rhs!(nnls_merging, lhs_matrix, scaling, lhs_ncols)

    scaler_el = 1.0/ref_k_elastic
    scaler_ion = 1.0/ref_k_ion

    @inbounds for col in 1:lhs_ncols
        lhs_matrix[nnls_merging.n_moments_vel+1, col] *= scaler_el
        lhs_matrix[nnls_merging.n_moments_vel+2, col] *= scaler_ion
    end

    @inbounds nnls_merging.rhs_vector[nnls_merging.n_moments_vel+1] *= scaler_el
    @inbounds nnls_merging.rhs_vector[nnls_merging.n_moments_vel+2] *= scaler_ion
    return nothing
end

"""
    compute_post_merge_particles_nnls!(nnls_merging::NNLSMerge{3}, x::Vector{Float64}, particles::ParticleVector{3},
                                       pia, cell, species, lhs_ncols,
                                       vel_pos_matrix,
                                       max_err, w_threshold, work_index, column_norms)

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
* `vel_pos_matrix`: the matrix of size `6 x n_particles`, where `n_particles` is the pre-merge number of particles
    + any fictitious particles, where their velocities and positions will be stored
* `max_err`: maximum allowed value of the residual of the NNLS system
* `w_threshold`: the relative (w.r.t the total computational weight of the particles being merge) value of the computational weight below which particles are discarded
* `work_index`: the index of the `NNLSWorkspace` used to solve the NNLS system
* `column_norms`: the vector of the column-wise norms of the LHS matrix
"""
function compute_post_merge_particles_nnls!(nnls_merging::NNLSMerge{3}, x::Vector{Float64}, particles::ParticleVector{3},
                                            pia, cell, species, lhs_ncols,
                                            vel_pos_matrix,
                                            max_err, w_threshold, work_index, column_norms)
    @inbounds if nnls_merging.work[work_index].rnorm > max_err
        return -1
    end

    @inbounds indexer = pia.indexer[cell, species]

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

    @inbounds if nonzero >= indexer.n_local
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
            i = map_cont_index(indexer, curr_particle_index)
            curr_particle_index += 1

            p_i = particles[i]

            p_i.w = x[j] * nnls_merging.w_total# * column_norms[j]
            p_i.v = SVector{3, Float64}(vel_pos_matrix[1,j],
                                        vel_pos_matrix[2,j],
                                        vel_pos_matrix[3,j])
            
            # write positions or mean position if not specified as conserved momentt
            px = nnls_merging.pos_i_x > 0.0 ? vel_pos_matrix[4,j] : nnls_merging.x0[1]
            py = nnls_merging.pos_i_y > 0.0 ? vel_pos_matrix[5,j] : nnls_merging.x0[2]
            pz = nnls_merging.pos_i_z > 0.0 ? vel_pos_matrix[6,j] : nnls_merging.x0[3]

            p_i.x = SVector{3,Float64}(px, py, pz)
        end
    end

    old_count = indexer.n_local
    n_particles_to_delete = old_count - curr_particle_index

    # if we delete from particles in last cell AND we delete less particles than were in group 2
    # then continuity is not broken
    # !(A && B) == !A || !B
    @inbounds if !(cell == pia.n_cells) || (n_particles_to_delete > indexer.n_group2)
        pia.contiguous[species] = false
    end

    for _ in 1:n_particles_to_delete
        delete_particle_end!(particles, pia, cell, species)
    end

    return 1
end

"""
    compute_post_merge_particles_nnls!(nnls_merging::NNLSMerge{2}, x::Vector{Float64}, particles::ParticleVector{2},
                                       pia, cell, species, lhs_ncols,
                                       vel_pos_matrix,
                                       max_err, w_threshold, work_index, column_norms)

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
* `vel_pos_matrix`: the matrix of size `6 x n_particles`, where `n_particles` is the pre-merge number of particles
    + any fictitious particles, where their velocities and positions will be stored
* `max_err`: maximum allowed value of the residual of the NNLS system
* `w_threshold`: the relative (w.r.t the total computational weight of the particles being merge) value of the computational weight below which particles are discarded
* `work_index`: the index of the `NNLSWorkspace` used to solve the NNLS system
* `column_norms`: the vector of the column-wise norms of the LHS matrix
"""
function compute_post_merge_particles_nnls!(nnls_merging::NNLSMerge{2}, x::Vector{Float64}, particles::ParticleVector{2},
                                            pia, cell, species, lhs_ncols,
                                            vel_pos_matrix,
                                            max_err, w_threshold, work_index, column_norms)
    @inbounds if nnls_merging.work[work_index].rnorm > max_err
        return -1
    end

    @inbounds indexer = pia.indexer[cell, species]

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

    @inbounds if nonzero >= indexer.n_local
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
            i = map_cont_index(indexer, curr_particle_index)
            curr_particle_index += 1

            p_i = particles[i]

            p_i.w = x[j] * nnls_merging.w_total# * column_norms[j]
            p_i.v = SVector{3, Float64}(vel_pos_matrix[1,j],
                                        vel_pos_matrix[2,j],
                                        vel_pos_matrix[3,j])
            
            # write positions or mean position if not specified as conserved momentt
            px = nnls_merging.pos_i_x > 0.0 ? vel_pos_matrix[4,j] : nnls_merging.x0[1]
            py = nnls_merging.pos_i_y > 0.0 ? vel_pos_matrix[5,j] : nnls_merging.x0[2]

            p_i.x = SVector{2,Float64}(px, py)
        end
    end

    old_count = indexer.n_local
    n_particles_to_delete = old_count - curr_particle_index

    # if we delete from particles in last cell AND we delete less particles than were in group 2
    # then continuity is not broken
    # !(A && B) == !A || !B
    @inbounds if !(cell == pia.n_cells) || (n_particles_to_delete > indexer.n_group2)
        pia.contiguous[species] = false
    end

    for _ in 1:n_particles_to_delete
        delete_particle_end!(particles, pia, cell, species)
    end

    return 1
end

"""
    compute_post_merge_particles_nnls!(nnls_merging::NNLSMerge{1}, x::Vector{Float64}, particles::ParticleVector{1},
                                       pia, cell, species, lhs_ncols,
                                       vel_pos_matrix,
                                       max_err, w_threshold, work_index, column_norms)

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
* `vel_pos_matrix`: the matrix of size `6 x n_particles`, where `n_particles` is the pre-merge number of particles
    + any fictitious particles, where their velocities and positions will be stored
* `max_err`: maximum allowed value of the residual of the NNLS system
* `w_threshold`: the relative (w.r.t the total computational weight of the particles being merge) value of the computational weight below which particles are discarded
* `work_index`: the index of the `NNLSWorkspace` used to solve the NNLS system
* `column_norms`: the vector of the column-wise norms of the LHS matrix
"""
function compute_post_merge_particles_nnls!(nnls_merging::NNLSMerge{1}, x::Vector{Float64}, particles::ParticleVector{1},
                                            pia, cell, species, lhs_ncols,
                                            vel_pos_matrix,
                                            max_err, w_threshold, work_index, column_norms)
    @inbounds if nnls_merging.work[work_index].rnorm > max_err
        return -1
    end

    @inbounds indexer = pia.indexer[cell, species]

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

    @inbounds if nonzero >= indexer.n_local
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
            i = map_cont_index(indexer, curr_particle_index)
            curr_particle_index += 1

            p_i = particles[i]

            p_i.w = x[j] * nnls_merging.w_total# * column_norms[j]
            p_i.v = SVector{3, Float64}(vel_pos_matrix[1,j],
                                        vel_pos_matrix[2,j],
                                        vel_pos_matrix[3,j])
            
            # write positions or mean position if not specified as conserved momentt
            px = nnls_merging.pos_i_x > 0.0 ? vel_pos_matrix[4,j] : nnls_merging.x0[1]

            p_i.x = SVector{1,Float64}(px)
        end
    end

    old_count = indexer.n_local
    n_particles_to_delete = old_count - curr_particle_index

    # if we delete from particles in last cell AND we delete less particles than were in group 2
    # then continuity is not broken
    # !(A && B) == !A || !B
    @inbounds if !(cell == pia.n_cells) || (n_particles_to_delete > indexer.n_group2)
        pia.contiguous[species] = false
    end

    for _ in 1:n_particles_to_delete
        delete_particle_end!(particles, pia, cell, species)
    end

    return 1
end

"""
    compute_post_merge_particles_nnls!(nnls_merging::NNLSMerge{1}, x::Vector{Float64}, particles::ParticleVector{1},
                                       pia, cell, species, lhs_ncols,
                                       vel_pos_matrix,
                                       max_err, w_threshold, work_index, column_norms)

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
* `vel_pos_matrix`: the matrix of size `6 x n_particles`, where `n_particles` is the pre-merge number of particles
    + any fictitious particles, where their velocities and positions will be stored
* `max_err`: maximum allowed value of the residual of the NNLS system
* `w_threshold`: the relative (w.r.t the total computational weight of the particles being merge) value of the computational weight below which particles are discarded
* `work_index`: the index of the `NNLSWorkspace` used to solve the NNLS system
* `column_norms`: the vector of the column-wise norms of the LHS matrix
"""
function compute_post_merge_particles_nnls!(nnls_merging::NNLSMerge{0}, x::Vector{Float64}, particles::ParticleVector{0},
                                            pia, cell, species, lhs_ncols,
                                            vel_pos_matrix,
                                            max_err, w_threshold, work_index, column_norms)
    @inbounds if nnls_merging.work[work_index].rnorm > max_err
        return -1
    end

    @inbounds indexer = pia.indexer[cell, species]

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

    @inbounds if nonzero >= indexer.n_local
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
            i = map_cont_index(indexer, curr_particle_index)
            curr_particle_index += 1

            p_i = particles[i]

            p_i.w = x[j] * nnls_merging.w_total# * column_norms[j]
            p_i.v = SVector{3, Float64}(vel_pos_matrix[1,j],
                                        vel_pos_matrix[2,j],
                                        vel_pos_matrix[3,j])

            p_i.x = SVector{0,Float64}()
        end
    end

    old_count = indexer.n_local
    n_particles_to_delete = old_count - curr_particle_index

    # if we delete from particles in last cell AND we delete less particles than were in group 2
    # then continuity is not broken
    # !(A && B) == !A || !B
    @inbounds if !(cell == pia.n_cells) || (n_particles_to_delete > indexer.n_group2)
        pia.contiguous[species] = false
    end

    for _ in 1:n_particles_to_delete
        delete_particle_end!(particles, pia, cell, species)
    end

    return 1
end

"""
    merge_nnls_based!(rng, nnls_merging::NNLSMerge{D}, particles::ParticleVector{D}, pia, cell, species;
                      vref=1.0, scaling=:variance,
                      max_err=1e-11, iteration_mult=2, w_threshold=0.0) where D

Perform NNLS-based merging.
The NNLS system is scaled to improve numerical stability, the scaling algorithm is set by the `scaling` parameter.
Even if scaling is done using the computed variances, `vref` might be used in case those variances are small.

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
* `max_err`: maximum allowed value of the residual of the NNLS system
* `iteration_mult`: the number by which the number of columns of the NNLS system matrix is multiplied, this gives the maximum
    number of iterations of the NNLS algorithm
* `w_threshold`: any particles with a relative weight smaller than this value will be discarded (and the weight of the remaining particles re-scaled)

# Returns
If the residual exceeds `max_err` or
the number of non-zero (or smaller than `w_threshold`) elements in the solution vector is equal to the original number of particles,
`-1` is returned to signify a failure of the merging algorithm.

# References
* G. Oblapenko, M. Torrilhon, Moment-preserving particle merging via non-negative least squares.
    [arXiv preprint, 2026](https://doi.org/10.48550/arXiv.2604.00668).
"""
function merge_nnls_based!(rng, nnls_merging::NNLSMerge{D}, particles::ParticleVector{D}, pia, cell, species;
                           vref=1.0, scaling=:variance,
                           max_err=1e-11, iteration_mult=2, w_threshold=0.0) where D
    # create LHS matrix
    @inbounds lhs_ncols = pia.indexer[cell, species].n_local

    if (lhs_ncols > nnls_merging.lhs_matrix_ncols_end) || (lhs_ncols < nnls_merging.lhs_matrix_ncols_start) || length(nnls_merging.vel_pos_matrices) == 0
        indexer = nnls_merging.lhs_matrix_ncols_end - nnls_merging.lhs_matrix_ncols_start + 2

        # we either have pre-allocation [1,2,...[Workspace]]
        # [1,2,Workspace]
        # [Workspace] - we don't pre-allocate if range is exactly 1 particle
        if indexer == 2
            indexer = 1
        end
         
        @inbounds nnls_ws = nnls_merging.work[indexer]
        if (size(nnls_ws.QA, 1) != nnls_merging.n_total_conserved) || (size(nnls_ws.QA, 2) != lhs_ncols)
            nnls_ws.QA = zeros(nnls_merging.n_total_conserved, lhs_ncols)
            resize!(nnls_ws.x, lhs_ncols)
            resize!(nnls_ws.w, lhs_ncols)
            resize!(nnls_ws.idx, lhs_ncols)
        end
        if size(nnls_merging.vel_pos_matrix_scratch, 2) < lhs_ncols
            nnls_merging.vel_pos_matrix_scratch = zeros(3+D, lhs_ncols)
        end
        if length(nnls_merging.column_norms_scratch) < lhs_ncols
            resize!(nnls_merging.column_norms_scratch, lhs_ncols)
        end
        vel_pos_matrix = nnls_merging.vel_pos_matrix_scratch
        column_norms = nnls_merging.column_norms_scratch
    else
        indexer = lhs_ncols - nnls_merging.lhs_matrix_ncols_start + 1
        @inbounds vel_pos_matrix = nnls_merging.vel_pos_matrices[indexer]
        @inbounds column_norms = nnls_merging.column_norms[indexer]
    end
    nnls_merging.vref = vref
    nnls_merging.inv_vref = 1.0 / vref

    # create LHS matrix and fill RHS vector using existing particles
    col_index = compute_lhs_and_rhs!(nnls_merging, nnls_merging.work[indexer].QA, 
                                     vel_pos_matrix, particles, pia, cell, species)


    @inbounds scale_lhs_rhs!(nnls_merging, nnls_merging.work[indexer].QA, scaling, lhs_ncols)

    @inbounds scale_columns!(nnls_merging.work[indexer].QA, column_norms)
    
    nnls_merging.work[indexer].Qb .= nnls_merging.rhs_vector
    @inbounds solve!(nnls_merging.work[indexer], iteration_mult * size(nnls_merging.work[indexer].QA, 2), 1e-14)
    @inbounds return compute_post_merge_particles_nnls!(nnls_merging, nnls_merging.work[indexer].x, particles, pia, cell, species,
                                              lhs_ncols, vel_pos_matrix,
                                              max_err, w_threshold, indexer, column_norms)
end

"""
    merge_nnls_based_rate_preserving!(rng, nnls_merging::NNLSMerge{D},
                                      interaction, electron_neutral_interactions, computed_cs,
                                      particles::ParticleVector{D}, pia, cell, species, neutral_species_index,
                                      ref_cs_elastic, ref_cs_ion; scaling=:variance,
                                      vref=1.0,  max_err=1e-11,
                                      iteration_mult=2,
                                      extend::CSExtend=CSExtendConstant) where D

Perform NNLS-based merging of electrons that conserves approximate elastic scattering and electron-impact ionization rates.
The NNLS system is scaled to improve numerical stability, the scaling algorithm is set by the `scaling` parameter.
Even if scaling is done using the computed variances, `vref` might be used in case those variances are small.

The reference velocity is also used in conjunction with the reference cross-sections to scale the parts of
the NNLS matrix and RHS corresponding to conservation of electron-neutral collision rates. The reference rate is computed
as ``\\sigma_{r,ref} v_{ref}``, where ``\\sigma_{r,ref}`` is the reference process cross-section. In case `scaling==:variance`,
the reference velocity for the computation of reference rates is computed as ``v_{ref} = \\sqrt{E_x^2 + E_y^2 + E_z^2}``,
where ``E_x``, ``E_y``, ``E_z`` are the variances of the velocity in the corresponding directions.

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
* `ref_cs_elastic`: the reference elastic scattering cross-section used to scale the rates
* `ref_cs_ion`: the reference electron-impact ionization cross-section used to scale the rates

# Keyword arguments
* `vref`: the reference velocity used to scale the velocities in the case of `scaling=:vref`
* `scaling`: how to scale entries in the LHS and RHS of the NNLS system - either based on
    the reference velocity `vref` (`scaling=:vref`)
    or on the computed variances in each direction (`scaling=:variance`)
* `max_err`: maximum allowed value of the residual of the NNLS system
* `iteration_mult`: the number by which the number of columns of the NNLS system matrix is multiplied, this gives the maximum
    number of iterations of the NNLS algorithm
* `w_threshold`: any particles with a relative weight smaller than this value will be discarded (and the weight of the remaining particles re-scaled)
* `extend`: enum of `CSExtend` type that sets how out-of-range energy values are treated when computing cross-sections

# Returns
If the residual exceeds `max_err` or
the number of non-zero (or smaller than `w_threshold`) elements in the solution vector is equal to the original number of particles,
`-1` is returned to signify a failure of the merging algorithm.

# References
* G. Oblapenko, M. Torrilhon, Moment-preserving particle merging via non-negative least squares.
    [arXiv preprint, 2026](https://doi.org/10.48550/arXiv.2604.00668).
"""
function merge_nnls_based_rate_preserving!(rng, nnls_merging::NNLSMerge{D},
                                           interaction, electron_neutral_interactions, computed_cs,
                                           particles::ParticleVector{D}, pia, cell, species, neutral_species_index,
                                           ref_cs_elastic, ref_cs_ion; vref=1.0, scaling=:variance,
                                           max_err=1e-11,
                                           iteration_mult=2, w_threshold=0.0,
                                           extend::CSExtend=CSExtendConstant) where D

    # create LHS matrix
    @inbounds lhs_ncols = pia.indexer[cell, species].n_local

    if (lhs_ncols > nnls_merging.lhs_matrix_ncols_end) || (lhs_ncols < nnls_merging.lhs_matrix_ncols_start) || length(nnls_merging.vel_pos_matrices) == 0
        indexer = nnls_merging.lhs_matrix_ncols_end - nnls_merging.lhs_matrix_ncols_start + 2

        # we either have pre-allocation [1,2,...[Workspace]]
        # [1,2,Workspace]
        # [Workspace] - we don't pre-allocate if range is exactly 1 particle
        if indexer == 2
            indexer = 1
        end

        @inbounds nnls_ws = nnls_merging.work[indexer]
        if (size(nnls_ws.QA, 1) != nnls_merging.n_total_conserved) || (size(nnls_ws.QA, 2) != lhs_ncols)
            nnls_ws.QA = zeros(nnls_merging.n_total_conserved, lhs_ncols)
            resize!(nnls_ws.x, lhs_ncols)
            resize!(nnls_ws.w, lhs_ncols)
            resize!(nnls_ws.idx, lhs_ncols)
        end
        if size(nnls_merging.vel_pos_matrix_scratch, 2) < lhs_ncols
            nnls_merging.vel_pos_matrix_scratch = zeros(3+D, lhs_ncols)
        end
        if length(nnls_merging.column_norms_scratch) < lhs_ncols
            resize!(nnls_merging.column_norms_scratch, lhs_ncols)
        end
        vel_pos_matrix = nnls_merging.vel_pos_matrix_scratch
        column_norms = nnls_merging.column_norms_scratch
    else
        indexer = lhs_ncols - nnls_merging.lhs_matrix_ncols_start + 1
        @inbounds vel_pos_matrix = nnls_merging.vel_pos_matrices[indexer]
        @inbounds column_norms = nnls_merging.column_norms[indexer]
    end
    nnls_merging.vref = vref
    nnls_merging.inv_vref = 1.0 / vref

    # create LHS matrix and fill RHS vector using existing particles
    @inbounds col_index = compute_lhs_and_rhs_rate_preserving!(nnls_merging, nnls_merging.work[indexer].QA,
                                                     vel_pos_matrix,
                                                     interaction[species,neutral_species_index],
                                                     electron_neutral_interactions, computed_cs, 
                                                     particles, pia, cell, species, neutral_species_index, extend)
    
    if scaling==:variance
        vr_tmp = sqrt(nnls_merging.Ev[1]^2 + nnls_merging.Ev[2]^2 + nnls_merging.Ev[3]^2)
        ref_k_elastic = ref_cs_elastic * vr_tmp
        ref_k_ion = ref_cs_ion * vr_tmp
    else
        ref_k_elastic = ref_cs_elastic * vref
        ref_k_ion = ref_cs_ion * vref
    end

    @inbounds scale_lhs_rhs_rate_preserving!(nnls_merging, nnls_merging.work[indexer].QA, ref_k_elastic, ref_k_ion, scaling, lhs_ncols)

    @inbounds scale_columns!(nnls_merging.work[indexer].QA, column_norms)

    @inbounds nnls_merging.work[indexer].Qb .= nnls_merging.rhs_vector
    @inbounds solve!(nnls_merging.work[indexer], iteration_mult * size(nnls_merging.work[indexer].QA, 2), 1e-14)
    @inbounds return compute_post_merge_particles_nnls!(nnls_merging, nnls_merging.work[indexer].x, particles, pia, cell, species,
                                              lhs_ncols, vel_pos_matrix,
                                              max_err, w_threshold, indexer, column_norms)
end

"""
    merge_nnls_based_rate_preserving!(rng, nnls_merging::NNLSMerge{D},
                                      interaction, electron_neutral_interactions, computed_cs,
                                      particles::ParticleVector{D}, particles_neutral::ParticleVector{D}, pia, cell, species, neutral_species_index,
                                      ref_cs_elastic, ref_cs_ion; vref=1.0, scaling=:variance,
                                      max_err=1e-11,
                                      iteration_mult=2, w_threshold=0.0,
                                      extend::CSExtend=CSExtendConstant) where D

Perform NNLS-based merging of electrons that conserves **exact** elastic scattering and electron-impact ionization rates
for one specific neutral species.
The NNLS system is scaled to improve numerical stability, the scaling algorithm is set by the `scaling` parameter.
Even if scaling is done using the computed variances, `vref` might be used in case those variances are small.

The reference velocity is also used in conjunction with the reference cross-sections to scale the parts of
the NNLS matrix and RHS corresponding to conservation of electron-neutral collision rates. The reference rate is computed
as ``\\sigma_{r,ref} v_{ref}``, where ``\\sigma_{r,ref}`` is the reference process cross-section. In case `scaling==:variance`,
the reference velocity for the computation of reference rates is computed as ``v_{ref} = \\sqrt{E_x^2 + E_y^2 + E_z^2}``,
where ``E_x``, ``E_y``, ``E_z`` are the variances of the velocity in the corresponding directions.

# Positional arguments
* `rng`: the random number generator instance
* `nnls_merging`: the `NNLSMerge` instance
* `interaction`: the `Interaction` instance describing the electron-neutral interaction being considered
* `electron_neutral_interactions`:  the `ElectronNeutralInteractions` instance storing the tabulated cross-section
    data used to compute the rates
* `computed_cs`: the vector of `ComputedCrossSection` instances in which the computed values will be stored
* `particles`: the `ParticleVector` instance containing the particles to be merged
* `particles_neutral`: the `ParticleVector` instance containing the neutral collision partner particles (they are not affected by the merge)
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the grid cell in which particles are being merged
* `species`: the index of the species being merged
* `neutral_species_index`: the index of the neutral species which is the collision partner in the electron-neutral
    collisions for which approximate rates are being preserved.
* `ref_cs_elastic`: the reference elastic scattering cross-section used to scale the rates
* `ref_cs_ion`: the reference electron-impact ionization cross-section used to scale the rates

# Keyword arguments
* `vref`: the reference velocity used to scale the velocities in the case of `scaling=:vref`
* `scaling`: how to scale entries in the LHS and RHS of the NNLS system - either based on
    the reference velocity `vref` (`scaling=:vref`)
    or on the computed variances in each direction (`scaling=:variance`)
* `max_err`: maximum allowed value of the residual of the NNLS system
* `iteration_mult`: the number by which the number of columns of the NNLS system matrix is multiplied, this gives the maximum
    number of iterations of the NNLS algorithm
* `w_threshold`: any particles with a relative weight smaller than this value will be discarded (and the weight of the remaining particles re-scaled)
* `extend`: enum of `CSExtend` type that sets how out-of-range energy values are treated when computing cross-sections

# Returns
If the residual exceeds `max_err` or
the number of non-zero (or smaller than `w_threshold`) elements in the solution vector is equal to the original number of particles,
`-1` is returned to signify a failure of the merging algorithm.

# References
* G. Oblapenko, M. Torrilhon, Moment-preserving particle merging via non-negative least squares.
    [arXiv preprint, 2026](https://doi.org/10.48550/arXiv.2604.00668).
"""
function merge_nnls_based_rate_preserving!(rng, nnls_merging::NNLSMerge{D},
                                           interaction, electron_neutral_interactions, computed_cs,
                                           particles::ParticleVector{D}, particles_neutral::ParticleVector{D}, pia, cell, species, neutral_species_index,
                                           ref_cs_elastic, ref_cs_ion; vref=1.0, scaling=:variance,
                                           max_err=1e-11,
                                           iteration_mult=2, w_threshold=0.0,
                                           extend::CSExtend=CSExtendConstant) where D

    # create LHS matrix
    @inbounds lhs_ncols = pia.indexer[cell, species].n_local

    if (lhs_ncols > nnls_merging.lhs_matrix_ncols_end) || (lhs_ncols < nnls_merging.lhs_matrix_ncols_start) || length(nnls_merging.vel_pos_matrices) == 0
        indexer = nnls_merging.lhs_matrix_ncols_end - nnls_merging.lhs_matrix_ncols_start + 2

        # we either have pre-allocation [1,2,...[Workspace]]
        # [1,2,Workspace]
        # [Workspace] - we don't pre-allocate if range is exactly 1 particle
        if indexer == 2
            indexer = 1
        end

        @inbounds nnls_ws = nnls_merging.work[indexer]
        if (size(nnls_ws.QA, 1) != nnls_merging.n_total_conserved) || (size(nnls_ws.QA, 2) != lhs_ncols)
            nnls_ws.QA = zeros(nnls_merging.n_total_conserved, lhs_ncols)
            resize!(nnls_ws.x, lhs_ncols)
            resize!(nnls_ws.w, lhs_ncols)
            resize!(nnls_ws.idx, lhs_ncols)
        end
        if size(nnls_merging.vel_pos_matrix_scratch, 2) < lhs_ncols
            nnls_merging.vel_pos_matrix_scratch = zeros(3+D, lhs_ncols)
        end
        if length(nnls_merging.column_norms_scratch) < lhs_ncols
            resize!(nnls_merging.column_norms_scratch, lhs_ncols)
        end
        vel_pos_matrix = nnls_merging.vel_pos_matrix_scratch
        column_norms = nnls_merging.column_norms_scratch
    else
        indexer = lhs_ncols - nnls_merging.lhs_matrix_ncols_start + 1
        @inbounds vel_pos_matrix = nnls_merging.vel_pos_matrices[indexer]
        @inbounds column_norms = nnls_merging.column_norms[indexer]
    end
    nnls_merging.vref = vref
    nnls_merging.inv_vref = 1.0 / vref

    # create LHS matrix and fill RHS vector using existing particles
    @inbounds col_index = compute_lhs_and_rhs_rate_preserving!(nnls_merging, nnls_merging.work[indexer].QA,
                                                     vel_pos_matrix,
                                                     interaction[species,neutral_species_index],
                                                     electron_neutral_interactions, computed_cs, 
                                                     particles, particles_neutral, pia, cell, species, neutral_species_index, extend)
    
    if scaling==:variance
        vr_tmp = sqrt(nnls_merging.Ev[1]^2 + nnls_merging.Ev[2]^2 + nnls_merging.Ev[3]^2)
        ref_k_elastic = ref_cs_elastic * vr_tmp
        ref_k_ion = ref_cs_ion * vr_tmp
    else
        ref_k_elastic = ref_cs_elastic * vref
        ref_k_ion = ref_cs_ion * vref
    end

    @inbounds scale_lhs_rhs_rate_preserving!(nnls_merging, nnls_merging.work[indexer].QA, ref_k_elastic, ref_k_ion, scaling, lhs_ncols)

    @inbounds scale_columns!(nnls_merging.work[indexer].QA, column_norms)

    @inbounds @inbounds nnls_merging.work[indexer].Qb .= nnls_merging.rhs_vector
    solve!(nnls_merging.work[indexer], iteration_mult * size(nnls_merging.work[indexer].QA, 2), 1e-14)
    @inbounds return compute_post_merge_particles_nnls!(nnls_merging, nnls_merging.work[indexer].x, particles, pia, cell, species,
                                              lhs_ncols, vel_pos_matrix,
                                              max_err, w_threshold, indexer, column_norms)
end
end