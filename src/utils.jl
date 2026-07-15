function weighted_percentile_interpolated(values, weights, quantiles=0.5)
    i = sortperm(values)
    c = cumsum(weights[i])

    q = searchsorted(c, quantiles * c[end])
    
    if q.start == q.stop
        return values[i[q.start]]
    else
        if q.start == 0
            return values[i[q.stop]]
        elseif q.stop == 0
            return values[i[q.start]]
        else
            return 0.5 * (values[i[q.start]] + values[i[q.stop]])
        end
    end

    #TODO: figure out why this is producing weird results
end

function weighted_median_with_interpolation(values, weights)
    # assumes length(values) == length(weights) >= 2
    n = length(weights)
    total = sum(weights)
    # sort indices by values
    idx = sortperm(values)
    c = 0.0
    half = total * 0.5
    @inbounds for i in 1:n
        ii = idx[i]
        c += weights[ii]
        if abs(c - half) < 1e-14
            jj = idx[i+1]
            scale = 1.0 / (weights[ii] + weights[jj])
            return scale * (values[ii] * weights[ii] + values[jj] * weights[jj])
        elseif c > half
            # first element is already > 50% of the weight
            if i == 1
                return values[ii]
            else
                jj = idx[i-1]
                scale = 1.0 / (weights[ii] + weights[jj])
                return scale * (values[ii] * weights[ii] + values[jj] * weights[jj])
            end
        end
    end
    # should not get here
    return values[idx[end]]
end

"""
    variance_scaling(var_before::SVector{N,Float64}, var_post::SVector{N,Float64})

Scale the the components of a variance vector `var_before` by the square root of the ratio of the components of it and
a vector `var_post`, checking for zero/negative values in `var_post` and non-finite values in the scaling factor.

A component of `var_post` is treated as zero if it is not larger than `variance_scaling_rel_tol` times the
corresponding component of `var_before`, and a scaling factor of 1.0 is returned for it. A plain `var_post > 0.0`
check is not enough: when all post-merge particles collapse onto the mean (e.g. an N:1 merge of a single octree bin),
`var_post` is zero only in exact arithmetic, and is in practice a round-off residue of order `eps()^2 * var_before`.
Scaling by `sqrt(var_before / var_post)` would then amplify that residue by ~`1/eps()` and destroy the conservation
of the mean, instead of restoring a variance that cannot be restored by scaling in the first place.

# Positional arguments
* `var_before`: variance vector components of which to scale
* `var_post`: variance vector components by which to scale
"""
@inline function variance_scaling(var_before::SVector{N,Float64}, var_post::SVector{N,Float64}) where N
    map(var_before, var_post) do vb, vp
        sf = vp > variance_scaling_rel_tol * vb ? sqrt(vb / vp) : 1.0
        isfinite(sf) ? sf : 1.0
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
    scale_columns!(matrix, column_norms)

Scale the columns of the matrix to have unit L2 norm and store the inverse of the original norm in a vector.

# Positional arguments
* `matrix`: the matrix of the LHS
* `column_norms`: the vector in which to store the computed inverses of the original column-wise norms
"""
function scale_columns!(matrix, column_norms)
    m = size(matrix, 1)
    n = size(matrix, 2)
    @inbounds for i in 1:n
        nn = 0.0
        for j in 1:m
            nn = nn + matrix[j,i] * matrix[j,i]
        end
        nn = sqrt(nn)
        nn = nn > 1e-15 ? nn : 1.0
        column_norms[i] = 1.0 / nn
        for j in 1:m
            matrix[j,i] *= column_norms[i]
        end
    end
end