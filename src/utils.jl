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

"""
    scale_columns!(matrix, column_norms)

Scale the columns of the matrix to have unit L2 norm and store the inverse of the original norm in a vector.

# Positional arguments
* `matrix`: the matrix of the LHS
* `ncols`: number of columns in the matrix
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