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
end