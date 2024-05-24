function multi_index_moments(n)
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