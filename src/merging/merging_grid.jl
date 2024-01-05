using StaticArrays

mutable struct GridCell
    w_total::Float64
    v_mean::MVector{3,Float64}
    v_std_sq::MVector{3,Float64}
    x_mean::MVector{3,Float64}
    x_std_sq::MVector{3,Float64}
end

