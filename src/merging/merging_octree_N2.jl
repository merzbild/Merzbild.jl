using StaticArrays

mutable struct Octree
    max_Nbins::Int64
    Nbins::Int64
    bins::GridCell
end