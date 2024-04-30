using StaticArrays

@enum OctreeMergeSplit MidSplit=1 MeanSplit=2 MedianSplit=3

# struct for N:2 merge
mutable struct OctreeN2
    max_Nbins::Int64
    Nbins::Int64
    bins::GridCell

    # particles in bins[i] have indices in particle_index_buffer[bin_start[i]:bin_end[i]]
    bin_start::Vector{Int64}
    bin_end::Vector{Int64}

    # this stores sorted particle indices
    particle_index_buffer::Vector{Int64}

    split::OctreeMergeSplit
end