@muladd begin

"""
    OctreeBinSplit OctreeBinMidSplit=1 OctreeBinMeanSplit=2 OctreeBinMedianSplit=3
    
Enum defining how the velocity along which the bin is split is chosen.
# Possible values:
* `OctreeBinMidSplit`: the bin is split along the middle velocity
* `OctreeBinMeanSplit`: the bin is split along the mean velocity of the particles in the bin
* `OctreeBinMedianSplit`: the bin is split along the median velocity of the particles in the bin
"""
@enum OctreeBinSplit OctreeBinMidSplit=1 OctreeBinMeanSplit=2 OctreeBinMedianSplit=3

"""
    OctreeInitBin OctreeInitBinMinMaxVel=1 OctreeInitBinMinMaxVelSym=2 OctreeInitBinC=3

Enum defining how the bounds of the initial bin are computed.
# Possible values:
* `OctreeInitBinMinMaxVel`: the minimum and maximum velocities of the particles being merged are used to compute the bounds
* `OctreeInitBinMinMaxVelSym`: the minimum and maximum velocities of the particles being merged are used to compute the bounds,
    but the bounds are then symmetrized in each velocity direction: `[-max(abs(min_v), abs(max_v)), max(abs(min_v), abs(max_v))]`
* `OctreeInitBinC`: the initial bounds are set to `[-c, c]`` in each direction, where `c`` speed of light 
"""
@enum OctreeInitBin OctreeInitBinMinMaxVel=1 OctreeInitBinMinMaxVelSym=2 OctreeInitBinC=3

"""
    OctreeBinBounds OctreeBinBoundsInherit=1 OctreeBinBoundsRecompute=2

Enum defining how the bounds of a split sub-octant bin are computed.
# Possible values:
* `OctreeBinBoundsInherit`: the splitting velocity and the appropriate bounds of the parent bin are inherited
* `OctreeBinBoundsRecompute`: the bounds are recomputed based on the particles in the bin
"""
@enum OctreeBinBounds OctreeBinBoundsInherit=1 OctreeBinBoundsRecompute=2

"""
    OctreeCell

Struct holding computed bin properties needed for refinement of an octree bin.

# Fields
* `np`: number of particles in cell
* `w`: total computational weight of particles in cell
* `v_min`: vector of the per-component lower bounds of the velocities in the cell
* `v_max`: vector of the per-component upper bounds of the velocities in the cell
* `depth`: level of refinement the cell is at
* `can_be_refined`: whether the cell can be refined further
"""
mutable struct OctreeCell
    # this holds only the data needed for refinement
    np::Int64
    w::Float64
    
    v_min::SVector{3,Float64}
    v_max::SVector{3,Float64}

    depth::Int64
    can_be_refined::Bool
end

"""
    OctreeFullCell{D,M}

Struct holding computed bin properties required to merge the particles in a bin down
to M particles.

# Fields
* `v_mean`: mean velocity of particles in cell
* `v_std_sq`: variance of velocity of particles in cell
* `x_mean`: mean position of particles in cell
* `x_std_sq`: variance of position of particles in cell
* `particle_indices`: indices of particles in the cell
* `weights`: the post-merge weights to assign to the particles in the cell
* `velocities`: the post-merge velocities to assign to the particles in the cell
* `positions`: the post-merge positions to assign to the particles in the cell
"""
mutable struct OctreeFullCell{D,M}
    v_mean::SVector{3,Float64}
    v_std_sq::SVector{3,Float64}
    x_mean::SVector{D,Float64}
    x_std_sq::SVector{D,Float64}
    particle_indices::MVector{M, Int64}
    weights::MVector{M, Float64}
    velocities::MVector{M, SVector{3,Float64}} # two dimensional SVector, first index means which particle in the bin, then is it's velocity
    positions::MVector{M, SVector{D,Float64}}
end

# struct for N:2 merge
"""
    OctreeMerge{D,M}

Struct for N:M Octree merging for particles with D-dimensional position vectors.

# Fields
* `max_Nbins`: maximum possible number of bins
* `Nbins`: number of bins currently used
* `bins`: Vector of `OctreeCell` instances used to compute the properties required for bin refinement
* `full_bins`: Vector of `OctreeFullCell{D,M}` instances used to compute the post-merge particles in each bin
* `n_particles`: total number of particles being merged
* `bin_start`: denotes start of indices of particles in bin `i` in the `particle_indexes_sorted` array
* `bin_end`: denotes end of indices of particles in bin `i` in the `particle_indexes_sorted` array
* `particle_indexes_sorted`: Vector of particle indices of the particles being merged
* `particle_octants`: Vector of particle octants for each particle used during radix sort
* `particles_sort_output`: Vector of integer indices used to store particle indices during radix sort
* `particle_in_bin_counter`: `MVector` of size 8, stores the number of particles in each bin
* `nonempty_counter`: `MVector` of size 8, stores the number of particles in each non-empty bin
    (the octants to which these bins correspond to are in `nonempty_bins`)
* `nonempty_bins`: `MVector` of size 8, a sequential list of non-empty octants
* `ndens_counter`: `MVector` of size 8, used in bin splitting, stores number density in each (non-empty) bin
* `bin_bounds_compute`: enum of `OctreeBinBounds` type defining whether bin bounds are fully defined
    by the parent bin and splitting velocity (`vel_middle`), or whether they are recomputed for each new sub-octant bin
* `split`: enum of `OctreeBinSplit` type defining how bins are split
* `vel_middle`: used to store the velocity along which a bin is split into octants
* `v_min_parent`: used in bin splitting to store the vector of the per-component lower bounds of the velocities in the cell
* `v_max_parent`: used in bin splitting to store the vector of the per-component upper bounds of the velocities in the cell
* `direction_vec3`: used to store randomly sampled direction signs, of length 3
* `direction_vecD`: used to store randomly sampled direction signs, of length D
* `init_bin_bounds`: enum of `OctreeInitBin` type defining how the bounds of the top-level bin are set
* `max_depth`: maximum allowed depth of a bin
* `total_post_merge_np`: used to keep track of number of post-merge particles
* `v_mean`: used to store mean velocity for conservative N:1 merging
* `x_mean`: used to store mean position for conservative N:1 merging
* `v_var_before`: used to store pre-merge variance of velocity for conservative N:1 merging
* `x_var_before`: used to store pre-merge variance of position for conservative N:1 merging
* `v_mean_post`: used to store post-merge mean velocity for conservative N:1 merging
* `x_mean_post`: used to store post-merge mean position for conservative N:1 merging
* `v_var_post`: used to store post-merge variance of velocity for conservative N:1 merging
* `x_var_post`: used to store post-merge variance of position for conservative N:1 merging
"""
mutable struct OctreeMerge{D,M}
    max_Nbins::Int64
    Nbins::Int64  # actual bins computed
    bins::Vector{OctreeCell}
    full_bins::Vector{OctreeFullCell{D,M}}
    n_particles::Int64  # particles being sorted

    # particles in bins[i] have indices in particle_index_buffer[bin_start[i]:bin_end[i]]
    bin_start::Vector{Int64}
    bin_end::Vector{Int64}

    # this stores particle indices in the cell we're merging in, default size is 8192
    # in case particles don't fit it's increase to length(particles) + DELTA_PARTICLES
    particle_indexes_sorted::Vector{Int64}

    # used to store particle octants during radix sort, default size is 8192
    # in case particles don't fit it's increase to length(particles) + DELTA_PARTICLES
    particle_octants::Vector{Int64}

    # used to store particle indices during radix sort, default size is 8192
    # in case particles don't fit it's increase to length(particles) + DELTA_PARTICLES
    particles_sort_output::Vector{Int64}


    # count how many particles in each bin, used for radix sort
    particle_in_bin_counter::MVector{8, Int64}

    # this stores # of particle in each non-empty bin sequentially (without knowledge of which bin this belongs to)
    nonempty_counter::MVector{8, Int64}

    # a sequential list of non-empty octants
    nonempty_bins::MVector{8, Int64}

    # this stores number density in each non-empty bin 
    ndens_counter::MVector{8, Float64}

    bin_bounds_compute::OctreeBinBounds
    split::OctreeBinSplit
    vel_middle::SVector{3,Float64}  # defines how we split octant
    v_min_parent::SVector{3,Float64}  # used in splitting
    v_max_parent::SVector{3,Float64}

    direction_vec3::SVector{3,Float64}
    direction_vecD::SVector{D,Float64}

    init_bin_bounds::OctreeInitBin

    max_depth::Int64
    total_post_merge_np::Int64 # used to keep track of number of post-merge particles

    # N:1 case --> Pre-allocated fields for compute_new_particles! to avoid allocations
    v_mean::SVector{3, Float64}
    x_mean::SVector{D, Float64}
    v_var_before::SVector{3, Float64}
    x_var_before::SVector{D, Float64}
    v_mean_post::SVector{3, Float64}
    x_mean_post::SVector{D, Float64}
    v_var_post::SVector{3, Float64}
    x_var_post::SVector{D, Float64}
end

"""
    fill_bins(Nbins)
    
Fill the octree bins structs with zero data, used as a utility function for initialization.

Positional arguments:
* `Nbins`: number of `OctreeCell` bins to create

Returns:
An array of `Nbins` `OctreeCell` instances filled with zeros.
"""
function fill_bins(Nbins)
    return [OctreeCell(0, 0.0, SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0), 0, true) for i in 1:Nbins]
end

"""
    fill_full_bins(::Val{D}, ::Val{M}, Nbins)

Fill the octree bins full structs with zero data, used as a utility function for initialization.

Positional arguments:
* `D`: dimension of position vectors of particles to be merged
* `M`: number of post-merge particles in each bin
* `Nbins`: number of `OctreeFullCell` bins to create

Returns:
An array of `Nbins` `OctreeFullCell` instances filled with zeros.
"""
function fill_full_bins(::Val{D}, ::Val{M}, Nbins) where {D, M}
    return [OctreeFullCell{D,M}(SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0),
                           zero(SVector{D,Float64}), zero(SVector{D,Float64}),
                           MVector{M, Int64}(zeros(Int64, M)),
                           MVector{M, Float64}(zeros(Float64, M)), # weights
                           MVector{M, SVector{3,Float64}}([SVector{3,Float64}(0.0, 0.0, 0.0) for _ in 1:M]),
                           MVector{M, SVector{D,Float64}}([zero(SVector{D,Float64}) for _ in 1:M]))
                           for _ in 1:Nbins]
end

"""
    OctreeMerge{D,M}(split::OctreeBinSplit; init_bin_bounds=OctreeInitBinMinMaxVel, bin_bounds_compute=OctreeBinBoundsInherit,
              max_Nbins=4096, max_depth=10)
    
Create an Octree N:M merging instance for particles with D-dimensional position vectors.

Positional arguments:
* `split`: a enum of `OctreeBinSplit` type which tells how to split a bin into sub-bins

Keyword arguments:
* `init_bin_bounds`: a enum of `OctreeInitBin` type which defines how the bounds of the top-level bin are set
* `bin_bounds_compute`: a enum of `OctreeBinBounds` type which defines whether the bounds of sub-bins are recomputed
    based on the minimum/maximum velocities of the particles in those sub-bins, or the bounds are inherited from the
    bin that was split
* `max_Nbins`: maximum number of bins allowed (this only counts leaf-level bins)
* `max_depth`: maximum depth of a sub-bin starting from the top-level bin containing all particles (which has a depth of 0)

Returns:
`OctreeMerge` instance with everything set to 0.
"""
function OctreeMerge{D,M}(split::OctreeBinSplit; init_bin_bounds=OctreeInitBinMinMaxVel, bin_bounds_compute=OctreeBinBoundsInherit,
                 max_Nbins=4096, max_depth=10) where {D, M}
    return OctreeMerge{D,M}(max_Nbins, 0, fill_bins(max_Nbins),
                            fill_full_bins(Val(D), Val(M), max_Nbins), 0,
                            zeros(max_Nbins), zeros(max_Nbins),  # bin_start, bin_end
                            zeros(8192), zeros(8192), zeros(8192),
                            MVector{8, Int64}(0, 0, 0, 0, 0, 0, 0, 0),
                            MVector{8, Int64}(0, 0, 0, 0, 0, 0, 0, 0),
                            MVector{8, Int64}(0, 0, 0, 0, 0, 0, 0, 0),
                            MVector{8, Float64}(0, 0, 0, 0, 0, 0, 0, 0),
                            bin_bounds_compute, split,
                            SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0),
                            zero(SVector{3,Float64}),
                            zero(SVector{D,Float64}),
                            init_bin_bounds, max_depth, 0,    
                            SVector{3,Float64}(0.0, 0.0, 0.0),  # v_mean
                            zero(SVector{D,Float64}),  # x_mean
                            SVector{3,Float64}(0.0, 0.0, 0.0),  # v_var_before
                            zero(SVector{D,Float64}),  # x_var_before
                            SVector{3,Float64}(0.0, 0.0, 0.0),  # v_mean_post
                            zero(SVector{D,Float64}),  # x_mean_post
                            SVector{3,Float64}(0.0, 0.0, 0.0),  # v_var_post
                            zero(SVector{D,Float64})   # x_var_post
)
end


"""
    OctreeMerge(split::OctreeBinSplit; init_bin_bounds=OctreeInitBinMinMaxVel, bin_bounds_compute=OctreeBinBoundsInherit,
              max_Nbins=4096, max_depth=10)
    
Create an Octree N:2 merging instance for particles with 3-dimensional position vectors.

Positional arguments:
* `split`: a enum of `OctreeBinSplit` type which tells how to split a bin into sub-bins

Keyword arguments:
* `init_bin_bounds`: a enum of `OctreeInitBin` type which defines how the bounds of the top-level bin are set
* `bin_bounds_compute`: a enum of `OctreeBinBounds` type which defines whether the bounds of sub-bins are recomputed
    based on the minimum/maximum velocities of the particles in those sub-bins, or the bounds are inherited from the
    bin that was split
* `max_Nbins`: maximum number of bins allowed (this only counts leaf-level bins)
* `max_depth`: maximum depth of a sub-bin starting from the top-level bin containing all particles (which has a depth of 0)

Returns:
`OctreeMerge` instance with everything set to 0.
"""
OctreeMerge(split::OctreeBinSplit; init_bin_bounds=OctreeInitBinMinMaxVel, bin_bounds_compute=OctreeBinBoundsInherit,
              max_Nbins=4096, max_depth=10) = OctreeMerge{3,2}(split::OctreeBinSplit; init_bin_bounds=init_bin_bounds, bin_bounds_compute=bin_bounds_compute,
                 max_Nbins=max_Nbins, max_depth=max_depth)

"""
    clear_octree!(octree)
    
Reset octree before doing a new merge.

# Positional arguments:
* `octree`: the `OctreeMerge` instance
"""
function clear_octree!(octree)
    octree.Nbins = 0
end

"""
    resize_octree_buffers!(octree, n_particles)

Check and resize octree buffers if needed to accommodate a larger number of particles. The size of the buffers
is set to `n_particles + DELTA_PARTICLES` if the their sizes are smaller than `n_particles`.

# Positional arguments:
* `octree`: the `OctreeMerge` instance
* `n_particles`: the current number of particles being merged
"""
function resize_octree_buffers!(octree, n_particles)
    if (length(octree.particle_indexes_sorted) < n_particles)
        resize!(octree.particle_indexes_sorted, n_particles + DELTA_PARTICLES)
    end
    if (length(octree.particle_octants) < n_particles)
        resize!(octree.particle_octants, n_particles + DELTA_PARTICLES)
    end
    if (length(octree.particles_sort_output) < n_particles)
        resize!(octree.particles_sort_output, n_particles + DELTA_PARTICLES)
    end
end

"""
    compute_octant(particle_v, v_middle)

Compute octant of particle velocity relative to a `v_middle``
The order of the octants is:

    1. - - -
    2. + - -
    3. - + -
    4. + + -
    5. - - +
    6. + - +
    7. - + +
    8. + + +

# Positional arguments:
* `particle_v`: the velocity of the particle
* `v_middle`: the velocity relative to which the octant is computed

# Returns:
The octant number
"""
@inline function compute_octant(particle_v, v_middle)
    @inbounds return 1 + (particle_v[1] > v_middle[1]) + 2 * (particle_v[2] > v_middle[2]) + 4 * (particle_v[3] > v_middle[3])
end

"""
    bin_bounds_inherit!(octree, bin_id, v_min_parent, v_max_parent, v_middle, octant)

Compute new bin bounds of one of the 8 octant sub-bins inheriting bounds of parent bin. 
The order of the octants is

    1. - - -
    2. + - -
    3. - + -
    4. + + -
    5. - - +
    6. + - +
    7. - + +
    8. + + +

# Positional arguments:
* `octree`: the `OctreeMerge` instance
* `bin_id`: the index of the bin for which the velocity bounds are being recomputed
* `v_min_parent`: the velocity vector of the lower bound of the velocities of the parent bin
* `v_max_parent`: the velocity vector of the upper bound of the velocities of the parent bin
* `v_middle`: the velocity across which the split is being performed
* `octant`: the octant of the parent bin to which this bin corresponds
"""
function bin_bounds_inherit!(octree, bin_id, v_min_parent, v_max_parent, v_middle, octant)
    bin = octree.bins[bin_id]

    if octant == 1
        @inbounds bin.v_min = v_min_parent
        @inbounds bin.v_max = v_middle
    elseif octant == 2
        @inbounds bin.v_min = SVector{3,Float64}(v_middle[1], v_min_parent[2], v_min_parent[3])
        @inbounds bin.v_max = SVector{3,Float64}(v_max_parent[1], v_middle[2], v_middle[3])
    elseif octant == 3
        @inbounds bin.v_min = SVector{3,Float64}(v_min_parent[1], v_middle[2], v_min_parent[3])
        @inbounds bin.v_max = SVector{3,Float64}(v_middle[1], v_max_parent[2], v_middle[3])
    elseif octant == 4
        @inbounds bin.v_min = SVector{3,Float64}(v_middle[1], v_middle[2], v_min_parent[3])
        @inbounds bin.v_max = SVector{3,Float64}(v_max_parent[1], v_max_parent[2], v_middle[3])
    elseif octant == 5
        @inbounds bin.v_min = SVector{3,Float64}(v_min_parent[1], v_min_parent[2], v_middle[3])
        @inbounds bin.v_max = SVector{3,Float64}(v_middle[1], v_middle[2], v_max_parent[3])
    elseif octant == 6
        @inbounds bin.v_min = SVector{3,Float64}(v_middle[1], v_min_parent[2], v_middle[3])
        @inbounds bin.v_max = SVector{3,Float64}(v_max_parent[1], v_middle[2], v_max_parent[3])
    elseif octant == 7
        @inbounds bin.v_min = SVector{3,Float64}(v_min_parent[1], v_middle[2], v_middle[3])
        @inbounds bin.v_max = SVector{3,Float64}(v_middle[1], v_max_parent[2], v_max_parent[3])
    else
        @inbounds bin.v_min = v_middle
        @inbounds bin.v_max = v_max_parent
    end
end

"""
    bin_bounds_recompute!(octree, bin_id, bs, be, particles)

Recompute bin bounds based on particle velocities by setting them to the smallest and largest
velocities of the particles in each velocity direction.

# Positional arguments:
* `octree`: the `OctreeMerge` instance
* `bin_id`: the index of the bin for which the velocity bounds are being recomputed
* `bs`: index of the first particle in bin
* `be`: index of the last particle in bin
* `particles`: the `ParticleVector` instance of the particles to be merged
"""
function bin_bounds_recompute!(octree, bin_id, bs, be, particles::ParticleVector{D}) where D
    # compute bin bounds based on the particles in the bin
    minvx = 9_299_792_458.0  # speed of light + 9e9
    minvy = 9_299_792_458.0
    minvz = 9_299_792_458.0

    maxvx = -9_299_792_458.0
    maxvy = -9_299_792_458.0
    maxvz = -9_299_792_458.0

    @inbounds for i in bs:be
        pin = octree.particle_indexes_sorted[i]
        
        pv = particles[pin].v

        if (pv[1] < minvx)
            minvx = pv[1]
        end
        if (pv[1] > maxvx)
            maxvx = pv[1]
        end

        if (pv[2] < minvy)
            minvy = pv[2]
        end
        if (pv[2] > maxvy)
            maxvy = pv[2]
        end

        if (pv[3] < minvz)
            minvz = pv[3]
        end
        if (pv[3] > maxvz)
            maxvz = pv[3]
        end
    end

    @inbounds octree.bins[bin_id].v_min = SVector{3, Float64}(minvx, minvy, minvz)
    @inbounds octree.bins[bin_id].v_max = SVector{3, Float64}(maxvx, maxvy, maxvz)
end

"""
    compute_v_mean!(octree, bs, be, particles)

Compute mean velocity of particles in a bin.

# Positional arguments:
* `octree`: the `OctreeMerge` instance
* `bs`: index of the first particle in bin
* `be`: index of the last particle in bin
* `particles`: the `ParticleVector` instance of the particles to be merged
"""
@inline function compute_v_mean!(octree, bs, be, particles::ParticleVector{D}) where D
    n_tot = 0.0
    v_mean = SVector{3,Float64}(0.0, 0.0, 0.0)
    @inbounds for i in bs:be
        p = particles[octree.particle_indexes_sorted[i]]
        n_tot += p.w
        v_mean = v_mean + p.w * p.v
    end
    octree.vel_middle = v_mean / n_tot
end

"""
    compute_v_median!(octree, bs, be, particles)

Compute median velocity of particles in a bin. NOTE: allocates memory and is probably not fully correct!

# Positional arguments:
* `octree`: the `OctreeMerge` instance
* `bs`: index of the first particle in bin
* `be`: index of the last particle in bin
* `particles`: the `ParticleVector` instance of the particles to be merged
"""
function compute_v_median!(octree, bs, be, particles::ParticleVector{D}) where D
    # https://rcoh.me/posts/linear-time-median-finding/
    # https://en.wikipedia.org/wiki/Weighted_median
    @inbounds w_vec = [particles[octree.particle_indexes_sorted[i]].w for i in bs:be]
    w_vec = w_vec ./ sum(w_vec)

    @inbounds vx_vec = [particles[octree.particle_indexes_sorted[i]].v[1] for i in bs:be]
    @inbounds vy_vec = [particles[octree.particle_indexes_sorted[i]].v[2] for i in bs:be]
    @inbounds vz_vec = [particles[octree.particle_indexes_sorted[i]].v[3] for i in bs:be]

    octree.vel_middle = SVector{3, Float64}(weighted_median_with_interpolation(vx_vec, w_vec),
                                            weighted_median_with_interpolation(vy_vec, w_vec),
                                            weighted_median_with_interpolation(vz_vec, w_vec))                                  
end

function bin_bounds_recompute_and_v_mean!(octree, bin_id, bs, be, particles::ParticleVector{D}) where D
    # do everything in 1 pass over the particles
    # TODO
end

"""
    get_new_bin_id(i, bin_id, Nbins)

Get index of a newly created bin once a bin with index `bin_id` is split into 8 sub-bins. The
`bin_id` index is re-used for the 1-st sub-bin, and the other 7 sub-bins are tacked onto the end of the list of bins
(so they have indices Nbins + 1, Nbins + 2, ..., where Nbins was the total number of octree bins before the split).

# Positional arguments
* `i`: index of sub-octant (ranging from 1 to 8) of the split bin for which to return a new bin index
* `bin_id`: index of bin being split
* `Nbins`: total number of octree bins before the split

# Returns
Index of a newly created bin corresponding to a bin created from sub-octant `i` of bin `bin_id`.
"""
@inline function get_new_bin_id(i, bin_id, Nbins)
    return i == 1 ? bin_id : Nbins + i - 1
end

"""
    split_bin!(octree::OctreeMerge{D,M}, bin_id, particles::ParticleVector{D})

Sort particles into sub-bins of a bin with index `bin_id` (by splitting it into octants),
keeping track of which sub-bins particles end up in.
Also sets the velocity bounds of the new bins.

# Positional arguments:
* `octree`: the `OctreeMerge` instance
* `bin_id`: octree bin index
* `particles`: the `ParticleVector` instance of the particles to be merged
"""
function split_bin!(octree::OctreeMerge{D,M}, bin_id, particles::ParticleVector{D}) where {D,M}
    n_nonempty_bins = 0
    
    # the bin we're splitting
    bin0 = octree.bins[bin_id]

    particle_in_bin_counter = octree.particle_in_bin_counter
    ndens_counter = octree.ndens_counter
    nonempty_bins = octree.nonempty_bins
    nonempty_counter = octree.nonempty_counter

    fill!(particle_in_bin_counter, 0)
    fill!(ndens_counter, 0.0)
    fill!(nonempty_bins, 0)
    fill!(nonempty_counter, 0)
    
    @inbounds current_depth = bin0.depth
    @inbounds bs = octree.bin_start[bin_id]
    @inbounds be = octree.bin_end[bin_id]

    # if we compute bin bounds based on particle velocities, we need to actually iterate through particles
    # and compute
    # if not, we can do it later directly for the sub-bins
    if (octree.bin_bounds_compute == OctreeBinBoundsRecompute)
        bin_bounds_recompute!(octree, bin_id, bs, be, particles)
    end

    if (octree.split == OctreeBinMidSplit)
        @inbounds octree.vel_middle = 0.5 * (bin0.v_min + bin0.v_max)
    elseif (octree.split == OctreeBinMeanSplit)
        compute_v_mean!(octree, bs, be, particles) # octree.vel_middle = octree.bins[bin_id].v_mean
    elseif (octree.split == OctreeBinMedianSplit)
        compute_v_median!(octree, bs, be, particles)
    end

    vel_middle = octree.vel_middle
    particle_octants = octree.particle_octants
    particle_indexes_sorted = octree.particle_indexes_sorted
    @inbounds for i in bs:be
        p = particles[particle_indexes_sorted[i]]
        oct = compute_octant(p.v, vel_middle)
        particle_in_bin_counter[oct] += 1
        particle_octants[i-bs+1] = oct
        ndens_counter[oct] += p.w
    end

    n_eb = 0
    @inbounds if (particle_in_bin_counter[1] > 0)
        n_nonempty_bins += 1
        n_eb += 1
        @inbounds nonempty_counter[n_eb] = particle_in_bin_counter[1]
        @inbounds nonempty_bins[n_eb] = 1
    end

    @inbounds for i in 2:8
        if (particle_in_bin_counter[i] > 0)
            n_nonempty_bins += 1
            n_eb += 1
            nonempty_counter[n_eb] = particle_in_bin_counter[i]
            nonempty_bins[n_eb] = i
        end

        particle_in_bin_counter[i] += particle_in_bin_counter[i-1]
    end

    # first bin - we change nothing for the start
    @inbounds octree.bin_end[bin_id] = octree.bin_start[bin_id] + nonempty_counter[1] - 1
    
    # new bins will point to the sorted particles, but are not contiguous
    # first bin is in the old place and the new ones are tacked on
    @inbounds for i in 2:n_nonempty_bins
        bi = get_new_bin_id(i, bin_id, octree.Nbins)
        bim1 = get_new_bin_id(i-1, bin_id, octree.Nbins)
        octree.bin_start[bi] = octree.bin_end[bim1] + 1
        octree.bin_end[bi] = octree.bin_start[bi] + nonempty_counter[i] - 1
    end

    # we had a bin that would've produced M particles
    # now we replaced ith with n_nonempty_bins bins that each produce 1 to M particles
    octree.total_post_merge_np -= M
    if (octree.bin_bounds_compute == OctreeBinBoundsInherit)
        @inbounds octree.v_min_parent = bin0.v_min
        @inbounds octree.v_max_parent = bin0.v_max

        # iterate over non-empty bins and inherit parent bin bounds + split around middle velocity
        @inbounds for i in 1:n_nonempty_bins
            bi = get_new_bin_id(i, bin_id, octree.Nbins)

            bin = octree.bins[bi]
            neb = nonempty_bins[i]

            bin_bounds_inherit!(octree, bi,
                                octree.v_min_parent, octree.v_max_parent,
                                octree.vel_middle, neb)
            bin.np = nonempty_counter[i]
            bin.w = ndens_counter[neb]
            bin.depth = current_depth + 1

            octree.total_post_merge_np += get_bin_post_merge_np(octree, bi)
            # octree.bins[bin_id + i - 1].post_merge_np = get_bin_post_merge_np(octree, bin_id + i - 1)
            if (bin.np > M) && (bin.depth < octree.max_depth)
                bin.can_be_refined = true
            else
                bin.can_be_refined = false
            end
        end
    else
        # still need to fill out info on number of particles and total weight
        # will recompute bin bounds if we do next round of refinement
        @inbounds for i in 1:n_nonempty_bins
            bi = get_new_bin_id(i, bin_id, octree.Nbins)
            bin = octree.bins[bi]
            bin.np = nonempty_counter[i]
            bin.w = ndens_counter[nonempty_bins[i]]
            bin.depth = current_depth + 1

            # we had a bin that would've produced 2 particles
            # now we replaced ith with n_nonempty_bins bins that each produce 1 or 2 particles
            octree.total_post_merge_np += get_bin_post_merge_np(octree, bi)
            if (bin.np > M) && (bin.depth < octree.max_depth)
                bin.can_be_refined = true
            else
                bin.can_be_refined = false
            end
        end
    end

    particles_sort_output = octree.particles_sort_output
    @inbounds for i in bs:be
        pin = particle_indexes_sorted[i]
        j = particle_octants[i - bs + 1]
        particles_sort_output[particle_in_bin_counter[j]] = pin
        particle_in_bin_counter[j] -= 1
    end

    # write sorted indices
    unsafe_copyto!(particle_indexes_sorted, bs, particles_sort_output, 1, be - bs + 1)

    octree.Nbins += n_nonempty_bins - 1
end

"""
    compute_bin_props!(octree::OctreeMerge{D,M}, bin_id, particles::ParticleVector{D})

Compute properties in a bin required for merging: total computational weight, mean velocity and position,
standard deviation of particle velocities and positions.

# Positional arguments:
* `octree`: the `OctreeMerge` instance
* `bin_id`: octree bin index
* `particles`: the `ParticleVector` instance of the particles to be merged
"""
function compute_bin_props!(octree::OctreeMerge{D,M}, bin_id, particles::ParticleVector{D}) where {D,M}
    @inbounds bs = octree.bin_start[bin_id]
    @inbounds be = octree.bin_end[bin_id]

    @inbounds bin = octree.bins[bin_id]
    @inbounds full_bin = octree.full_bins[bin_id]

    @inbounds if (bin.w) == 0
        # we can discard any particles in the bin
        bin.np = 0
        return 
    end

    # store indices of the first M particles in octree bin so that we
    # have somewhere to write post-merge data
    @inbounds actual_n = get_bin_post_merge_np(octree, bin_id)
    # for general M particles - mutate the MVector in place to avoid allocating a new one
    @inbounds for i in 1:M
        full_bin.particle_indices[i] = i <= actual_n ? octree.particle_indexes_sorted[bs + i - 1] : 0
    end

    # if only M or fewer particles in bin then we don't need to compute any properties
    @inbounds if (octree.bins[bin_id].np <= M)
        return
    end

    v_mean = SVector{3, Float64}(0.0, 0.0, 0.0)
    v_std_sq = SVector{3, Float64}(0.0, 0.0, 0.0)

    x_mean = zero(SVector{D, Float64})
    x_std_sq = zero(SVector{D, Float64})

    inv_w = 1.0 / octree.bins[bin_id].w

    @inbounds for i in bs:be
        p_i = particles[octree.particle_indexes_sorted[i]]
        v_mean = v_mean + p_i.w * p_i.v
        x_mean = x_mean + p_i.w * p_i.x
    end
    v_mean = v_mean * inv_w
    x_mean = x_mean * inv_w

    @inbounds for i in bs:be
        p_i = particles[octree.particle_indexes_sorted[i]]
        v_std_sq = v_std_sq + p_i.w * (p_i.v - v_mean).^2
        x_std_sq = x_std_sq + p_i.w * (p_i.x - x_mean).^2
    end
    v_std_sq = v_std_sq * inv_w
    x_std_sq = x_std_sq * inv_w

    full_bin.v_mean = v_mean
    full_bin.v_std_sq = v_std_sq

    full_bin.x_mean = x_mean
    full_bin.x_std_sq = x_std_sq
end

"""
    get_bin_post_merge_np(octree::OctreeMerge{D,M}, bin_id)

Get number of post-merge particles in a bin: M if the number of particles in the bin is >= M, otherwise
the number of particles in the bin is returned.

# Positional arguments:
* `octree`: the `OctreeMerge` instance
* `bin_id`: octree bin index

# Returns
The number of post-merge particles in a single octree bin (0, 1, or 2).
"""
@inline function get_bin_post_merge_np(octree::OctreeMerge{D,M}, bin_id) where {D,M}
    # how many particles will we get after merging in the bin:
    # M if np >= M
    # np otherwise
    @inbounds return octree.bins[bin_id].np >= M ? M : octree.bins[bin_id].np
end

"""
    compute_new_particles!(rng, octree::OctreeMerge{D,2}, particles::ParticleVector{D}, pia, cell, species)

Compute post-merge particles with particles based on octree bin properties without checking or setting particle locations
(for spatially homogeneous merging).
N:2 merging in each bin.

# Positional arguments:
* `rng`: the random number generator instance
* `octree`: the `OctreeMerge` instance
* `particles`: the `ParticleVector` instance of the particles to be merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the cell index
* `species`: the species index
"""
function compute_new_particles!(rng, octree::OctreeMerge{D,2}, particles::ParticleVector{D}, pia, cell, species) where D
    # given computed Octree, create new particles instead of the old ones
    
    Nbins = octree.Nbins
    @inbounds for bin_id in 1:Nbins
        bin = octree.bins[bin_id]
        loc_np = bin.np
        full_bin = octree.full_bins[bin_id]
        if (loc_np > 2)
            full_bin.weights[1] = 0.5 * bin.w
            full_bin.weights[2] = full_bin.weights[1]

            full_bin.v_std_sq = sqrt.(full_bin.v_std_sq)
            # full_bin.x_std_sq = sqrt.(full_bin.x_std_sq)
            
            octree.direction_vec3 = @SVector rand(rng, direction_signs, 3)
            full_bin.velocities[1] = full_bin.v_mean + octree.direction_vec3 .* full_bin.v_std_sq
            full_bin.velocities[2] = full_bin.v_mean - octree.direction_vec3 .* full_bin.v_std_sq

            # octree.direction_vec = @SVector rand(rng, direction_signs, 3)
            # octree.full_bins[bin_id].x1 = octree.full_bins[bin_id].x_mean + octree.direction_vec .* octree.full_bins[bin_id].x_std_sq
            # octree.full_bins[bin_id].x2 = octree.full_bins[bin_id].x_mean - octree.direction_vec .* octree.full_bins[bin_id].x_std_sq
        elseif (loc_np == 2)
            # get the particle indices we saved and just write data based on them
            i = full_bin.particle_indices[1]
            full_bin.weights[1] = particles[i].w
            full_bin.velocities[1] = particles[i].v
            # full_bin.positions[1] = particles[i].x

            i = full_bin.particle_indices[2]
            full_bin.weights[2] = particles[i].w
            full_bin.velocities[2] = particles[i].v
            # full_bin.positions[2] = particles[i].x
        elseif (loc_np == 1)
            # get the particle indices we saved and just write data based on them
            i = full_bin.particle_indices[1]
            full_bin.weights[1] = particles[i].w
            full_bin.velocities[1] = particles[i].v
            # full_bin.positions[1] = particles[i].x
        end
    end

    curr_particle_index = write_back_to_particles!(octree, particles, pia, cell, species)

    @inbounds indexer = pia.indexer[cell,species]
    @inbounds old_count = indexer.n_local
    n_particles_to_delete = old_count - curr_particle_index

    # if we delete from particles in last cell AND we delete less particles than were in group 2
    # then continuity is not broken
    # !(A && B) == !A || !B
    @inbounds if !(cell == pia.n_cells) || (n_particles_to_delete > pia.indexer[cell,species].n_group2)
        pia.contiguous[species] = false
    end

    for _ in 1:n_particles_to_delete
        delete_particle_end!(particles, pia, cell, species)
    end
end

"""
    compute_new_particles!(rng, octree::OctreeMerge{D,2}, particles::ParticleVector{D}, pia, cell, species, grid::Grid1DUniform)

Compute post-merge particles particles based on octree bin properties; placing out-of-domain particles back into the domain.
N:2 merging in each bin.

# Positional arguments:
* `rng`: the random number generator instance
* `octree`: the `OctreeMerge` instance
* `particles`: the `ParticleVector` instance of the particles to be merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the cell index
* `species`: the species index
* `grid`: the `Grid1DUniform` grid
"""
function compute_new_particles!(rng, octree::OctreeMerge{D,2}, particles::ParticleVector{D}, pia, cell, species, grid::Grid1DUniform) where D
    # given computed Octree, create new particles instead of the old ones
    
    @inbounds for bin_id in 1:octree.Nbins
        bin = octree.bins[bin_id]
        loc_np = bin.np
        full_bin = octree.full_bins[bin_id]
        if (loc_np > 2)
            full_bin.weights[1] = 0.5 * bin.w
            full_bin.weights[2] = full_bin.weights[1]

            full_bin.v_std_sq = sqrt.(full_bin.v_std_sq)
            full_bin.x_std_sq = sqrt.(full_bin.x_std_sq)
            
            octree.direction_vec3 = @SVector rand(rng, direction_signs, 3)
            full_bin.velocities[1] = full_bin.v_mean + octree.direction_vec3 .* full_bin.v_std_sq
            full_bin.velocities[2] = full_bin.v_mean - octree.direction_vec3 .* full_bin.v_std_sq

            octree.direction_vecD = @SVector rand(rng, direction_signs, D)
            full_bin.positions[1] = full_bin.x_mean + octree.direction_vecD .* full_bin.x_std_sq
            full_bin.positions[2] = full_bin.x_mean - octree.direction_vecD .* full_bin.x_std_sq
        elseif (loc_np == 2)
            # get the particle indices we saved and just write data based on them
            i = full_bin.particle_indices[1]
            full_bin.weights[1] = particles[i].w
            full_bin.velocities[1] = particles[i].v
            full_bin.positions[1] = particles[i].x

            i = full_bin.particle_indices[2]
            full_bin.weights[2] = particles[i].w
            full_bin.velocities[2] = particles[i].v
            full_bin.positions[2] = particles[i].x
        elseif (loc_np == 1)
            # get the particle indices we saved and just write data based on them
            i = full_bin.particle_indices[1]
            full_bin.weights[1] = particles[i].w
            full_bin.velocities[1] = particles[i].v
            full_bin.positions[1] = particles[i].x
        end
    end

    curr_particle_index = write_back_to_particles!(octree, particles, pia, cell, species, grid)

    @inbounds indexer = pia.indexer[cell,species]
    @inbounds old_count = indexer.n_local
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
end

"""
    compute_new_particles!(rng, octree::OctreeMerge{D,1}, particles::ParticleVector{D}, pia, cell, species)

Compute post-merge particles with particles based on octree bin properties without checking or setting particle locations
(for spatially homogeneous merging).
N:1 merging in each bin.
It computes properties of all particles. And after the first merge, it scales the particles' velocities based on these properties to achieve conservation.

# Positional arguments:
* `rng`: the random number generator instance
* `octree`: the `OctreeMerge` instance
* `particles`: the `ParticleVector` instance of the particles to be merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the grid cell in which particles are being merged
* `species`: the index of the species being merged
"""
function compute_new_particles!(rng, octree::OctreeMerge{D,1}, particles::ParticleVector{D}, pia, cell, species) where D
    # --- 1. Compute statistics before merging ---
    # bin-level properties cannot be used here: compute_bin_props! skips bins with np <= M,
    # so full_bins[bin_id].v_mean is stale for single-particle bins when M == 1
    w_total = 0.0
    octree.v_mean = SVector{3,Float64}(0.0, 0.0, 0.0)

    @inbounds for bin_id in 1:octree.Nbins
        bs, be = octree.bin_start[bin_id], octree.bin_end[bin_id]
        for i in bs:be
            p_i = particles[octree.particle_indexes_sorted[i]]
            w_total += p_i.w
            octree.v_mean = octree.v_mean + p_i.w * p_i.v
        end
    end

    octree.v_mean = octree.v_mean / w_total

    octree.v_var_before = SVector{3,Float64}(0.0, 0.0, 0.0)

    @inbounds for bin_id in 1:octree.Nbins
        bs, be = octree.bin_start[bin_id], octree.bin_end[bin_id]
        for i in bs:be
            p_i = particles[octree.particle_indexes_sorted[i]]
            octree.v_var_before = octree.v_var_before + p_i.w * (p_i.v - octree.v_mean).^2
        end
    end

    octree.v_var_before = octree.v_var_before / w_total

    # --- 2. Generate post-merge particle data into full_bins ---
    @inbounds for bin_id in 1:octree.Nbins
        loc_np = octree.bins[bin_id].np
        fb = octree.full_bins[bin_id]
        if (loc_np > 1)
            # N:1 merge
            fb.weights[1] = octree.bins[bin_id].w
            fb.velocities[1] = fb.v_mean
        elseif (loc_np == 1)
            # no merge, just copy over
            indices = fb.particle_indices
            fb.weights[1] = particles[indices[1]].w
            fb.velocities[1] = particles[indices[1]].v
        end
    end

    # --- 3. Compute post-merge statistics ---
    # mean and total weight are conserved by construction
    octree.v_var_post = SVector{3,Float64}(0.0, 0.0, 0.0)

    @inbounds for bin_id in 1:octree.Nbins
        if octree.bins[bin_id].np > 0
            fb = octree.full_bins[bin_id]
            octree.v_var_post += fb.weights[1] * (fb.velocities[1] - octree.v_mean).^2
        end
    end

    octree.v_var_post = octree.v_var_post / w_total

    # --- 4. Compute scaling factors and scale data in full_bins ---
    scaling_factor_v = variance_scaling(octree.v_var_before, octree.v_var_post)

    @inbounds for bin_id in 1:octree.Nbins
        if octree.bins[bin_id].np > 0
            fb = octree.full_bins[bin_id]
            fb.velocities[1] = (fb.velocities[1] - octree.v_mean) .* scaling_factor_v + octree.v_mean
        end
    end

    # --- 5. write back to particles ---
    curr_particle_index = write_back_to_particles!(octree, particles, pia, cell, species)

    # --- 6. delete extra particles ---
    @inbounds old_count = pia.indexer[cell, species].n_local
    n_particles_to_delete = old_count - curr_particle_index

    @inbounds if !(cell == size(pia.indexer)[1]) || (n_particles_to_delete > pia.indexer[cell, species].n_group2)
        pia.contiguous[species] = false
    end

    for _ in 1:n_particles_to_delete
        delete_particle_end!(particles, pia, cell, species)
    end
end

"""
    compute_new_particles!(rng, octree::OctreeMerge{D,1}, particles::ParticleVector{D}, pia, cell, species, grid::Grid1DUniform)

Compute post-merge particles with particles based on octree bin properties, N:1 merging in each bin; placing out-of-domain particles back into the domain..
It computes properties of all particles. And after the first merge, it scales the particles based on these properties to achieve conservation.

# Positional arguments:
* `rng`: the random number generator instance
* `octree`: the `OctreeMerge` instance
* `particles`: the `ParticleVector` instance of the particles to be merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the grid cell in which particles are being merged
* `species`: the index of the species being merged
"""
function compute_new_particles!(rng, octree::OctreeMerge{D,1}, particles::ParticleVector{D}, pia, cell, species, grid::Grid1DUniform) where D
    # --- 1. Compute statistics before merging ---
    # bin-level properties cannot be used here: compute_bin_props! skips bins with np <= M,
    # so full_bins[bin_id].v_mean is stale for single-particle bins when M == 1
    w_total = 0.0
    octree.v_mean = SVector{3,Float64}(0.0, 0.0, 0.0)
    octree.x_mean = zero(SVector{D,Float64})

    @inbounds for bin_id in 1:octree.Nbins
        bs, be = octree.bin_start[bin_id], octree.bin_end[bin_id]
        for i in bs:be
            p_i = particles[octree.particle_indexes_sorted[i]]
            w_total += p_i.w
            octree.v_mean = octree.v_mean + p_i.w * p_i.v
            octree.x_mean = octree.x_mean + p_i.w * p_i.x
        end
    end

    octree.v_mean = octree.v_mean / w_total
    octree.x_mean = octree.x_mean / w_total

    octree.v_var_before = SVector{3,Float64}(0.0, 0.0, 0.0)
    octree.x_var_before = zero(SVector{D,Float64})

    @inbounds for bin_id in 1:octree.Nbins
        bs, be = octree.bin_start[bin_id], octree.bin_end[bin_id]
        for i in bs:be
            p_i = particles[octree.particle_indexes_sorted[i]]
            octree.v_var_before = octree.v_var_before + p_i.w * (p_i.v - octree.v_mean).^2
            octree.x_var_before = octree.x_var_before + p_i.w * (p_i.x - octree.x_mean).^2
        end
    end

    octree.v_var_before = octree.v_var_before / w_total
    octree.x_var_before = octree.x_var_before / w_total

    # --- 2. Generate post-merge particle data into full_bins ---
    @inbounds for bin_id in 1:octree.Nbins
        loc_np = octree.bins[bin_id].np
        fb = octree.full_bins[bin_id]
        if (loc_np > 1)
            # N:1 merge
            fb.weights[1] = octree.bins[bin_id].w
            fb.velocities[1] = fb.v_mean
            fb.positions[1] = fb.x_mean
        elseif (loc_np == 1)
            # no merge, just copy over
            indices = fb.particle_indices
            fb.weights[1] = particles[indices[1]].w
            fb.velocities[1] = particles[indices[1]].v
            fb.positions[1] = particles[indices[1]].x
        end
    end

    # --- 3. Compute post-merge statistics ---
    # mean and total weight are conserved by construction
    octree.v_var_post = SVector{3,Float64}(0.0, 0.0, 0.0)
    octree.x_var_post = zero(SVector{D,Float64})

    @inbounds for bin_id in 1:octree.Nbins
        if octree.bins[bin_id].np > 0
            fb = octree.full_bins[bin_id]
            octree.v_var_post += fb.weights[1] * (fb.velocities[1] - octree.v_mean).^2
            octree.x_var_post += fb.weights[1] * (fb.positions[1] - octree.x_mean).^2
        end
    end

    octree.v_var_post = octree.v_var_post / w_total
    octree.x_var_post = octree.x_var_post / w_total

    # --- 4. Compute scaling factors and scale data in full_bins ---
    scaling_factor_v = variance_scaling(octree.v_var_before, octree.v_var_post)
    scaling_factor_x = variance_scaling(octree.x_var_before, octree.x_var_post)

    @inbounds for bin_id in 1:octree.Nbins
        if octree.bins[bin_id].np > 0
            fb = octree.full_bins[bin_id]
            fb.velocities[1] = (fb.velocities[1] - octree.v_mean) .* scaling_factor_v + octree.v_mean
            fb.positions[1] = (fb.positions[1] - octree.x_mean) .* scaling_factor_x + octree.x_mean
        end
    end

    # --- 5. write back to particles ---
    curr_particle_index = write_back_to_particles!(octree, particles, pia, cell, species, grid)

    # --- 6. delete extra particles ---
    @inbounds old_count = pia.indexer[cell, species].n_local
    n_particles_to_delete = old_count - curr_particle_index

    @inbounds if !(cell == size(pia.indexer)[1]) || (n_particles_to_delete > pia.indexer[cell, species].n_group2)
        pia.contiguous[species] = false
    end

    for _ in 1:n_particles_to_delete
        delete_particle_end!(particles, pia, cell, species)
    end
end

"""
    write_back_to_particles!(octree::OctreeMerge{D,M}, particles::ParticleVector{D}, pia, cell, species) where {D,M}

Write back merged particles to the particle array (0D case - position of particles not set).
It will only be used in the function compute_new_particles! without a grid argument, which means no boundary handling is needed.

# Positional arguments:
* `particles`: the `ParticleVector` instance containing the particles to be merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the grid cell in which particles are being merged
* `species`: the index of the species being merged
* `octree`: the `OctreeMerge` instance

# Returns
The current particle index after writing.
"""
function write_back_to_particles!(octree::OctreeMerge{D,M}, particles::ParticleVector{D}, pia, cell, species) where {D,M}
    curr_particle_index = 0
    @inbounds indexer = pia.indexer[cell,species]
    @inbounds for bin_id in 1:octree.Nbins
        loc_np = octree.bins[bin_id].np
        full_bin = octree.full_bins[bin_id]
        if (loc_np >= M)
            @inbounds for j in 1:M
                i = map_cont_index(indexer, curr_particle_index)
                curr_particle_index += 1
                particles[i].w = full_bin.weights[j]
                particles[i].v = full_bin.velocities[j]
            end
        elseif (loc_np > 0)
            @inbounds for j in 1:loc_np
                i = map_cont_index(indexer, curr_particle_index)
                curr_particle_index += 1
                particles[i].w = full_bin.weights[j]
                particles[i].v = full_bin.velocities[j]
            end
        end
    end
    return curr_particle_index
end


"""
    write_back_to_particles!(octree::OctreeMerge{D,M}, particles::ParticleVector{D}, pia, cell, species, grid::Grid1DUniform) where N

Write back merged particles to the particle array with boundary clamping (1D case).
Particles positioned outside the grid boundaries (`grid.min_x`, `grid.max_x`) are clamped to the boundary.
It will only be used in the function compute_new_particles! with grid::Grid1DUniform argument.

# Positional arguments:
* `particles`: the `ParticleVector` instance containing the particles to be merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the grid cell in which particles are being merged
* `species`: the index of the species being merged
* `octree`: the `OctreeMerge` instance
* `grid`: the `Grid1DUniform` instance

# Returns
The current particle index after writing.
"""
function write_back_to_particles!(octree::OctreeMerge{D,M}, particles::ParticleVector{D}, pia, cell, species, grid::Grid1DUniform) where {D,M}
    curr_particle_index = 0
    @inbounds indexer = pia.indexer[cell,species]
    @inbounds for bin_id in 1:octree.Nbins
        loc_np = octree.bins[bin_id].np
        full_bin = octree.full_bins[bin_id]
        if (loc_np >= M) 
            @inbounds for j in 1:M
                i = map_cont_index(indexer, curr_particle_index)
                curr_particle_index += 1
                particles[i].w = full_bin.weights[j]
                particles[i].v = full_bin.velocities[j]
                val = full_bin.positions[j]
                clamped_x = clamp(val[1], grid.min_x, grid.max_x)

                particles[i].x = set_x(val, clamped_x)
            end
        elseif (loc_np > 0)
            @inbounds for j in 1:loc_np
                i = map_cont_index(indexer, curr_particle_index)
                curr_particle_index += 1
                particles[i].w = full_bin.weights[j]
                particles[i].v = full_bin.velocities[j]
                particles[i].x = full_bin.positions[j]
            end
        end
    end
    return curr_particle_index
end


"""
    init_octree!(octree::OctreeMerge{D,M}, particles::ParticleVector{D}, pia, cell, species)

Initialize the top bin in an octree by copying particle indices and setting bin bounds.

# Positional arguments
* `octree`: the `OctreeMerge` instance
* `particles`: the `ParticleVector` instance containing the particles to be merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the grid cell in which particles are being merged
* `species`: the index of the species being merged
"""
function init_octree!(octree::OctreeMerge{D,M}, particles::ParticleVector{D}, pia, cell, species) where {D,M}
    octree.Nbins = 1
    @inbounds octree.particle_indexes_sorted[1:pia.indexer[cell,species].n_group1] = pia.indexer[cell,species].start1:pia.indexer[cell,species].end1

    @inbounds if (pia.indexer[cell,species].n_group2 > 0)
        @inbounds octree.particle_indexes_sorted[pia.indexer[cell,species].n_group1+1:pia.indexer[cell,species].n_local] = pia.indexer[cell,species].start2:pia.indexer[cell,species].end2
    end

    @inbounds octree.n_particles = pia.indexer[cell,species].n_local
    @inbounds octree.bins[1].depth = 0
    @inbounds octree.bin_start[1] = 1
    @inbounds octree.bin_end[1] = octree.n_particles

    @inbounds octree.bins[1].np = octree.n_particles
    @inbounds octree.bins[1].w = 1e50  # just do a large estimate, we will be splitting this bin anyway

    octree.total_post_merge_np = get_bin_post_merge_np(octree, 1)
    @inbounds if (octree.bins[1].np > M) && (octree.max_depth > 0)
        @inbounds octree.bins[1].can_be_refined = true
    else
        @inbounds octree.bins[1].can_be_refined = false
    end

    if (octree.init_bin_bounds == OctreeInitBinC)
        @inbounds octree.bins[1].v_min = SVector{3, Float64}(-299_792_458.0, -299_792_458.0, -299_792_458.0)  # speed of light
        @inbounds octree.bins[1].v_max = SVector{3, Float64}(299_792_458.0, 299_792_458.0, 299_792_458.0)
    else
        bin_bounds_recompute!(octree, 1, 1, octree.n_particles, particles)
        if (octree.init_bin_bounds == OctreeInitBinMinMaxVelSym)
            @inbounds maxvx = max(abs(octree.bins[1].v_min[1]), abs(octree.bins[1].v_max[1]))
            @inbounds maxvy = max(abs(octree.bins[1].v_min[2]), abs(octree.bins[1].v_max[2]))
            @inbounds maxvz = max(abs(octree.bins[1].v_min[3]), abs(octree.bins[1].v_max[3]))

            @inbounds octree.bins[1].v_min = SVector{3, Float64}(-maxvx, -maxvy, -maxvz)
            @inbounds octree.bins[1].v_max = SVector{3, Float64}(maxvx, maxvy, maxvz)
        end
    end
end


"""
    compute_octree!(octree::OctreeMerge{D,M}, particles::ParticleVector{D}, target_np)

Perform refinement of an octree for N:2 merging until target number of particles reached or nothing left to refine.

# Positional arguments
* `octree`: the `OctreeMerge` instance
* `particles`: the `ParticleVector` instance containing the particles to be merged
* `target_np`: the target post-merge number of particles; the post-merge number of particles will not exceed this value
    but may be not exactly equal to it 
"""
function compute_octree!(octree::OctreeMerge{D,M}, particles::ParticleVector{D}, target_np) where {D,M}
    while true
        refine_id = -1
        max_w = -1.0
        @inbounds for bin_id in 1:octree.Nbins
            bin = octree.bins[bin_id]
            # find bin with largest number density, needs to be refine-able
            if ((bin.w > max_w) && bin.can_be_refined)
               max_w = bin.w
               refine_id = bin_id
           end
        end

        if refine_id == -1
            # found no bin to refine, i.e. ran out of particles
            break
        elseif (octree.total_post_merge_np + 7*M > target_np)
            # refining a bin can produce up to 8*M particles, so we don't do it if threshold exceeded
            # but if we have a bin it has potentially M particles already, so refinement increases count
            # only by 7*M
            break
        else
            split_bin!(octree, refine_id, particles)
        end

        if octree.Nbins + 7 > octree.max_Nbins
            # reach max number of bins possible
            break
        end
    end

    if octree.Nbins == 1
        @inbounds octree.bins[1].w = 0.0
        
        @inbounds for ii in 1:octree.n_particles
            octree.bins[1].w += particles[octree.particle_indexes_sorted[ii]].w
        end
    end
    
    for bin_id in 1:octree.Nbins
        compute_bin_props!(octree, bin_id, particles)
    end
end

"""
    merge_octree!(rng, octree::OctreeMerge{D,M}, particles::ParticleVector{D}, pia, cell, species, target_np)

Perform octree N:M merging without checking whether particle positions end up outside of the simulation domain.

# Positional arguments
* `rng`: the random number generator instance
* `octree`: the `OctreeMerge` instance
* `particles`: the `ParticleVector` instance containing the particles to be merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the grid cell in which particles are being merged
* `species`: the index of the species being merged
* `target_np`: the target post-merge number of particles; the post-merge number of particles will not exceed this value
    but may be not exactly equal to it

# References
* R.S. Martin, J.-L. Cambier, Octree particle management for DSMC and PIC simulations.
    [J. Comput. Phys., 2016](https://doi.org/10.1016/j.jcp.2016.01.020).
"""
function merge_octree!(rng, octree::OctreeMerge{D,M}, particles::ParticleVector{D}, pia, cell, species, target_np) where {D,M}
    clear_octree!(octree)
    @inbounds resize_octree_buffers!(octree, pia.indexer[cell,species].n_local)
    init_octree!(octree, particles, pia, cell, species)
    compute_octree!(octree, particles, target_np)
    compute_new_particles!(rng, octree, particles, pia, cell, species)
end

"""
    merge_octree!(rng, octree::OctreeMerge{D,M}, particles::ParticleVector{D}, pia, cell, species, target_np, grid::Grid1DUniform)

Perform octree N:M merging, checking whether particle positions end up outside of the simulation domain, and placing them back into the domain
if needed.

# Positional arguments
* `rng`: the random number generator instance
* `octree`: the `OctreeMerge` instance
* `particles`: the `ParticleVector` instance containing the particles to be merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the grid cell in which particles are being merged
* `species`: the index of the species being merged
* `target_np`: the target post-merge number of particles; the post-merge number of particles will not exceed this value
    but may be not exactly equal to it
* `grid`: the `Grid1DUniform` grid

* R.S. Martin, J.-L. Cambier, Octree particle management for DSMC and PIC simulations.
    [J. Comput. Phys., 2016](https://doi.org/10.1016/j.jcp.2016.01.020).
"""
function merge_octree!(rng, octree::OctreeMerge{D,M}, particles::ParticleVector{D}, pia, cell, species, target_np, grid::Grid1DUniform) where {D,M}
    clear_octree!(octree)
    @inbounds resize_octree_buffers!(octree, pia.indexer[cell,species].n_local)
    init_octree!(octree, particles, pia, cell, species)
    compute_octree!(octree, particles, target_np)
    compute_new_particles!(rng, octree, particles, pia, cell, species, grid)
end

end