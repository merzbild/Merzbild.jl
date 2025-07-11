using StaticArrays

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

    depth::Int32
    can_be_refined::Bool
end

"""
    OctreeFullCell

Struct holding computed bin properties required to merge the particles in a bin.

# Fields
* `v_mean`: mean velocity of particles in cell
* `v_std_sq`: variance of velocity of particles in cell
* `x_mean`: mean position of particles in cell
* `x_std_sq`: variance of position of particles in cell
* `particle_index1`: index of first particle in the cell (first as in the first particle that is found to
    belong to the cell)
* `particle_index1`: index of second particle in the cell (second as in the second particle that is found to
    belong to the cell)
* `w1`: the post-merge weight to assign to the first particle in the cell
* `w2`: the post-merge weight to assign to the second particle in the cell
* `v1`: the post-merge velocity to assign to the first particle in the cell
* `v2`: the post-merge velocity to assign to the second particle in the cell
* `x1`: the post-merge position to assign to the first particle in the cell
* `x2`: the post-merge position to assign to the second particle in the cell
"""
mutable struct OctreeFullCell
    v_mean::SVector{3,Float64}
    v_std_sq::SVector{3,Float64}
    x_mean::SVector{3,Float64}
    x_std_sq::SVector{3,Float64}
    particle_index1::Int64
    particle_index2::Int64

    w1::Float64
    w2::Float64
    v1::SVector{3,Float64}  # these are for the post-merge quantities
    v2::SVector{3,Float64}
    x1::SVector{3,Float64}
    x2::SVector{3,Float64}
end

# struct for N:2 merge
"""
    OctreeN2Merge

Struct for N:2 Octree merging.

# Fields
* `max_Nbins`: maximum possible number of bins
* `Nbins`: number of bins currently used
* `bins`: Vector of `OctreeCell` instances used to compute the properties required for bin refinement
* `full_bins`: Vector of `OctreeFullCell` instances used to compute the post-merge particles in each bin
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
* `direction_vec`: used to store randomly sampled direction signs
* `init_bin_bounds`: enum of `OctreeInitBin` type defining how the bounds of the top-level bin are set
* `max_depth`: maximum allowed depth of a bin
* `total_post_merge_np`: used to keep track of number of post-merge particles
"""
mutable struct OctreeN2Merge
    max_Nbins::Int32
    Nbins::Int32  # actual bins computed
    bins::Vector{OctreeCell}
    full_bins::Vector{OctreeFullCell}
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

    direction_vec::SVector{3,Float64}

    init_bin_bounds::OctreeInitBin

    max_depth::Int32
    total_post_merge_np::Int64 # used to keep track of number of post-merge particles
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
    fill_full_bins(Nbins)

Fill the octree bins full structs with zero data, used as a utility function for initialization.

Positional arguments:
* `Nbins`: number of `OctreeFullCell` bins to create

Returns:
An array of `Nbins` `OctreeFullCell` instances filled with zeros.
"""
function fill_full_bins(Nbins)
    return [OctreeFullCell(SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0),
                           SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0),
                           0, 0, 0.0, 0.0, 
                           SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0),
                           SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0)) for i in 1:Nbins]
end

"""
    OctreeN2Merge(split::OctreeBinSplit; init_bin_bounds=OctreeInitBinMinMaxVel, bin_bounds_compute=OctreeBinBoundsInherit,
              max_Nbins=4096, max_depth=10)
    
Create an Octree N:2 merging instance.

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
`OctreeN2Merge` instance with everything set to 0.
"""
OctreeN2Merge(split::OctreeBinSplit; init_bin_bounds=OctreeInitBinMinMaxVel, bin_bounds_compute=OctreeBinBoundsInherit,
              max_Nbins=4096, max_depth=10) = OctreeN2Merge(max_Nbins, 0, fill_bins(max_Nbins), fill_full_bins(max_Nbins), 0,
                                                            zeros(max_Nbins), zeros(max_Nbins),  # bin_start, bin_end
                                                            zeros(8192), zeros(8192), zeros(8192),
                                                            MVector{8, Int64}(0, 0, 0, 0, 0, 0, 0, 0),
                                                            MVector{8, Int64}(0, 0, 0, 0, 0, 0, 0, 0),
                                                            MVector{8, Int64}(0, 0, 0, 0, 0, 0, 0, 0),
                                                            MVector{8, Float64}(0, 0, 0, 0, 0, 0, 0, 0),
                                                            bin_bounds_compute, split,
                                                            SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0),
                                                            SVector{3,Float64}(0.0, 0.0, 0.0),
                                                            init_bin_bounds, max_depth, 0)

"""
    clear_octree!(octree)
    
Reset octree before doing a new merge.

# Positional arguments:
* `octree`: the `OctreeN2Merge` instance
"""
function clear_octree!(octree)
    octree.Nbins = 0
end

"""
    resize_octree_buffers!(octree, n_particles)

Check and resize octree buffers if needed to accommodate a larger number of particles. The size of the buffers
is set to `n_particles + DELTA_PARTICLES` if the their sizes are smaller than `n_particles`.

# Positional arguments:
* `octree`: the `OctreeN2Merge` instance
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

    oct = 1
    if (particle_v[1] > v_middle[1])
        oct += 1
    end
    if (particle_v[2] > v_middle[2])
        oct += 2
    end
    if (particle_v[3] > v_middle[3])
        oct += 4
    end
    return oct
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
* `octree`: the `OctreeN2Merge` instance
* `bin_id`: the index of the bin for which the velocity bounds are being recomputed
* `v_min_parent`: the velocity vector of the lower bound of the velocities of the parent bin
* `v_max_parent`: the velocity vector of the upper bound of the velocities of the parent bin
* `v_middle`: the velocity across which the split is being performed
* `octant`: the octant of the parent bin to which this bin corresponds
"""
function bin_bounds_inherit!(octree, bin_id, v_min_parent, v_max_parent, v_middle, octant)

    
    if octant == 1
        octree.bins[bin_id].v_min = v_min_parent
        octree.bins[bin_id].v_max = v_middle
    elseif octant == 2
        octree.bins[bin_id].v_min = SVector{3,Float64}(v_middle[1], v_min_parent[2], v_min_parent[3])
        octree.bins[bin_id].v_max = SVector{3,Float64}(v_max_parent[1], v_middle[2], v_middle[3])
    elseif octant == 3
        octree.bins[bin_id].v_min = SVector{3,Float64}(v_min_parent[1], v_middle[2], v_min_parent[3])
        octree.bins[bin_id].v_max = SVector{3,Float64}(v_middle[1], v_max_parent[2], v_middle[3])
    elseif octant == 4
        octree.bins[bin_id].v_min = SVector{3,Float64}(v_middle[1], v_middle[2], v_min_parent[3])
        octree.bins[bin_id].v_max = SVector{3,Float64}(v_max_parent[1], v_max_parent[2], v_middle[3])
    elseif octant == 5
        octree.bins[bin_id].v_min = SVector{3,Float64}(v_min_parent[1], v_min_parent[2], v_middle[3])
        octree.bins[bin_id].v_max = SVector{3,Float64}(v_middle[1], v_middle[2], v_max_parent[3])
    elseif octant == 6
        octree.bins[bin_id].v_min = SVector{3,Float64}(v_middle[1], v_min_parent[2], v_middle[3])
        octree.bins[bin_id].v_max = SVector{3,Float64}(v_max_parent[1], v_middle[2], v_max_parent[3])
    elseif octant == 7
        octree.bins[bin_id].v_min = SVector{3,Float64}(v_min_parent[1], v_middle[2], v_middle[3])
        octree.bins[bin_id].v_max = SVector{3,Float64}(v_middle[1], v_max_parent[2], v_max_parent[3])
    else
        octree.bins[bin_id].v_min = v_middle
        octree.bins[bin_id].v_max = v_max_parent
    end
end

"""
    bin_bounds_recompute!(octree, bin_id, bs, be, particles)

Recompute bin bounds based on particle velocities by setting them to the smallest and largest
velocities of the particles in each velocity direction.

# Positional arguments:
* `octree`: the `OctreeN2Merge` instance
* `bin_id`: the index of the bin for which the velocity bounds are being recomputed
* `bs`: index of the first particle in bin
* `be`: index of the last particle in bin
* `particles`: the `ParticleVector` instance of the particles to be merged
"""
function bin_bounds_recompute!(octree, bin_id, bs, be, particles)
    # compute bin bounds based on the particles in the bin
    minvx = 9_299_792_458.0  # speed of light + 9e9
    minvy = 9_299_792_458.0
    minvz = 9_299_792_458.0

    maxvx = -9_299_792_458.0
    maxvy = -9_299_792_458.0
    maxvz = -9_299_792_458.0

    # for pi in octree.particle_indexes_sorted[bs:be]
    for i in bs:be
        pi = octree.particle_indexes_sorted[i]
        if (particles[pi].v[1] < minvx)
            minvx = particles[pi].v[1]
        elseif (particles[pi].v[1] > maxvx)
            maxvx = particles[pi].v[1]
        end

        if (particles[pi].v[2] < minvy)
            minvy = particles[pi].v[2]
        elseif (particles[pi].v[2] > maxvy)
            maxvy = particles[pi].v[2]
        end

        if (particles[pi].v[3] < minvz)
            minvz = particles[pi].v[3]
        elseif (particles[pi].v[3] > maxvz)
            maxvz = particles[pi].v[3]
        end

        octree.bins[bin_id].v_min = SVector{3, Float64}(minvx, minvy, minvz)
        octree.bins[bin_id].v_max = SVector{3, Float64}(maxvx, maxvy, maxvz)
    end
end

"""
    compute_v_mean!(octree, bs, be, particles)

Compute mean velocity of particles in a bin.

# Positional arguments:
* `octree`: the `OctreeN2Merge` instance
* `bs`: index of the first particle in bin
* `be`: index of the last particle in bin
* `particles`: the `ParticleVector` instance of the particles to be merged
"""
@inline function compute_v_mean!(octree, bs, be, particles)
    n_tot = 0.0
    for i in bs:be
        pi = octree.particle_indexes_sorted[i]
        n_tot += particles[pi].w
        octree.vel_middle = octree.vel_middle + particles[pi].w * particles[pi].v
    end
    octree.vel_middle = octree.vel_middle / n_tot
end

"""
    compute_v_median!(octree, bs, be, particles)

Compute median velocity of particles in a bin. NOTE: allocates memory and is probably not fully correct!

# Positional arguments:
* `octree`: the `OctreeN2Merge` instance
* `bs`: index of the first particle in bin
* `be`: index of the last particle in bin
* `particles`: the `ParticleVector` instance of the particles to be merged
"""
function compute_v_median!(octree, bs, be, particles)
    # https://rcoh.me/posts/linear-time-median-finding/
    # https://en.wikipedia.org/wiki/Weighted_median
    w_vec = [particles[octree.particle_indexes_sorted[i]].w for i in bs:be]
    w_vec = w_vec ./ sum(w_vec)

    vx_vec = [particles[octree.particle_indexes_sorted[i]].v[1] for i in bs:be]
    vy_vec = [particles[octree.particle_indexes_sorted[i]].v[2] for i in bs:be]
    vz_vec = [particles[octree.particle_indexes_sorted[i]].v[3] for i in bs:be]

    octree.vel_middle = SVector{3, Float64}(weighted_percentile_interpolated(vx_vec, w_vec),
                                            weighted_percentile_interpolated(vy_vec, w_vec),
                                            weighted_percentile_interpolated(vz_vec, w_vec))                                  
end

function bin_bounds_recompute_and_v_mean!(octree, bin_id, bs, be, particles)
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
    split_bin!(octree, bin_id, particles)

Sort particles into sub-bins of a bin with index `bin_id` (by splitting it into octants),
keeping track of which sub-bins particles end up in.
Also sets the velocity bounds of the new bins.

# Positional arguments:
* `octree`: the `OctreeN2Merge` instance
* `bin_id`: octree bin index
* `particles`: the `ParticleVector` instance of the particles to be merged
"""
function split_bin!(octree, bin_id, particles)
    octree.particle_in_bin_counter .= 0 # reset counter
    octree.ndens_counter .= 0.0
    octree.nonempty_bins .= 0
    octree.nonempty_counter .= 0
    current_depth = octree.bins[bin_id].depth

    n_nonempty_bins = 0
    bs = octree.bin_start[bin_id]
    be = octree.bin_end[bin_id]

    # if we compute bin bounds based on particle velocities, we need to actually iterate through particles
    # and compute
    # if not, we can do it later directly for the sub-bins
    if (octree.bin_bounds_compute == OctreeBinBoundsRecompute)
        bin_bounds_recompute!(octree, bin_id, bs, be, particles)
    end

    if (octree.split == OctreeBinMidSplit)
        octree.vel_middle = 0.5 * (octree.bins[bin_id].v_min + octree.bins[bin_id].v_max)
    elseif (octree.split == OctreeBinMeanSplit)
        compute_v_mean!(octree, bs, be, particles) # octree.vel_middle = octree.bins[bin_id].v_mean
    elseif (octree.split == OctreeBinMedianSplit)
        compute_v_median!(octree, bs, be, particles)
    end

    for i in bs:be
        pi = octree.particle_indexes_sorted[i]
        oct = compute_octant(particles[pi].v, octree.vel_middle)
        octree.particle_in_bin_counter[oct] += 1
        octree.particle_octants[i-bs+1] = oct
        octree.ndens_counter[oct] += particles[pi].w
    end

    n_eb = 0
    if (octree.particle_in_bin_counter[1] > 0)
        n_nonempty_bins += 1
        n_eb += 1
        octree.nonempty_counter[n_eb] = octree.particle_in_bin_counter[1]
        octree.nonempty_bins[n_eb] = 1
    end

    for i in 2:8
        if (octree.particle_in_bin_counter[i] > 0)
            n_nonempty_bins += 1
            n_eb += 1
            octree.nonempty_counter[n_eb] = octree.particle_in_bin_counter[i]
            octree.nonempty_bins[n_eb] = i
        end

        octree.particle_in_bin_counter[i] += octree.particle_in_bin_counter[i-1]
    end

    # first bin - we change nothing for the start
    octree.bin_end[bin_id] = octree.bin_start[bin_id] + octree.nonempty_counter[1] - 1
    
    # new bins will point to the sorted particles, but are not contiguous
    # first bin is in the old place and the new ones are tacked on
    for i in 2:n_nonempty_bins
        bi = get_new_bin_id(i, bin_id, octree.Nbins)
        bim1 = get_new_bin_id(i-1, bin_id, octree.Nbins)
        octree.bin_start[bi] = octree.bin_end[bim1] + 1
        octree.bin_end[bi] = octree.bin_start[bi] + octree.nonempty_counter[i] - 1
    end

    # we had a bin that would've produced 2 particles
    # now we replaced ith with n_nonempty_bins bins that each produce 1 or 2 particles
    octree.total_post_merge_np -= 2
    if (octree.bin_bounds_compute == OctreeBinBoundsInherit)
        octree.v_min_parent = octree.bins[bin_id].v_min
        octree.v_max_parent = octree.bins[bin_id].v_max

        # iterate over non-empty bins and inherit parent bin bounds + split around middle velocity
        for i in 1:n_nonempty_bins
            bi = get_new_bin_id(i, bin_id, octree.Nbins)
            bin_bounds_inherit!(octree, bi,
                                octree.v_min_parent, octree.v_max_parent,
                                octree.vel_middle, octree.nonempty_bins[i])
            octree.bins[bi].np = octree.nonempty_counter[i]
            octree.bins[bi].w = octree.ndens_counter[octree.nonempty_bins[i]]
            octree.bins[bi].depth = current_depth + 1

            octree.total_post_merge_np += get_bin_post_merge_np(octree, bi)
            # octree.bins[bin_id + i - 1].post_merge_np = get_bin_post_merge_np(octree, bin_id + i - 1)
            if (octree.bins[bi].np > 2) && (octree.bins[bi].depth < octree.max_depth)
                octree.bins[bi].can_be_refined = true
            else
                octree.bins[bi].can_be_refined = false
            end
        end
    else
        # still need to fill out info on number of particles and total weight
        # will recompute bin bounds if we do next round of refinement
        for i in 1:n_nonempty_bins
            bi = get_new_bin_id(i, bin_id, octree.Nbins)
            octree.bins[bi].np = octree.nonempty_counter[i]
            octree.bins[bi].w = octree.ndens_counter[octree.nonempty_bins[i]]
            octree.bins[bi].depth = current_depth + 1

            # we had a bin that would've produced 2 particles
            # now we replaced ith with n_nonempty_bins bins that each produce 1 or 2 particles
            octree.total_post_merge_np += get_bin_post_merge_np(octree, bi)
            if (octree.bins[bi].np > 2) && (octree.bins[bi].depth < octree.max_depth)
                octree.bins[bi].can_be_refined = true
            else
                octree.bins[bi].can_be_refined = false
            end
        end
    end

    for i in bs:be
        pi = octree.particle_indexes_sorted[i]
        j = octree.particle_octants[i - bs + 1]
        octree.particles_sort_output[octree.particle_in_bin_counter[j]] = pi
        octree.particle_in_bin_counter[j] -= 1
    end

    # write sorted indices
    # octree.particle_indexes_sorted[bs:be] = octree.particles_sort_output[1:be-bs+1]
    for i in bs:be
        octree.particle_indexes_sorted[i] = octree.particles_sort_output[i-bs+1]
    end

    octree.Nbins += n_nonempty_bins - 1
end

"""
    compute_bin_props!(octree, bin_id, particles)

Compute properties in a bin required for merging: total computational weight, mean velocity and position,
standard deviation of particle velocities and positions.

# Positional arguments:
* `octree`: the `OctreeN2Merge` instance
* `bin_id`: octree bin index
* `particles`: the `ParticleVector` instance of the particles to be merged
"""
function compute_bin_props!(octree, bin_id, particles)
    bs = octree.bin_start[bin_id]
    be = octree.bin_end[bin_id]

    octree.full_bins[bin_id].v_mean = SVector{3, Float64}(0.0, 0.0, 0.0)
    octree.full_bins[bin_id].v_std_sq = SVector{3, Float64}(0.0, 0.0, 0.0)

    octree.full_bins[bin_id].x_mean = SVector{3, Float64}(0.0, 0.0, 0.0)
    octree.full_bins[bin_id].x_std_sq = SVector{3, Float64}(0.0, 0.0, 0.0)

    # store indices of the first 1/2 particles in octree bin so that we
    # have somewhere to write post-merge data
    # but also so that we don't do unnecessary merging (2:2, 1:2)
    if (octree.bins[bin_id].np == 1)
        octree.full_bins[bin_id].particle_index1 = octree.particle_indexes_sorted[bs]
    elseif (octree.bins[bin_id].np >= 2)
        octree.full_bins[bin_id].particle_index1 = octree.particle_indexes_sorted[bs]
        octree.full_bins[bin_id].particle_index2 = octree.particle_indexes_sorted[bs + 1]
    end

    # if only 2 or fewer particles in bin then we don't need to compute any properties
    if (octree.bins[bin_id].np <= 2)
        return
    end

    # for pi in octree.particle_indexes_sorted[bs:be]
    for i in bs:be
        pi = octree.particle_indexes_sorted[i]
        octree.full_bins[bin_id].v_mean = octree.full_bins[bin_id].v_mean + particles[pi].w * particles[pi].v
        octree.full_bins[bin_id].x_mean = octree.full_bins[bin_id].x_mean + particles[pi].w * particles[pi].x
    end
    octree.full_bins[bin_id].v_mean = octree.full_bins[bin_id].v_mean / octree.bins[bin_id].w
    octree.full_bins[bin_id].x_mean = octree.full_bins[bin_id].x_mean / octree.bins[bin_id].w

    for i in bs:be
        pi = octree.particle_indexes_sorted[i]
        octree.full_bins[bin_id].v_std_sq = octree.full_bins[bin_id].v_std_sq +
                                            particles[pi].w * (particles[pi].v - octree.full_bins[bin_id].v_mean).^2
        octree.full_bins[bin_id].x_std_sq = octree.full_bins[bin_id].x_std_sq +
                                            particles[pi].w * (particles[pi].x - octree.full_bins[bin_id].x_mean).^2
    end
    octree.full_bins[bin_id].v_std_sq = octree.full_bins[bin_id].v_std_sq / octree.bins[bin_id].w
    octree.full_bins[bin_id].x_std_sq = octree.full_bins[bin_id].x_std_sq / octree.bins[bin_id].w
end

"""
    get_bin_post_merge_np(octree, bin_id)

Get number of post-merge particles in a bin: 2 if the number of particles in the bin is >= 2, otherwise
the number of particles in the bin (0 or 1) is returned.

# Positional arguments:
* `octree`: the `OctreeN2Merge` instance
* `bin_id`: octree bin index

# Returns
The number of post-merge particles in a single octree bin (0, 1, or 2).
"""
@inline function get_bin_post_merge_np(octree, bin_id)
    # how many particles will we get after merging in the bin:
    # 2 if np >= 2
    # np otherwise (0 or 1)
    return octree.bins[bin_id].np >= 2 ? 2 : octree.bins[bin_id].np
end

"""
    compute_new_particles!(rng, octree::OctreeN2Merge, particles, pia, cell, species)

Compute post-merge particles particles based on octree bin properties without checking particle locations.
So particles may end up outside of the domain.

# Positional arguments:
* `rng`: the random number generator instance
* `octree`: the `OctreeN2Merge` instance
* `particles`: the `ParticleVector` instance of the particles to be merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the cell index
* `species`: the species index
"""
function compute_new_particles!(rng, octree::OctreeN2Merge, particles, pia, cell, species)
    # given computed Octree, create new particles instead of the old ones
    
    for bin_id in 1:octree.Nbins
        if (octree.bins[bin_id].np > 2)
            octree.full_bins[bin_id].w1 = 0.5 * octree.bins[bin_id].w
            octree.full_bins[bin_id].w2 = octree.full_bins[bin_id].w1

            octree.full_bins[bin_id].v_std_sq = sqrt.(octree.full_bins[bin_id].v_std_sq)
            octree.full_bins[bin_id].x_std_sq = sqrt.(octree.full_bins[bin_id].x_std_sq)
            
            octree.direction_vec = @SVector rand(rng, direction_signs, 3)
            octree.full_bins[bin_id].v1 = octree.full_bins[bin_id].v_mean + octree.direction_vec .* octree.full_bins[bin_id].v_std_sq
            octree.full_bins[bin_id].v2 = octree.full_bins[bin_id].v_mean - octree.direction_vec .* octree.full_bins[bin_id].v_std_sq

            octree.direction_vec = @SVector rand(rng, direction_signs, 3)
            octree.full_bins[bin_id].x1 = octree.full_bins[bin_id].x_mean + octree.direction_vec .* octree.full_bins[bin_id].x_std_sq
            octree.full_bins[bin_id].x2 = octree.full_bins[bin_id].x_mean - octree.direction_vec .* octree.full_bins[bin_id].x_std_sq
        elseif (octree.bins[bin_id].np == 2)
            # get the particle indices we saved and just write data based on them
            i = octree.full_bins[bin_id].particle_index1
            octree.full_bins[bin_id].w1 = particles[i].w
            octree.full_bins[bin_id].v1 = particles[i].v
            octree.full_bins[bin_id].x1 = particles[i].x

            i = octree.full_bins[bin_id].particle_index2
            octree.full_bins[bin_id].w2 = particles[i].w
            octree.full_bins[bin_id].v2 = particles[i].v
            octree.full_bins[bin_id].x2 = particles[i].x
        elseif (octree.bins[bin_id].np == 1)
            # get the particle indices we saved and just write data based on them
            i = octree.full_bins[bin_id].particle_index1
            octree.full_bins[bin_id].w1 = particles[i].w
            octree.full_bins[bin_id].v1 = particles[i].v
            octree.full_bins[bin_id].x1 = particles[i].x
        end
    end

    curr_particle_index = 0
    for bin_id in 1:octree.Nbins
        if (octree.bins[bin_id].np >= 2)
            i = map_cont_index(pia.indexer[cell,species], curr_particle_index)
            curr_particle_index += 1
            particles[i].w = octree.full_bins[bin_id].w1
            particles[i].v = octree.full_bins[bin_id].v1
            particles[i].x = octree.full_bins[bin_id].x1

            i = map_cont_index(pia.indexer[cell,species], curr_particle_index)
            curr_particle_index += 1
            particles[i].w = octree.full_bins[bin_id].w2
            particles[i].v = octree.full_bins[bin_id].v2
            particles[i].x = octree.full_bins[bin_id].x2
        elseif (octree.bins[bin_id].np == 1)
            i = map_cont_index(pia.indexer[cell,species], curr_particle_index)
            curr_particle_index += 1
            particles[i].w = octree.full_bins[bin_id].w1
            particles[i].v = octree.full_bins[bin_id].v1
            particles[i].x = octree.full_bins[bin_id].x1
        end
    end

    old_count = pia.indexer[cell,species].n_local
    n_particles_to_delete = old_count - curr_particle_index

    # if we delete from particles in last cell AND we delete less particles than were in group 2
    # then continuity is not broken
    # !(A && B) == !A || !B
    if !(cell == size(pia.indexer)[1]) || (n_particles_to_delete > pia.indexer[cell,species].n_group2)
        pia.contiguous[species] = false
    end

    for _ in 1:n_particles_to_delete
        delete_particle_end!(particles, pia, cell, species)
    end
end


"""
    compute_new_particles!(rng, octree::OctreeN2Merge, particles, pia, cell, species, grid)

Compute post-merge particles particles based on octree bin properties; placing out-of-domain particles back into the domain.

# Positional arguments:
* `rng`: the random number generator instance
* `octree`: the `OctreeN2Merge` instance
* `particles`: the `ParticleVector` instance of the particles to be merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the cell index
* `species`: the species index
* `grid`: the `Grid1DUniform` grid
"""
function compute_new_particles!(rng, octree::OctreeN2Merge, particles, pia, cell, species, grid::Grid1DUniform)
    # given computed Octree, create new particles instead of the old ones
    
    for bin_id in 1:octree.Nbins
        if (octree.bins[bin_id].np > 2)
            octree.full_bins[bin_id].w1 = 0.5 * octree.bins[bin_id].w
            octree.full_bins[bin_id].w2 = octree.full_bins[bin_id].w1

            octree.full_bins[bin_id].v_std_sq = sqrt.(octree.full_bins[bin_id].v_std_sq)
            octree.full_bins[bin_id].x_std_sq = sqrt.(octree.full_bins[bin_id].x_std_sq)
            
            octree.direction_vec = @SVector rand(rng, direction_signs, 3)
            octree.full_bins[bin_id].v1 = octree.full_bins[bin_id].v_mean + octree.direction_vec .* octree.full_bins[bin_id].v_std_sq
            octree.full_bins[bin_id].v2 = octree.full_bins[bin_id].v_mean - octree.direction_vec .* octree.full_bins[bin_id].v_std_sq

            octree.direction_vec = @SVector rand(rng, direction_signs, 3)
            octree.full_bins[bin_id].x1 = octree.full_bins[bin_id].x_mean + octree.direction_vec .* octree.full_bins[bin_id].x_std_sq
            octree.full_bins[bin_id].x2 = octree.full_bins[bin_id].x_mean - octree.direction_vec .* octree.full_bins[bin_id].x_std_sq
        elseif (octree.bins[bin_id].np == 2)
            # get the particle indices we saved and just write data based on them
            i = octree.full_bins[bin_id].particle_index1
            octree.full_bins[bin_id].w1 = particles[i].w
            octree.full_bins[bin_id].v1 = particles[i].v
            octree.full_bins[bin_id].x1 = particles[i].x

            i = octree.full_bins[bin_id].particle_index2
            octree.full_bins[bin_id].w2 = particles[i].w
            octree.full_bins[bin_id].v2 = particles[i].v
            octree.full_bins[bin_id].x2 = particles[i].x
        elseif (octree.bins[bin_id].np == 1)
            # get the particle indices we saved and just write data based on them
            i = octree.full_bins[bin_id].particle_index1
            octree.full_bins[bin_id].w1 = particles[i].w
            octree.full_bins[bin_id].v1 = particles[i].v
            octree.full_bins[bin_id].x1 = particles[i].x
        end
    end

    curr_particle_index = 0
    for bin_id in 1:octree.Nbins
        if (octree.bins[bin_id].np >= 2)
            i = map_cont_index(pia.indexer[cell,species], curr_particle_index)
            curr_particle_index += 1
            particles[i].w = octree.full_bins[bin_id].w1
            particles[i].v = octree.full_bins[bin_id].v1

            if (octree.full_bins[bin_id].x1[1] < grid.min_x)
                particles[i].x = [grid.min_x, octree.full_bins[bin_id].x1[2], octree.full_bins[bin_id].x1[3]]
            elseif (octree.full_bins[bin_id].x1[1] > grid.max_x)
                particles[i].x = [grid.max_x, octree.full_bins[bin_id].x1[2], octree.full_bins[bin_id].x1[3]]
            else
                particles[i].x = octree.full_bins[bin_id].x1
            end

            i = map_cont_index(pia.indexer[cell,species], curr_particle_index)
            curr_particle_index += 1
            particles[i].w = octree.full_bins[bin_id].w2
            particles[i].v = octree.full_bins[bin_id].v2

            if (octree.full_bins[bin_id].x2[1] < grid.min_x)
                particles[i].x = [grid.min_x, octree.full_bins[bin_id].x2[2], octree.full_bins[bin_id].x2[3]]
            elseif (octree.full_bins[bin_id].x2[1] > grid.max_x)
                particles[i].x = [grid.max_x, octree.full_bins[bin_id].x2[2], octree.full_bins[bin_id].x2[3]]
            else
                particles[i].x = octree.full_bins[bin_id].x2
            end
        elseif (octree.bins[bin_id].np == 1)
            i = map_cont_index(pia.indexer[cell,species], curr_particle_index)
            curr_particle_index += 1
            particles[i].w = octree.full_bins[bin_id].w1
            particles[i].v = octree.full_bins[bin_id].v1
            particles[i].x = octree.full_bins[bin_id].x1
        end
    end

    old_count = pia.indexer[cell,species].n_local
    n_particles_to_delete = old_count - curr_particle_index

    # if we delete from particles in last cell AND we delete less particles than were in group 2
    # then continuity is not broken
    # !(A && B) == !A || !B
    if !(cell == size(pia.indexer)[1]) || (n_particles_to_delete > pia.indexer[cell,species].n_group2)
        pia.contiguous[species] = false
    end

    for _ in 1:n_particles_to_delete
        delete_particle_end!(particles, pia, cell, species)
    end
end

"""
    init_octree!(octree, particles, pia, cell, species)

Initialize the top bin in an octree by copying particle indices and setting bin bounds.

# Positional arguments
* `octree`: the `OctreeN2Merge` instance
* `particles`: the `ParticleVector` instance containing the particles to be merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the grid cell in which particles are being merged
* `species`: the index of the species being merged
"""
function init_octree!(octree, particles, pia, cell, species)
    octree.Nbins = 1
    octree.particle_indexes_sorted[1:pia.indexer[cell,species].n_group1] = pia.indexer[cell,species].start1:pia.indexer[cell,species].end1

    if (pia.indexer[cell,species].n_group2 > 0)
        octree.particle_indexes_sorted[pia.indexer[cell,species].n_group1+1:pia.indexer[cell,species].n_local] = pia.indexer[cell,species].start2:pia.indexer[cell,species].end2
    end

    octree.n_particles = pia.indexer[cell,species].n_local
    octree.bins[1].depth = 0
    octree.bin_start[1] = 1
    octree.bin_end[1] = octree.n_particles

    octree.bins[1].np = octree.n_particles
    octree.bins[1].w = 1e50  # just do a large estimate, we will be splitting this bin anyway

    octree.total_post_merge_np = get_bin_post_merge_np(octree, 1)
    if (octree.bins[1].np > 2) && (octree.max_depth > 0)
        octree.bins[1].can_be_refined = true
    else
        octree.bins[1].can_be_refined = false
    end

    if (octree.init_bin_bounds == OctreeInitBinC)
        octree.bins[1].v_min = SVector{3, Float64}(-299_792_458.0, -299_792_458.0, -299_792_458.0)  # speed of light
        octree.bins[1].v_max = SVector{3, Float64}(299_792_458.0, 299_792_458.0, 299_792_458.0)
    else
        bin_bounds_recompute!(octree, 1, 1, octree.n_particles, particles)
        if (octree.init_bin_bounds == OctreeInitBinMinMaxVelSym)
            maxvx = max(abs(octree.bins[1].v_min[1]), abs(octree.bins[1].v_max[1]))
            maxvy = max(abs(octree.bins[1].v_min[2]), abs(octree.bins[1].v_max[2]))
            maxvz = max(abs(octree.bins[1].v_min[3]), abs(octree.bins[1].v_max[3]))

            octree.bins[1].v_min = SVector{3, Float64}(-maxvx, -maxvy, -maxvz)
            octree.bins[1].v_max = SVector{3, Float64}(maxvx, maxvy, maxvz)
        end
    end
end

"""
    merge_octree_N2_based!(rng, octree, particles, pia, cell, species, target_np)

Perform N:2 merging without checking whether particle positions end up outside of the simulation domain.

# Positional arguments
* `rng`: the random number generator instance
* `octree`: the `OctreeN2Merge` instance
* `particles`: the `ParticleVector` instance containing the particles to be merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the grid cell in which particles are being merged
* `species`: the index of the species being merged
* `target_np`: the target post-merge number of particles; the post-merge number of particles will not exceed this value
    but may be not exactly equal to it
"""
function merge_octree_N2_based!(rng, octree, particles, pia, cell, species, target_np)
    clear_octree!(octree)
    resize_octree_buffers!(octree, pia.indexer[cell,species].n_local)
    init_octree!(octree, particles, pia, cell, species)

    while true
        refine_id = -1
        max_w = -1.0
        for bin_id in 1:octree.Nbins
            # find bin with largest number density, needs to be refine-able
            if ((octree.bins[bin_id].w > max_w) && octree.bins[bin_id].can_be_refined)
               max_w = octree.bins[bin_id].w
               refine_id = bin_id
           end
        end

        if refine_id == -1
            # found no bin to refine, i.e. ran out of particles
            break
        elseif (octree.total_post_merge_np + 14 > target_np)
            # refining a bin can produce up to 16 particles, so we don't do it if threshold exceeded
            # but if we have a bin it has potentially 2 particles already, so refinement increases count
            # only by 14
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
        octree.bins[1].w = 0.0
        for pi in octree.particle_indexes_sorted[1:octree.n_particles]
            octree.bins[1].w += particles[pi].w
        end
    end
    
    for bin_id in 1:octree.Nbins
        compute_bin_props!(octree, bin_id, particles)
    end
    compute_new_particles!(rng, octree, particles, pia, cell, species)
end


"""
    merge_octree_N2_based!(rng, octree, particles, pia, cell, species, target_np, grid)

Perform N:2 merging, checking whether particle positions end up outside of the simulation domain, and pushing them back into the domain
if needed.

# Positional arguments
* `rng`: the random number generator instance
* `octree`: the `OctreeN2Merge` instance
* `particles`: the `ParticleVector` instance containing the particles to be merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the grid cell in which particles are being merged
* `species`: the index of the species being merged
* `target_np`: the target post-merge number of particles; the post-merge number of particles will not exceed this value
    but may be not exactly equal to it
* `grid`: the `Grid1DUniform` grid
"""
function merge_octree_N2_based!(rng, octree, particles, pia, cell, species, target_np, grid::Grid1DUniform)
    clear_octree!(octree)
    resize_octree_buffers!(octree, pia.indexer[cell,species].n_local)
    init_octree!(octree, particles, pia, cell, species)

    while true
        refine_id = -1
        max_w = -1.0
        for bin_id in 1:octree.Nbins
            # find bin with largest number density, needs to be refine-able
            if ((octree.bins[bin_id].w > max_w) && octree.bins[bin_id].can_be_refined)
               max_w = octree.bins[bin_id].w
               refine_id = bin_id
           end
        end

        if refine_id == -1
            # found no bin to refine, i.e. ran out of particles
            break
        elseif (octree.total_post_merge_np + 14 > target_np)
            # refining a bin can produce up to 16 particles, so we don't do it if threshold exceeded
            # but if we have a bin it has potentially 2 particles already, so refinement increases count
            # only by 14
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
        octree.bins[1].w = 0.0
        for pi in octree.particle_indexes_sorted[1:octree.n_particles]
            octree.bins[1].w += particles[pi].w
        end
    end
    
    for bin_id in 1:octree.Nbins
        compute_bin_props!(octree, bin_id, particles)
    end
    compute_new_particles!(rng, octree, particles, pia, cell, species, grid)
end