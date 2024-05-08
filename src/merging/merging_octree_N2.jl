using StaticArrays

@enum OctreeBinSplit OctreeBinMidSplit=1 OctreeBinMeanSplit=2 OctreeBinMedianSplit=3

# initial bin bounds:
# 1) min/max velocities of particles
# 2) min/max but symmetrized (-max(abs(min), abs(max)), -max(abs(min), abs(max)))
# 3) speed of light
@enum OctreeInitBin  OctreeInitBinMinMaxVel=1 OctreeInitBinMinMaxVelSym=2 OctreeInitBinC=3

mutable struct OctreeCell
    np::Int64
    w::Float64
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

    v_min::SVector{3,Float64}
    v_max::SVector{3,Float64}

    depth::Int64
end

# struct for N:2 merge
mutable struct OctreeN2
    max_Nbins::Int32
    Nbins::Int32  # actual bins computed
    bins::Vector{OctreeCell}
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

    split::OctreeBinSplit
    vel_middle::SVector{3,Float64}  # defines how we split octant

    init_bin_bounds::OctreeInitBin

    max_depth::Int32
end

function fill_bins(Nbins)
    return [OctreeCell(0, 0.0, SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0),
            0, 0, 0.0, 0.0, SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0),
            SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0), 0) for i in 1:Nbins]
end

function create_merging_octree(split::OctreeBinSplit; init_bin_bounds=OctreeInitBinMinMaxVel, max_Nbins=4096, max_depth=10)
    return OctreeN2(max_Nbins, 0, fill_bins(max_Nbins), 0,
                    zeros(max_Nbins), zeros(max_Nbins),  # bin_start, bin_end
                    zeros(8192), zeros(8192), zeros(8192),
                    MVector{8, Int64}(0, 0, 0, 0, 0, 0, 0, 0),
                    split, SVector{3,Float64}(0.0, 0.0, 0.0), init_bin_bounds, max_depth)
end

function clear_octree!(octree)
    octree.Nbins = 0
end

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

function compute_octant(particle_v, v_middle)

    # octants order:
    # - - -
    # + - -
    # - + -
    # + + -
    # - - +
    # + - +
    # - + +
    # + + +

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

function median_vel()
    # https://rcoh.me/posts/linear-time-median-finding/
    # https://en.wikipedia.org/wiki/Weighted_median
    nothing
end

function split_bin!(octree, bin_id, particles)
    octree.particle_in_bin_counter .= 0 # reset counter

    n_nonempty_bins = 0
    bs = octree.bin_start[bin_id]
    be = octree.bin_end[bin_id]

    if (octree.split == OctreeBinMidSplit)
        octree.vel_middle = 0.5 * (octree.bins[bin_id].v_min + octree.bins[bin_id].v_max)
    elseif (octree.split == OctreeBinMeanSplit)
        octree.vel_middle = octree.bins[bin_id].v_mean
    elseif (octree.split == OctreeBinMedianSplit)
        octree.vel_middle = SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    for (i, pi) in enumerate(octree.particle_indexes_sorted[bs:be])
        oct = compute_octant(particles[pi].v, octree.vel_middle)
        octree.particle_in_bin_counter[oct] += 1
        octree.particle_octants[i] = oct
    end

    if (octree.particle_in_bin_counter[1] > 0)
        n_nonempty_bins += 1
    end

    for i in 2:8
        if (octree.particle_in_bin_counter[i] > 0)
            n_nonempty_bins += 1
        end

        octree.particle_in_bin_counter[i] += octree.particle_in_bin_counter[i-1]
    end

    octree.Nbins += n_nonempty_bins - 1

    # shift bins around to accomodate the new non-empty bins we have
    octree.bin_start[bin_id+n_nonempty_bins:octree.Nbins+n_nonempty_bins-1] .= octree.bin_start[bin_id+1:octree.Nbins]
    octree.bin_end[bin_id+n_nonempty_bins:octree.Nbins+n_nonempty_bins-1] .= octree.bin_end[bin_id+1:octree.Nbins]
    octree.bins[bin_id+n_nonempty_bins:octree.Nbins+n_nonempty_bins-1] .= octree.bins[bin_id+1:octree.Nbins]

    # TODO: compute props

    # for (i, pi) in enumerate(octree.particle_indexes_sorted[be:-1:bs])
    for (i, pi) in enumerate(octree.particle_indexes_sorted[bs:be])
        j = octree.particle_octants[i]
        octree.particles_sort_output[octree.particle_in_bin_counter[j]] = pi
        octree.particle_in_bin_counter[j] -= 1
    end

    # write sorted indices
    octree.particle_indexes_sorted[bs:be] = octree.particles_sort_output[1:be-bs+1]
end

function compute_new_particles!()
    # given computed Octree, create new particles instead of the old ones
    nothing
end

# initialize first bin with all the particles - set number of bins to 1, copy over particle indices
function init_octree!(cell, species, octree, pia)
    octree.Nbins = 1
    octree.particle_indexes_sorted[1:pia.indexer[cell,species].n_group1] = pia.indexer[cell,species].start1:pia.indexer[cell,species].end1

    if (pia.indexer[cell,species].n_group2 > 0)
        octree.particle_indexes_sorted[pia.indexer[cell,species].n_group1+1:pia.indexer[cell,species].n_local] = pia.indexer[cell,species].start2:pia.indexer[cell,species].end2
    end

    octree.n_particles = pia.indexer[cell,species].n_local
    octree.bins[1].depth = 0
    octree.bin_start[1] = 1
    octree.bin_end[1] = octree.n_particles

    if (octree.init_bin_bounds == OctreeInitBinC)
        octree.bins[1].v_min = SVector{3, Float64}(-299_792_458.0, -299_792_458.0, -299_792_458.0)
        octree.bins[1].v_max = SVector{3, Float64}(299_792_458.0, 299_792_458.0, 299_792_458.0)
    else
        # TODO: compute min/max velocities
        if (octree.init_bin_bounds == OctreeInitBinMinMaxVelSym)
            nothing # TODO: symmetrize
        end
    end
end

function merge_octree_N2_based!(cell, species, octree, particles, particle_indexer_array, target_np)
    clear_octree!(octree)
    resize_octree_buffers!(octree, particle_indexer_array.indexer[cell,species].n_local)
    init_octree!(cell, species, octree, particle_indexer_array)
    # loop until reached target
    # inner loop - choose bin to refine, refine, exit inner loop

    compute_new_particles!()
end