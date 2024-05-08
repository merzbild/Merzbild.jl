using StaticArrays

@enum OctreeBinSplit OctreeBinMidSplit=1 OctreeBinMeanSplit=2 OctreeBinMedianSplit=3

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
    max_Nbins::Int64
    Nbins::Int64  # actual bins computed
    bins::OctreeCell
    n_particles::Int64  # particles being sorted

    # particles in bins[i] have indices in particle_index_buffer[bin_start[i]:bin_end[i]]
    bin_start::Vector{Int64}
    bin_end::Vector{Int64}

    # this stores sorted particle indices
    particle_indexes_sorted::Vector{Int64}

    # to store particle octants during sort, default size is 8192
    # in case particles don't fit it's increase to length(particles) + DELTA_PARTICLES
    particle_octants::Vector{Int64}

    # used to store sorted particle indices
    # in case particles don't fit it's increase to length(particles) + DELTA_PARTICLES
    particles_sort_output::Vector{Int64}


    # count how many particles in each bin, used for radix sort
    particle_in_bin_counter::MVector{8, Int64}

    split::OctreeBinSplit
    max_depth::Int64

    vel_middle::SVector{3,Float64}  # defines how we split octant
end

function create_merging_octree()
    nothing
end

function clear_octree!()
    # clear out data from previous merge
    nothing
end

function resize_octree_buffers!(octree, n_particles)
    nothing
end

function split_bin!(octree, bin_id, particles)
    octree.particle_in_bin_counter .= 0 # reset counter

    n_nonempty_bins = 0
    bs = bin_start[bin_id]
    be = bin_end[bin_id]

    if (octree.split == OctreeBinMidSplit)
        octree.vel_middle = 0.5 * (octree.bins[bin_id].v_min + octree.bins[bin_id].v_max)
    elseif (octree.split == OctreeBinMeanSplit)
        octree.vel_middle = octree.bins[bin_id].v_mean
    elseif (octree.split == OctreeBinMedianSplit)
        octree.vel_middle = SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    for (i, pi) in enumerate(octree.particle_indexes_sorted[bin_start[bin_id]:bin_end[bin_id]])
        oct = compute_octant(particles[pi], octree.bins[bin_id], vel_middle)
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

    # shift bins around to accomodate the new non-empty bins we have
    octree.bin_start[bin_id+n_nonempty_bins:octree.Nbins+n_nonempty_bins-1] .= octree.bin_start[bin_id+1:octree.Nbins]
    octree.bin_end[bin_id+n_nonempty_bins:octree.Nbins+n_nonempty_bins-1] .= octree.bin_end[bin_id+1:octree.Nbins]
    octree.bins[bin_id+n_nonempty_bins:octree.Nbins+n_nonempty_bins-1] .= octree.bins[bin_id+1:octree.Nbins]

    # TODO: compute props

    for (i, pi) in enumerate(particle_indexes_sorted[be:-1:bs])
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

function merge_octree_N2_based!(target_np)
    clear_octree!()
    
    # loop until reached target
    # inner loop - choose bin to refine, refine, exit inner loop

    compute_new_particles!()
end