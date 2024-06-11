using StaticArrays

@enum OctreeBinSplit OctreeBinMidSplit=1 OctreeBinMeanSplit=2 OctreeBinMedianSplit=3

# initial bin bounds:
# 1) min/max velocities of particles
# 2) min/max but symmetrized (-max(abs(min), abs(max)), -max(abs(min), abs(max)))
# 3) speed of light
@enum OctreeInitBin  OctreeInitBinMinMaxVel=1 OctreeInitBinMinMaxVelSym=2 OctreeInitBinC=3

# how bounds of split bins are computed: using splitting velocity and bounds of bounding bin
# or using min/max velocities of particles in sub-bin
@enum OctreeBinBounds OctreeBinBoundsInherit=1 OctreeBinBoundsRecompute=2

mutable struct OctreeCell
    # this holds only the data needed for refinement
    np::Int64
    w::Float64
    
    v_min::SVector{3,Float64}
    v_max::SVector{3,Float64}

    depth::Int64
end

mutable struct OctreeFullCell
    # this holds the additional data needed for actual merging
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
mutable struct OctreeN2
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
end

function fill_bins(Nbins)
    return [OctreeCell(0, 0.0, SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0), 0) for i in 1:Nbins]
end

function fill_full_bins(Nbins)
    return [OctreeFullCell(SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0),
                           SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0),
                           0, 0, 0.0, 0.0, 
                           SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0),
                           SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0)) for i in 1:Nbins]
end

function create_merging_octree(split::OctreeBinSplit; init_bin_bounds=OctreeInitBinMinMaxVel, bin_bounds_compute=OctreeBinBoundsInherit,
                               max_Nbins=4096, max_depth=10)
    return OctreeN2(max_Nbins, 0, fill_bins(max_Nbins), fill_full_bins(max_Nbins), 0,
                    zeros(max_Nbins), zeros(max_Nbins),  # bin_start, bin_end
                    zeros(8192), zeros(8192), zeros(8192),
                    MVector{8, Int64}(0, 0, 0, 0, 0, 0, 0, 0),
                    MVector{8, Int64}(0, 0, 0, 0, 0, 0, 0, 0),
                    MVector{8, Int64}(0, 0, 0, 0, 0, 0, 0, 0),
                    MVector{8, Float64}(0, 0, 0, 0, 0, 0, 0, 0),
                    bin_bounds_compute, split,
                    SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0), SVector{3,Float64}(0.0, 0.0, 0.0),
                    SVector{3,Float64}(0.0, 0.0, 0.0),
                    init_bin_bounds, max_depth)
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

function bin_bounds_inherit!(octree, bin_id, v_min_parent, v_max_parent, v_middle, octant)
    # TODO: TEST
    # inherit bin bounds based on parent
    # octants order:
    # - - -
    # + - -
    # - + -
    # + + -
    # - - +
    # + - +
    # - + +
    # + + +
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

function bin_bounds_recompute!(octree, bin_id, bs, be, particles)
    # compute bin bounds based on the particles in the bin
    minvx = 9_299_792_458.0  # speed of light + 9e9
    minvy = 9_299_792_458.0
    minvz = 9_299_792_458.0

    maxvx = -9_299_792_458.0
    maxvy = -9_299_792_458.0
    maxvz = -9_299_792_458.0

    for pi in octree.particle_indexes_sorted[bs:be]
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

function compute_v_mean!(octree, bs, be, particles)
    n_tot = 0.0
    for pi in octree.particle_indexes_sorted[bs:be]
        n_tot += particles[pi].w
        octree.vel_middle = octree.vel_middle + particles[pi].w * particles[pi].v
    end
    octree.vel_middle = octree.vel_middle / n_tot
end

function bin_bounds_recompute_and_v_mean!(octree, bin_id, bs, be, particles)
    # do everything in 1 pass over the particles
end

function split_bin!(octree, bin_id, particles)
    octree.particle_in_bin_counter .= 0 # reset counter
    octree.ndens_counter .= 0.0
    octree.nonempty_bins .= 0
    octree.nonempty_counter .= 0

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
        # TODO
        octree.vel_middle = SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    for (i, pi) in enumerate(octree.particle_indexes_sorted[bs:be])
        oct = compute_octant(particles[pi].v, octree.vel_middle)
        octree.particle_in_bin_counter[oct] += 1
        octree.particle_octants[i] = oct
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

    octree.Nbins += n_nonempty_bins - 1

    # shift bins around to accomodate the new non-empty bins we have
    octree.bin_start[bin_id+n_nonempty_bins:octree.Nbins+n_nonempty_bins-1] .= octree.bin_start[bin_id+1:octree.Nbins]
    octree.bin_end[bin_id+n_nonempty_bins:octree.Nbins+n_nonempty_bins-1] .= octree.bin_end[bin_id+1:octree.Nbins]

    # TODO: loop instead of deepcopy? would that be faster?
    octree.bins[bin_id+n_nonempty_bins:octree.Nbins+n_nonempty_bins-1] .= deepcopy(octree.bins[bin_id+1:octree.Nbins])

    # first bin - we change nothing for the start
    octree.bin_end[bin_id] = octree.bin_start[bin_id] + octree.nonempty_counter[1] - 1
    
    for i in 1:n_nonempty_bins-1
        octree.bin_start[bin_id+i] = octree.bin_end[bin_id+i-1] + 1
        octree.bin_end[bin_id+i] = octree.bin_start[bin_id+i] + octree.nonempty_counter[i+1] - 1
    end

    if (octree.bin_bounds_compute == OctreeBinBoundsInherit)
        octree.v_min_parent = octree.bins[bin_id].v_min
        octree.v_max_parent = octree.bins[bin_id].v_max

        # iterate over non-empty bins and inherit parent bin bounds + split around middle velocity
        for i in 1:n_nonempty_bins
            bin_bounds_inherit!(octree, bin_id + i - 1,
                                octree.v_min_parent, octree.v_max_parent,
                                octree.vel_middle, octree.nonempty_bins[i])
            octree.bins[bin_id + i - 1].np = octree.nonempty_counter[i]
            octree.bins[bin_id + i - 1].w = octree.ndens_counter[octree.nonempty_bins[i]]
        end
    else
        # still need to fill out info on number of particles and total weight
        # will recompute bin bounds if we do next round of refinement
        for i in 1:n_nonempty_bins
            octree.bins[bin_id + i - 1].np = octree.nonempty_counter[i]
            octree.bins[bin_id + i - 1].w = octree.ndens_counter[octree.nonempty_bins[i]]
        end
    end

    # for (i, pi) in enumerate(octree.particle_indexes_sorted[be:-1:bs])
    for (i, pi) in enumerate(octree.particle_indexes_sorted[bs:be])
        j = octree.particle_octants[i]
        octree.particles_sort_output[octree.particle_in_bin_counter[j]] = pi
        octree.particle_in_bin_counter[j] -= 1
    end

    # write sorted indices
    octree.particle_indexes_sorted[bs:be] = octree.particles_sort_output[1:be-bs+1]
end

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

    for pi in octree.particle_indexes_sorted[bs:be]
        octree.full_bins[bin_id].v_mean = octree.full_bins[bin_id].v_mean + particles[pi].w * particles[pi].v
        octree.full_bins[bin_id].x_mean = octree.full_bins[bin_id].x_mean + particles[pi].w * particles[pi].x
    end
    octree.full_bins[bin_id].v_mean = octree.full_bins[bin_id].v_mean / octree.bins[bin_id].w
    octree.full_bins[bin_id].x_mean = octree.full_bins[bin_id].x_mean / octree.bins[bin_id].w

    for pi in octree.particle_indexes_sorted[bs:be]
        octree.full_bins[bin_id].v_std_sq = octree.full_bins[bin_id].v_std_sq +
                                            particles[pi].w * (particles[pi].v - octree.full_bins[bin_id].v_mean).^2
        octree.full_bins[bin_id].x_std_sq = octree.full_bins[bin_id].x_std_sq +
                                            particles[pi].w * (particles[pi].x - octree.full_bins[bin_id].x_mean).^2
    end
    octree.full_bins[bin_id].v_std_sq = octree.full_bins[bin_id].v_std_sq / octree.bins[bin_id].w
    octree.full_bins[bin_id].x_std_sq = octree.full_bins[bin_id].x_std_sq / octree.bins[bin_id].w
end

function get_bin_post_merge_np(octree, bin_id)
    # how many particles will we get after merging in the bin:
    # 2 if np >= 2
    # np otherwise (0 or 1)
    return octree.bins[bin_id].np >= 2 ? 2 : octree.bins[bin_id].np
end

function compute_new_particles!(cell, species, octree::OctreeN2, particles, particle_indexer_array)
    # given computed Octree, create new particles instead of the old ones
    
    for bin_id in 1:octree.Nbins
        if (octree.bins[bin_id].np > 2)
            octree.full_bins[bin_id].w1 = 0.5 * octree.bins[bin_id].w
            octree.full_bins[bin_id].w2 = octree.full_bins[bin_id].w1

            octree.full_bins[bin_id].v_std_sq = sqrt.(octree.full_bins[bin_id].v_std_sq)
            octree.full_bins[bin_id].x_std_sq = sqrt.(octree.full_bins[bin_id].x_std_sq)
            
            octree.direction_vec = @SVector rand(direction_signs, 3)
            octree.full_bins[bin_id].v1 = octree.full_bins[bin_id].v_mean + octree.direction_vec .* octree.full_bins[bin_id].v_std_sq
            octree.full_bins[bin_id].v2 = octree.full_bins[bin_id].v_mean - octree.direction_vec .* octree.full_bins[bin_id].v_std_sq

            octree.direction_vec = @SVector rand(direction_signs, 3)
            octree.full_bins[bin_id].x1 = octree.full_bins[bin_id].x_mean + octree.direction_vec .* octree.full_bins[bin_id].x_std_sq
            octree.full_bins[bin_id].x2 = octree.full_bins[bin_id].x_mean - octree.direction_vec .* octree.full_bins[bin_id].x_std_sq
        elseif (octree.bins[bin_id].np == 2)
            # get the particle indices we saved and just write data based on them
            i = merging_grid.cells[index].particle_index1
            octree.full_bins[bin_id].w1 = particles[species][i].w
            octree.full_bins[bin_id].v1 = particles[species][i].v
            octree.full_bins[bin_id].x1 = particles[species][i].x

            i = merging_grid.cells[index].particle_index2
            octree.full_bins[bin_id].w2 = particles[species][i].w
            octree.full_bins[bin_id].v2 = particles[species][i].v
            octree.full_bins[bin_id].x2 = particles[species][i].x
        elseif (octree.bins[bin_id].np == 1)
            # get the particle indices we saved and just write data based on them
            i = merging_grid.cells[index].particle_index1
            octree.full_bins[bin_id].w1 = particles[species][i].w
            octree.full_bins[bin_id].v1 = particles[species][i].v
            octree.full_bins[bin_id].x1 = particles[species][i].x
        end
    end

    curr_particle_index = 0
    for bin_id in 1:octree.Nbins
        if (octree.bins[bin_id].np >= 2)
            i = map_cont_index(particle_indexer_array.indexer[cell,species], curr_particle_index)
            curr_particle_index += 1
            particles[species][i].w = octree.full_bins[bin_id].w1
            particles[species][i].v = octree.full_bins[bin_id].v1
            particles[species][i].x = octree.full_bins[bin_id].x1

            i = map_cont_index(particle_indexer_array.indexer[cell,species], curr_particle_index)
            curr_particle_index += 1
            particles[species][i].w = octree.full_bins[bin_id].w2
            particles[species][i].v = octree.full_bins[bin_id].v2
            particles[species][i].x = octree.full_bins[bin_id].x2
        elseif (octree.bins[bin_id].np == 1)
            i = map_cont_index(particle_indexer_array.indexer[cell,species], curr_particle_index)
            curr_particle_index += 1
            particles[species][i].w = octree.full_bins[bin_id].w1
            particles[species][i].v = octree.full_bins[bin_id].v1
            particles[species][i].x = octree.full_bins[bin_id].x1
        end
    end

    update_particle_indexer_new_lower_count(species, cell, particle_indexer_array, curr_particle_index)
end

# initialize first bin with all the particles - set number of bins to 1, copy over particle indices
function init_octree!(cell, species, octree, particles, pia)
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

function merge_octree_N2_based!(cell, species, octree, particles, particle_indexer_array, target_np)
    clear_octree!(octree)
    resize_octree_buffers!(octree, particle_indexer_array.indexer[cell,species].n_local)
    init_octree!(cell, species, octree, particles[species], particle_indexer_array)

    while true
        total_np = 0
        refine_id = -1
        max_w = -1
        for bin_id in 1:octree.Nbins
            total_np += get_bin_post_merge_np(octree, bin_id)
            if ((octree.bins[bin_id].w > max_w) && (octree.bins[bin_id].np > 2))
                max_w = octree.bins[bin_id].w
                refine_id = bin_id
            end
        end

        if refine_id == -1
            # found no bin to refine, i.e. ran out of particles
            break
        elseif (total_np + 14 > target_np)
            # refining a bin can produce up to 16 particles, so we don't do it if threshold exceeded
            # but if we have a bin it has potentially 2 particles already, so refinement increases count
            # only by 14
            break
        else
            split_bin!(octree, refine_id, particles[species])
        end
    end

    if octree.Nbins == 1
        octree.bins[1].w = 0.0
        for pi in octree.particle_indexes_sorted[1:octree.n_particles]
            octree.bins[1].w += particles[species][pi].w
        end
    end
    
    for bin_id in 1:octree.Nbins
        compute_bin_props!(octree, bin_id, particles[species])
    end
    compute_new_particles!(cell, species, octree, particles, particle_indexer_array)
end