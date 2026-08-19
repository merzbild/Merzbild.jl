"""
    ChunkExchanger

Struct used to organize exchange of particles between independent ParticleVector instances for
chunked multi-threaded simulations. It is assumed that the cell indices within each chunk are contiguous,
i.e. `chunk[i+1] = chunk[i]+1`.

The indexing for the exchanged groups is stored as arrays of shape `(n_chunks, n_cells)`,
and `[chunk_id, cell]` describes the particles that belong to cell `cell` and came from
particle chunk `chunk_id`. Group 1 holds the particles that arrived via swapping, group 2 the
particles that arrived via pushing, i.e. those added to the end of the particle array.
The end of a group is not stored, as it is given by `start + n_group - 1`
(which is `-1` for an empty group, matching the convention used by `ParticleIndexer`).

The table is self-clearing: every entry written during an exchange is read exactly once,
in the same timestep, by the chunk owning the cell, and
[`sort_particles_after_exchange!`](@ref) zeroes it as it reads it.
It is therefore already empty at the start of each timestep and does not need to be
[`reset!`](@ref) between timesteps.

The `occ_lo`/`occ_hi` fields hold, per chunk, the first and last cell in which the chunk
holds any particles. They let [`exchange_particles!`](@ref) reject a pair of chunks that
cannot have anything to exchange without scanning the cells of either chunk. They are
updated by [`update_occupancy_bounds!`](@ref) and start out as the whole grid, i.e. a
simulation that never calls it is still correct, just slower.

# Fields
* `n_chunks`: number of chunks used in the simulation
* `n_cells`: number of grid cells in the simulation
* `start1`: index of the first particle of group 1
* `n_group1`: number of particles in group 1
* `start2`: index of the first particle of group 2
* `n_group2`: number of particles in group 2
* `occ_lo`: first cell in which a chunk holds particles
* `occ_hi`: last cell in which a chunk holds particles
"""
mutable struct ChunkExchanger
    n_chunks::Int64
    n_cells::Int64
    start1::Array{Int64, 2}  # n_chunks x n_cells
    n_group1::Array{Int64, 2}
    start2::Array{Int64, 2}
    n_group2::Array{Int64, 2}
    occ_lo::Vector{Int64}  # n_chunks
    occ_hi::Vector{Int64}
end

"""
    ChunkExchanger(chunks, n_cells)

Create a `ChunkExchanger` for `length(chunks)` chunks and `n_cells` cell.

# Positional arguments
* `chunks`: list of cell chunks
* `n_cells`: total number of cells in the simulation
"""
function ChunkExchanger(chunks, n_cells)
    n_chunks = length(chunks)
    return ChunkExchanger(n_chunks, n_cells,
                          zeros(Int64, n_chunks, n_cells), zeros(Int64, n_chunks, n_cells),
                          zeros(Int64, n_chunks, n_cells), zeros(Int64, n_chunks, n_cells),
                          ones(Int64, n_chunks), fill(n_cells, n_chunks))
end

"""
    update_occupancy_bounds!(chunk_exchanger, gridsort, pia, chunk_id, species)

Update the first and last cell in which chunk `chunk_id` holds particles of the given `species`,
so that [`exchange_particles!`](@ref) can reject pairs of chunks with nothing to exchange
without scanning any cells.

Must be called after [`sort_particles!`](@ref) and before [`exchange_particles!`](@ref),
as it relies on `gridsort.cell_counts` holding the prefix sum of the per-cell particle counts
left there by the sort. If not called, one should set `occ_lo` to 1 and `occ_hi` to `n_cells`
in each chunk.

# Positional arguments
* `chunk_exchanger`: the `ChunkExchanger` instance in which to store the bounds
* `gridsort`: the `GridSortInPlace` used to sort the chunk's particles
* `pia`: the `ParticleIndexerArray` instance associated with the chunk
* `chunk_id`: the chunk for which to update the bounds
* `species`: the particle species for which the bounds are computed
"""
function update_occupancy_bounds!(chunk_exchanger, gridsort, pia, chunk_id, species)
    @inbounds n_tot = pia.n_total[species]
    cell_counts = gridsort.cell_counts

    if n_tot == 0
        @inbounds chunk_exchanger.occ_lo[chunk_id] = 1
        @inbounds chunk_exchanger.occ_hi[chunk_id] = 0
    else
        # after sorting, cell_counts[cell+1] holds the number of particles in cells 1:cell-1,
        # so the first entry reaching 1 (n_tot) is 2 past the first (last) occupied cell
        @inbounds chunk_exchanger.occ_lo[chunk_id] = searchsortedfirst(cell_counts, 1) - 2
        @inbounds chunk_exchanger.occ_hi[chunk_id] = searchsortedfirst(cell_counts, n_tot) - 2
    end
    return nothing
end

"""
    reset!(chunk_exchanger, chunk_id)

Reset all indexing of the entries `[chunk_id,:]` of a `chunk_exchanger`.

A newly constructed `ChunkExchanger` is already empty, and
[`sort_particles_after_exchange!`](@ref) clears the entries it reads, so in a time loop
that sorts every chunk after every exchange this does not need to be called at all.
It is intended for re-using a `ChunkExchanger` whose indexing was left in an
unknown state, e.g. after an exchange that was not followed by a re-sort.

# Positional arguments
* `chunk_exchanger`: the `ChunkExchanger` instance
* `chunk_id`: the chunk for which to reset indexing
"""
function reset!(chunk_exchanger, chunk_id)
    @inbounds for i in 1:chunk_exchanger.n_cells
        chunk_exchanger.start1[chunk_id, i] = 0
        chunk_exchanger.n_group1[chunk_id, i] = 0
        chunk_exchanger.start2[chunk_id, i] = 0
        chunk_exchanger.n_group2[chunk_id, i] = 0
    end
end

"""
    push_particles_to_end!(pv_i::ParticleVector{D}, pv_j::ParticleVector{D}, s2, e2, offset) where D

Copy the particles `pv_i[s2+offset:e2+offset]` to the positions `s2:e2` of `pv_j`,
taking new particles in `pv_j` from its buffer, and adding the freed particles of `pv_i`
to its buffer. This does not update any `ParticleIndexer`/`ParticleIndexerArray` indexing.

# Positional arguments
* `pv_i`: the source `ParticleVector`
* `pv_j`: the destination `ParticleVector`
* `s2`: first index in `pv_j` to write to
* `e2`: last index in `pv_j` to write to
* `offset`: offset between the indices in `pv_j` and the indices of the particles in `pv_i`
"""
@inline function push_particles_to_end!(pv_i::ParticleVector{D}, pv_j::ParticleVector{D}, s2, e2, offset) where D
    # hoisted out of the loop: `pv_i`/`pv_j` are mutable, so these would otherwise be
    # re-loaded on every iteration, as would `nbuffer` be read-modify-written in place
    index_i = pv_i.index
    index_j = pv_j.index
    buffer_i = pv_i.buffer
    nbuffer_i = pv_i.nbuffer
    particles_i = pv_i.particles
    particles_j = pv_j.particles

    @inbounds for pid in s2:e2
        update_particle_buffer_new_particle!(pv_j, pid)

        true_index_i = index_i[pid + offset]
        p_j = particles_j[index_j[pid]]
        p_i = particles_i[true_index_i]

        p_j.w = p_i.w
        p_j.v = p_i.v
        p_j.x = p_i.x

        nbuffer_i += 1
        buffer_i[nbuffer_i] = true_index_i
    end
    pv_i.nbuffer = nbuffer_i
    return nothing
end

"""
    push_particles!(chunk_exchanger, particles_chunks::Vector{Vector{ParticleVector{D}}}, pia_chunks, species, i, j, offset_ij, s_ci_ij2, e_ci_ij) where D

Pushes particles of the specified `species` from chunk `i` to the end of chunk `j`,
handling partial (unfinished) swaps.

This function finalizes a particle redistribution process between two spatial
chunks in a parallel or chunked particle simulation.
It ensures all remaining particles from chunk `i` that need to be moved to chunk `j`
are properly transferred and indexed, even if the previous transfer process via swapping particles was interrupted mid-cell.
This also updates the `buffer` of the source `ParticleVector` as particles are removed from it.

# Positional arguments
* `chunk_exchanger`: the `ChunkExchanger` instance to track post-swap and post-push indices
* `particles_chunks`: Vector of Vector of `ParticleVector` (per chunk and per species, i.e.
    `particles_chunks[chunk_id][species]` is the correct order of access)
* `pia_chunks`: Vector of `ParticleIndexerArray` instances for each chunk
* `species`: the particle species being redistributed
* `i`: Source chunk index
* `j`: Destination chunk index
* `offset_ij`: number of particles already transferred from a cell where swapping
    was performed but incomplete (0 if no incomplete swapping)
* `s_ci_ij2`: cell index in chunk `i` where the swapping was interrupted, or where the push
    should start in case no swapping was performed
* `e_ci_ij`: final cell index to transfer from chunk `i` to chunk `j`
"""
function push_particles!(chunk_exchanger, particles_chunks::Vector{Vector{ParticleVector{D}}}, pia_chunks, species, i, j, offset_ij, s_ci_ij2, e_ci_ij) where D
    # push remaining particles from chunk i to end of chunk j
    # we take care of s_ci_ij2 separately
    # because we might have stopped the swapping process in the middle of the cell
    # this is the cell of particles in chunk i where we stopped swapping
    cell = s_ci_ij2

    ce_start2 = chunk_exchanger.start2
    ce_n_group2 = chunk_exchanger.n_group2

    @inbounds pia_chunk_j = pia_chunks[j]
    @inbounds indexer = pia_chunks[i].indexer[cell,species]

    @inbounds pv_i = particles_chunks[i][species]
    @inbounds pv_j = particles_chunks[j][species]

    # println("Push from $i to $j")
    # println("Starting from cell $s_ci_ij2, already pushed $offset_ij from there")
    # println("In total it has $(pia_chunks[i].indexer[cell,species].n_group1) particles that need to be pushed")
    # we couldn't swap all particles and stopped the swapping process in the middle of the cell
    if (indexer.n_group1) > 0 && (offset_ij < indexer.n_group1)
        n_leftover = indexer.n_group1 - offset_ij

        # println("Pushing from cell where swap was broken off, $n_leftover left")
        # we write to [i,cell], not [j,cell]
        # because otherwise we might overwrite this when we transfer from another chunk to j
        # so the exchanger tracks from which chunk the particles came
        @inbounds s2 = pia_chunk_j.n_total[species] + 1
        @inbounds pia_chunk_j.n_total[species] += n_leftover
        e2 = s2 + n_leftover - 1

        @inbounds ce_start2[i,cell] = s2
        @inbounds ce_n_group2[i,cell] = n_leftover

        offset = -s2 + indexer.start1 + offset_ij

        # println("will write to $s2:$e2 in chunk $j")

        push_particles_to_end!(pv_i, pv_j, s2, e2, offset)

        indexer.n_local = 0
        indexer.n_group1 = 0
        indexer.start1 = 0
        indexer.end1 = -1
    end
    @inbounds for cell in s_ci_ij2+1:e_ci_ij
        indexer = pia_chunks[i].indexer[cell,species]
        n_push = indexer.n_group1
        if n_push > 0
            s2 = pia_chunk_j.n_total[species] + 1
            pia_chunk_j.n_total[species] += n_push
            e2 = s2 + n_push - 1

            ce_start2[i,cell] = s2
            ce_n_group2[i,cell] = n_push

            offset = -s2 + indexer.start1

            push_particles_to_end!(pv_i, pv_j, s2, e2, offset)

            indexer.n_local = 0
            indexer.n_group1 = 0
            indexer.start1 = 0
            indexer.end1 = -1
        end
    end
end

"""
    update_swap_indexing!(chunk_exchanger, pia_chunks, species, i, j, s_ci_ij, e_ci_ij, s_ji, n_swap)

Updates index bookkeeping in `chunk_exchanger` for particles that were swapped from chunk `i` to chunk `j`.

After particles are swapped, this function records where in chunk
`j` those particles from chunk `i` were placed, cell by cell. It does not perform the actual swapping.
For correct bookkeeping it thus needs to be called twice, with all arguments dependent on `i` and `j`
symmetrically swapped.

# Positional arguments
* `chunk_exchanger`: the `ChunkExchanger` instance to track post-swap and post-push indices
* `pia_chunks`: Vector of `ParticleIndexerArray` instances for each chunk
* `species`: the particle species being redistributed
* `i`: source chunk index
* `j`: destination chunk index
* `s_ci_ij`: start cell index in chunk `i` for the swap
* `e_ci_ij`: end cell index in chunk `i` for the swap
* `s_ji`: starting index in chunk `j`'s particle array where swapped particles from chunk `i` were placed
* `n_swap`: total number of particles that were successfully swapped

# Returns
* `s_ci_ij2`: the last cell index in the iteration (i.e., where the last particle was swapped)
* `offset_ij`: number of particles swapped from `s_ci_ij2`, useful if the swap was interrupted mid-cell
"""
function update_swap_indexing!(chunk_exchanger, pia_chunks, species, i, j, s_ci_ij, e_ci_ij, s_ji, n_swap)
    # update indexing in chunk_exchanger after we've swapped particles
    # this is an update for particles moved from chunk i to chunk j
    # both are declared Int64 so that the returned tuple stays concretely typed on 32-bit
    # builds, where a bare `0` literal is an Int32 while the indexing fields are Int64
    s_ci_ij2::Int64 = s_ci_ij
    offset_ij::Int64 = 0

    ce_start1 = chunk_exchanger.start1
    ce_n_group1 = chunk_exchanger.n_group1

    # iterate over cells that belong to chunk j but where particles are present in chunk i
    @inbounds for cell in s_ci_ij:e_ci_ij
        # check out how many particles we actually can use (if there are any)
        indexer = pia_chunks[i].indexer[cell,species]
        offset = min(indexer.n_group1, n_swap)

        # we will swap all particles, so we can safely set this to 0/-1
        # if not, we use this data in the subsequent push of remaining particles
        # and reset indexing there 
        if offset == indexer.n_group1
            indexer.n_local = 0
            indexer.n_group1 = 0
            indexer.start1 = 0
            indexer.end1 = -1
        end

        if offset > 0
            # println("$i -> $j update in $cell from $s_ci_ij:$e_ci_ij: $s_ji / $offset")
            # we stored the starting index of where we started swapping particles
            # so the particles written to chunk j during the swap
            # will start at s_ji and continue
            # particles come from chunk i into cell that belongs to chunk j
            ce_start1[i,cell] = s_ji
            ce_n_group1[i,cell] = offset
            s_ji += offset - 1
            n_swap -= offset

            # store the cell where we currently are
            s_ci_ij2 = cell
            
            # this tells us how many particles from the group we actually managed
            # to swap and whether there are any remaining
            # this is the offset for particles sent from j to i
            offset_ij = offset

            if n_swap <= 0
                break
            end

            # so that start of next cell != end of current cell
            s_ji += 1
        end
    end
    return s_ci_ij2, offset_ij
end


"""
    exchange_particles!(chunk_exchanger, particles_chunks::Vector{Vector{ParticleVector{D}}}, pia_chunks, cell_chunks, species, i, j) where D

Redistribute particles between chunks `i` and `j` based on their spatial cell ownership.

This function ensures each particle resides in the chunk responsible for its current cell.
It performs symmetric swaps when possible, and pushes remaining particles if needed.
The indexing metadata (`chunk_exchanger` and `pia_chunks`) is updated accordingly,
and particles pushed to another chunk (not swapped) are added to the buffer for future re-use.
The particles before the start of the re-distribution need to be sorted,
so that no particles are indexed by the `start2:end2` part of a `ParticleIndexer`. 
After the operation, the `n_total[species]` value of `pia_chunks[chunk_id]`
will not include particles that were pushed to another chunk (the appropriate
`n_group1`, `start1`, `end1` values will be set to 0, 0, -1). However
indexing should not be relied on until particles are re-sorted, see (`sort_particles_after_exchange!`)[@ref].

# Positional arguments
* `chunk_exchanger`: the `ChunkExchanger` instance to track post-swap and post-push indices
* `particles_chunks`: Vector of Vector of `ParticleVector` (per chunk and per species, i.e.
    `particles_chunks[chunk_id][species]` is the correct order of access)
* `pia_chunks`: Vector of `ParticleIndexerArray` instances for each chunk
* `cell_chunks`: cell ownership list for each chunk, i.e. `cell_chunks[chunk_id]` is a list
    of cells belonging to chunk `chunk_id`; the cells within `cell_chunks[chunk_id]` should
    be ordered in increasing order and be continuous:
    i.e. `cell_chunks[chunk_id][i] == cell_chunks[chunk_id][i-1] + 1`
* `species`: the particle species being redistributed
* `i`: index of first chunk
* `j`: index of second chunk
"""
function exchange_particles!(chunk_exchanger, particles_chunks::Vector{Vector{ParticleVector{D}}}, pia_chunks, cell_chunks, species, i, j) where D
    # ChunkSplitters' chunk collections are only indexable by the native Int, so the chunk
    # ids are converted here and any wider integer type can be passed in
    chunk_i = Int(i)
    chunk_j = Int(j)

    # the cells of chunk j that chunk i could possibly hold particles for: the cells owned by j,
    # restricted to the cells in which chunk i holds anything at all. If this range is empty,
    # chunk i has nothing for chunk j and no cell has to be looked at
    @inbounds lo_ij = max(first(cell_chunks[chunk_j]), chunk_exchanger.occ_lo[i])
    @inbounds hi_ij = min(last(cell_chunks[chunk_j]), chunk_exchanger.occ_hi[i])

    # find how many particles need to be transferred from i to j
    # we find first index of particles in chunk i that belong to a cell
    # assigned to chunk j
    # the particle and cell indices are declared Int64 to match the `ParticleIndexer` fields
    # and the occupancy bounds: on 32-bit builds a bare `0` literal is an Int32 and the
    # variable would infer as Union{Int32,Int64}, which boxes and allocates
    s_ij::Int64 = 0
    s_ci_ij::Int64 = 0 # index of the cell
    @inbounds for cj in lo_ij:hi_ij
        st = pia_chunks[i].indexer[cj, species].start1
        if st > 0
            s_ij = st
            s_ci_ij = cj
            break
        end
    end

    # we find last index of particles in chunk i that belong to a cell
    # assigned to chunk j
    e_ij::Int64 = -1
    e_ci_ij::Int64 = 0 # index of the cell
    @inbounds for cj in hi_ij:-1:lo_ij
        et = pia_chunks[i].indexer[cj, species].end1
        if et > 0
            e_ij = et
            e_ci_ij = cj
            break
        end
    end

    np_from_i_to_j = e_ij - s_ij + 1

    # now we do the same, but for particles in chunk j
    # that should be transferred to chunk i
    @inbounds lo_ji = max(first(cell_chunks[chunk_i]), chunk_exchanger.occ_lo[j])
    @inbounds hi_ji = min(last(cell_chunks[chunk_i]), chunk_exchanger.occ_hi[j])

    s_ji::Int64 = 0
    s_ci_ji::Int64 = 0 # index of the cell
    @inbounds for ci in lo_ji:hi_ji
        st = pia_chunks[j].indexer[ci, species].start1
        if st > 0
            s_ji = st
            s_ci_ji = ci
            break
        end
    end

    e_ji::Int64 = -1
    e_ci_ji::Int64 = 0 # index of the cell
    @inbounds for ci in hi_ji:-1:lo_ji
        et = pia_chunks[j].indexer[ci, species].end1
        if et > 0
            e_ji = et
            e_ci_ji = ci
            break
        end
    end

    np_from_j_to_i = e_ji - s_ji + 1

    if np_from_i_to_j <= 0 && np_from_j_to_i <= 0
        return nothing
    end

    # compute whether we need to increase sizes of the particle vectors
    # how many particles does chunk i receive
    chunk_i_increase = np_from_j_to_i > 0 ? np_from_j_to_i : 0
    # how many particles does chunk i send away
    chunk_i_increase = np_from_i_to_j > 0 ? chunk_i_increase - np_from_i_to_j : chunk_i_increase

    # how many particles does chunk j receive
    chunk_j_increase = np_from_i_to_j > 0 ? np_from_i_to_j : 0
    # how many particles does chunk j send away
    chunk_j_increase = np_from_j_to_i > 0 ? chunk_j_increase - np_from_j_to_i : chunk_j_increase

    # println("$i -> $j ", s_ij, " ", e_ij)
    # println("$j -> $i ", s_ji, " ", e_ji)

    # println("$i -> $j #: $(np_from_i_to_j)")
    # println("$j -> $i #: $(np_from_j_to_i)")

    # println("$i ++: $(chunk_i_increase)")
    # println("$j ++: $(chunk_j_increase)")

    # println("$i ntot: $(pia_chunks[i].n_total[species])")
    # println("$j ntot: $(pia_chunks[j].n_total[species])")
    # println(chunk_i_increase)

    if (length(particles_chunks[i][species]) < pia_chunks[i].n_total[species]+chunk_i_increase)
        resize!(particles_chunks[i][species],
                length(particles_chunks[i][species])+chunk_i_increase+DELTA_PARTICLES)
    end

    if (length(particles_chunks[j][species]) < pia_chunks[j].n_total[species]+chunk_j_increase)
        resize!(particles_chunks[j][species],
                length(particles_chunks[j][species])+chunk_j_increase+DELTA_PARTICLES)
    end

    # swap particles that can be swapped
    n_swap = min(np_from_i_to_j, np_from_j_to_i)

    offset_ij::Int64 = 0
    offset_ji::Int64 = 0

    s_ci_ij2::Int64 = s_ci_ij
    s_ci_ji2::Int64 = s_ci_ji

    # println("n_swap = $n_swap")
    # now update chunk_exchanger indexing
    if n_swap > 0
        @inbounds pv_i = particles_chunks[i][species]
        @inbounds pv_j = particles_chunks[j][species]

        @inbounds for nsw in 1:n_swap
            swap_particles!(pv_i, pv_j, s_ij+nsw-1, s_ji+nsw-1)
            # println("Swapping $(s_ij+nsw-1) from $i with $(s_ji+nsw-1) from $j")
        end

        # update for particles sent from i to j
        s_ci_ij2, offset_ij = update_swap_indexing!(chunk_exchanger, pia_chunks, species,
                                                    i, j, s_ci_ij, e_ci_ij, s_ji, n_swap)
        # update for particles sent from j to i 
        s_ci_ji2, offset_ji = update_swap_indexing!(chunk_exchanger, pia_chunks, species,
                                                    j, i, s_ci_ji, e_ci_ji, s_ij, n_swap)

        np_from_i_to_j -= n_swap
        np_from_j_to_i -= n_swap
    end

    # move remaining particles that did not fit into the swap
    # only one case is possible, because n_swap = min(np_from_i_to_j, np_from_j_to_i)
    # so one of those will be 0
    if np_from_i_to_j > 0
        # println("pushing from $i to $j: $np_from_i_to_j to push")
        # push remaining particles from i to end of j
        push_particles!(chunk_exchanger, particles_chunks, pia_chunks, species, i, j, offset_ij, s_ci_ij2, e_ci_ij)
    elseif np_from_j_to_i > 0
        # println("pushing from $j to $i: $np_from_j_to_i to push")
        # push remaining particles from j to end of i
        push_particles!(chunk_exchanger, particles_chunks, pia_chunks, species, j, i, offset_ji, s_ci_ji2, e_ci_ji)
    end
end

"""
    exchange_particles!(chunk_exchanger, particles_chunks::Vector{Vector{ParticleVector{D}}}, pia_chunks, cell_chunks, species) where D

Redistribute particles between chunks based on their spatial cell ownership.

This function ensures each particle resides in the chunk responsible for its current cell.
It performs symmetric swaps when possible, and pushes remaining particles if needed.
The indexing metadata (`chunk_exchanger` and `pia_chunks`) is updated accordingly,
and particles pushed to another chunk (not swapped) are added to the buffer for future re-use.
The particles before the start of the re-distribution need to be sorted,
so that no particles are indexed by the `start2:end2` part of a `ParticleIndexer`. 
After the operation, the `n_total[species]` value of `pia_chunks[chunk_id]`
will not include particles that were pushed to another chunk (the appropriate
`n_group1`, `start1`, `end1` values will be set to 0, 0, -1). However
indexing should not be relied on until particles are re-sorted, see (`sort_particles_after_exchange!`)[@ref].

# Positional arguments
* `chunk_exchanger`: the `ChunkExchanger` instance to track post-swap and post-push indices
* `particles_chunks`: Vector of Vector of `ParticleVector` (per chunk and per species, i.e.
    `particles_chunks[chunk_id][species]` is the correct order of access)
* `pia_chunks`: Vector of `ParticleIndexerArray` instances for each chunk
* `cell_chunks`: cell ownership list for each chunk, i.e. `cell_chunks[chunk_id]` is a list
    of cells belonging to chunk `chunk_id`; the cells within `cell_chunks[chunk_id]` should
    be ordered in increasing order and be continuous:
    i.e. `cell_chunks[chunk_id][i] == cell_chunks[chunk_id][i-1] + 1`
* `species`: the particle species being redistributed
"""
function exchange_particles!(chunk_exchanger, particles_chunks::Vector{Vector{ParticleVector{D}}}, pia_chunks, cell_chunks, species) where D
    n_chunks = length(cell_chunks)
    @inbounds for i in 1:n_chunks-1
        for j in i+1:n_chunks
            exchange_particles!(chunk_exchanger, particles_chunks, pia_chunks, cell_chunks, species, i, j)
        end
    end
end

"""
    sort_particles_after_exchange!(chunk_exchanger, gridsort, particles::ParticleVector{D}, pia, cell_chunk, species) where D

Restore indexing of a `ParticleVector` and the associated `ParticleIndexerArray`
after particles have been swapped and pushed between chunks.

# Positional arguments
* `chunk_exchanger`: the `ChunkExchanger` instance used to track post-swap and post-push indices
* `gridsort`: The `GridSortInPlace` associated with the chunk
* `particles`: the `ParticleVector` for which to restore the indexing
* `pia`: the `ParticleIndexerArray` instances associated with the chunk
* `cell_chunk`: list or range of cells belonging to the chunk
* `species`: the particle species being for which the indexing is being restored
"""
function sort_particles_after_exchange!(chunk_exchanger, gridsort, particles::ParticleVector{D}, pia, cell_chunk, species) where D
    @inbounds n_tot = pia.n_total[species] 
    @inbounds if n_tot > length(gridsort.sorted_indices)
        resize!(gridsort.sorted_indices, n_tot + DELTA_PARTICLES)
    end
    
    ci = 0
    n_tot = 0
    offset = 0

    cell_counts = gridsort.cell_counts
    sorted_indices = gridsort.sorted_indices
    p_index = particles.index

    n_chunks = chunk_exchanger.n_chunks
    ce_start1 = chunk_exchanger.start1
    ce_n_group1 = chunk_exchanger.n_group1
    ce_start2 = chunk_exchanger.start2
    ce_n_group2 = chunk_exchanger.n_group2

    @inbounds for cell in cell_chunk
        indexer = pia.indexer[cell, species]
        cc = indexer.n_group1

        if cc > 0
            unsafe_copyto!(sorted_indices, ci+1, p_index, indexer.start1, cc)
            ci += cc
        end

        # for i in s1:e1
        #     ci += 1
        #     sorted_indices[ci] = p_index[i]
        # end

        for chunk_id in 1:n_chunks
            ng1 = ce_n_group1[chunk_id,cell]
            ng2 = ce_n_group2[chunk_id,cell]
            cc += ng1 + ng2

            # entries are cleared as they are read, so that the exchanger table stays
            # zeroed for the next timestep without a separate sweep over all cells
            if ng1 > 0
                unsafe_copyto!(sorted_indices, ci+1, p_index, ce_start1[chunk_id,cell], ng1)
                ci += ng1

                ce_start1[chunk_id,cell] = 0
                ce_n_group1[chunk_id,cell] = 0
            end

            if ng2 > 0
                unsafe_copyto!(sorted_indices, ci+1, p_index, ce_start2[chunk_id,cell], ng2)
                ci += ng2

                ce_start2[chunk_id,cell] = 0
                ce_n_group2[chunk_id,cell] = 0
            end
        end
        n_tot += cc

        indexer.n_group1 = cc
        indexer.n_local = cc

        if indexer.n_group1 > 0
            indexer.start1 = offset + 1
            indexer.end1 = offset + cc
        else
            indexer.start1 = 0
            indexer.end1 = -1
        end

        offset += cc

        cell_counts[cell] = cc

        indexer.n_group2 = 0
        indexer.start2 = 0
        indexer.end2 = -1
    end

    @inbounds pia.n_total[species] = n_tot
    @inbounds pia.index_last[species] = n_tot

    unsafe_copyto!(p_index, 1, sorted_indices, 1, n_tot)
end

"""
    generate_1_factorization(N_chunks)

Construct a `Vector{Vector{Tuple{Int,Int}}}` with the following properties:
* tuple elements `i` and `j` range from `1` to `N_chunks`
* tuples `(i,j)` and `(j,i)` are considered equivalent
* each tuple `(i,j)` (up to equivalency) appears in the result exactly once
* in each `Vector` of tuples all numbers are unique
* tuples `(i,i)` do not appear
* the `Vector`s of tuples are as long as possible.

This corresponds to a 1-factorization of a complete graph; the resulting
`Vector` of `Vector`s of `Tuple`s can be iterated over, and particle exchange
can be performed between chunks listed in the `Tuple`s in a given `Vector`
using multi-threaded, since each chunk appears at most once in a given vector.
This allows to multi-thread the particle exchange step.


# Positional arguments
* `N_chunks`: number of chunks

# Returns
A `Vector` of `Vector`s of `Tuple`s corresponding to the 1-factorization of a complete
graph with `N_chunks` vertices.
"""
function generate_1_factorization(N_chunks)
    pairs = [(i, j) for i in 1:N_chunks for j in i+1:N_chunks]

    list_of_lists::Vector{Vector{Tuple{Int,Int}}} = []

    for (i, j) in pairs
        assigned = false
        for inner_list in list_of_lists
            used_i = any(p -> i in p, inner_list)
            used_j = any(p -> j in p, inner_list)
            if !(used_i || used_j)
                push!(inner_list, (i, j))
                assigned = true
                break
            end
        end
        if !assigned
            push!(list_of_lists, [(i, j)])
        end
    end

    return list_of_lists
end
