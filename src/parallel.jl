"""
    ChunkExchanger

Struct used to organize exchange of particles between independent ParticleVector instances for
chunked multi-threaded simulations. It is assumed that the cell indices within each chunk are contiguous,
i.e. `chunk[i+1] = chunk[i]+1`.
The `n_local` field of `ParticleIndexer` is left unused.

Group 1 of `ParticleIndexer[chunk_id,cell]` holds the range of particles that
belong to cell `cell` that came from a particle chunk `chunk_id` via swapping.

Group 2 of `ParticleIndexer[chunk_id,cell]` holds the range of particles that
belong to cell `cell` that came from a particle chunk `chunk_id` via pushing, i.e.
they have been added to the end of the particle array.

# Fields
* `n_chunks`: number of chunks used in the simulation
* `n_cells`: number of grid cells in the simulation
* `indexer`: array of `ParticleIndexer` instances of shape `(n_chunks, n_cells)`
"""
mutable struct ChunkExchanger
    n_chunks::Int32
    n_cells::Int64
    indexer::Array{ParticleIndexer, 2}  # n_chunks x n_cells
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
    indexer = Array{ParticleIndexer, 2}(undef, (n_chunks, n_cells))

    for j in 1:n_cells
        for i in 1:n_chunks
            @inbounds indexer[i, j] = ParticleIndexer()
        end
    end
    return ChunkExchanger(n_chunks, n_cells, indexer)
end

"""
    reset!(chunk_exchanger, chunk_id)

Reset all indexing of `chunk_exchanger.indexer[chunk_id,:]`.

# Positional arguments
* `chunk_exchanger`: the `ChunkExchanger` instance
* `chunk_id`: the chunk for which to reset indexing
"""
function reset!(chunk_exchanger, chunk_id)
    for i in 1:chunk_exchanger.n_cells
        # reset only the things we actually need
        @inbounds chunk_exchanger.indexer[chunk_id, i].n_group1 = 0
        @inbounds chunk_exchanger.indexer[chunk_id, i].start1 = 0
        @inbounds chunk_exchanger.indexer[chunk_id, i].end1 = -1
        @inbounds chunk_exchanger.indexer[chunk_id, i].n_group2 = 0
        @inbounds chunk_exchanger.indexer[chunk_id, i].start2 = 0
        @inbounds chunk_exchanger.indexer[chunk_id, i].end2 = -1
    end
end

"""
    swap_particles!(pv1, pv2, i, j)

Swap particles `pv1[i]` and `pv2[j]` in two `ParticleVector` instances. This does not update any
associated indices or buffers.

# Positional arguments
* `pv1`: the first `ParticleVector`
* `pv2`: the second `ParticleVector`
* `i`: index of the particle in `pv1`
* `j`: index of the particle in `pv2`
"""
function swap_particles!(pv1, pv2, i, j)
    # TODO: check if inlining speeds things up or slows them down, doesn't seem to have an impact
    @inbounds tmp_w = pv1[i].w
    @inbounds tmp_v = pv1[i].v
    @inbounds tmp_x = pv1[i].x

    @inbounds pv1[i].w = pv2[j].w
    @inbounds pv1[i].v = pv2[j].v
    @inbounds pv1[i].x = pv2[j].x

    @inbounds pv2[j].w = tmp_w
    @inbounds pv2[j].v = tmp_v
    @inbounds pv2[j].x = tmp_x
end

"""
    push_particles!(chunk_exchanger, particles_chunks, pia_chunks, species, i, j, offset_ij, s_ci_ij2, e_ci_ij)

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
function push_particles!(chunk_exchanger, particles_chunks, pia_chunks, species, i, j, offset_ij, s_ci_ij2, e_ci_ij)
    # push remaining particles from chunk i to end of chunk j
    # we take care of s_ci_ij2 separately
    # because we might have stopped the swapping process in the middle of the cell
    # this is the cell of particles in chunk i where we stopped swapping
    cell = s_ci_ij2

    # println("Push from $i to $j")
    # println("Starting from cell $s_ci_ij2, already pushed $offset_ij from there")
    # println("In total it has $(pia_chunks[i].indexer[cell,species].n_group1) particles that need to be pushed")
    # we couldn't swap all particles and stopped the swapping process in the middle of the cell
    @inbounds if (pia_chunks[i].indexer[cell,species].n_group1) > 0 && (offset_ij < pia_chunks[i].indexer[cell,species].n_group1)
        n_leftover = pia_chunks[i].indexer[cell,species].n_group1 - offset_ij

        # println("Pushing from cell where swap was broken off, $n_leftover left")
        # we write to [i,cell], not [j,cell]
        # because otherwise we might overwrite this when we transfer from another chunk to j
        # so the exchanger tracks from which chunk the particles came
        @inbounds chunk_exchanger.indexer[i,cell].start2 = pia_chunks[j].n_total[species] + 1
        @inbounds pia_chunks[j].n_total[species] += n_leftover
        @inbounds chunk_exchanger.indexer[i,cell].end2 = pia_chunks[j].n_total[species]
        @inbounds chunk_exchanger.indexer[i,cell].n_group2 = n_leftover

        @inbounds s2 = chunk_exchanger.indexer[i,cell].start2
        @inbounds e2 = chunk_exchanger.indexer[i,cell].end2
        @inbounds  offset = -s2 + pia_chunks[i].indexer[cell,species].start1 + offset_ij

        # println("will write to $s2:$e2 in chunk $j")

        @inbounds for pid in s2:e2
            # println("Writing to $pid in $j using particle $(pid + offset) from $i")
            # move particle to chunk j
            update_particle_buffer_new_particle!(particles_chunks[j][species], pid)
            particles_chunks[j][species][pid].w = particles_chunks[i][species][pid + offset].w
            particles_chunks[j][species][pid].v = particles_chunks[i][species][pid + offset].v
            particles_chunks[j][species][pid].x = particles_chunks[i][species][pid + offset].x

            particles_chunks[i][species].nbuffer += 1
            particles_chunks[i][species].buffer[particles_chunks[i][species].nbuffer] =
                particles_chunks[i][species].index[pid + offset]
        end

        @inbounds pia_chunks[i].indexer[cell,species].n_local = 0
        @inbounds pia_chunks[i].indexer[cell,species].n_group1 = 0
        @inbounds pia_chunks[i].indexer[cell,species].start1 = 0
        @inbounds pia_chunks[i].indexer[cell,species].end1 = -1
    end
    @inbounds for cell in s_ci_ij2+1:e_ci_ij
        if pia_chunks[i].indexer[cell,species].n_group1 > 0
            chunk_exchanger.indexer[i,cell].start2 = pia_chunks[j].n_total[species] + 1
            pia_chunks[j].n_total[species] += pia_chunks[i].indexer[cell,species].n_group1
            chunk_exchanger.indexer[i,cell].end2 = pia_chunks[j].n_total[species]
            chunk_exchanger.indexer[i,cell].n_group2 = pia_chunks[i].indexer[cell,species].n_group1

            s2 = chunk_exchanger.indexer[i,cell].start2
            e2 = chunk_exchanger.indexer[i,cell].end2
            offset = -s2 + pia_chunks[i].indexer[cell,species].start1

            # write particles to chunk j
            for pid in s2:e2
                # move particle to chunk j
                update_particle_buffer_new_particle!(particles_chunks[j][species], pid)
                particles_chunks[j][species][pid].w = particles_chunks[i][species][pid + offset].w
                particles_chunks[j][species][pid].v = particles_chunks[i][species][pid + offset].v
                particles_chunks[j][species][pid].x = particles_chunks[i][species][pid + offset].x
                
                # update buffer in chunk i
                particles_chunks[i][species].nbuffer += 1
                particles_chunks[i][species].buffer[particles_chunks[i][species].nbuffer] =
                    particles_chunks[i][species].index[pid + offset]
            end
            pia_chunks[i].indexer[cell,species].n_local = 0
            pia_chunks[i].indexer[cell,species].n_group1 = 0
            pia_chunks[i].indexer[cell,species].start1 = 0
            pia_chunks[i].indexer[cell,species].end1 = -1
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
    s_ci_ij2 = s_ci_ij
    offset_ij = 0 

    # iterate over cells that belong to chunk j but where particles are present in chunk i
    @inbounds for cell in s_ci_ij:e_ci_ij
        # check out how many particles we actually can use (if there are any)
        offset = min(pia_chunks[i].indexer[cell,species].n_group1, n_swap)

        # we will swap all particles, so we can safely set this to 0/-1
        # if not, we use this data in the subsequent push of remaining particles
        # and reset indexing there 
        if offset == pia_chunks[i].indexer[cell,species].n_group1
            pia_chunks[i].indexer[cell,species].n_local = 0
            pia_chunks[i].indexer[cell,species].n_group1 = 0
            pia_chunks[i].indexer[cell,species].start1 = 0
            pia_chunks[i].indexer[cell,species].end1 = -1
        end

        if offset > 0
            # println("$i -> $j update in $cell from $s_ci_ij:$e_ci_ij: $s_ji / $offset")
            # we stored the starting index of where we started swapping particles
            # so the particles written to chunk j during the swap
            # will start at s_ji and continue
            # particles come from chunk i into cell that belongs to chunk j
            chunk_exchanger.indexer[i,cell].start1 = s_ji
            s_ji += offset - 1
            chunk_exchanger.indexer[i,cell].end1 = s_ji
            chunk_exchanger.indexer[i,cell].n_group1 = offset
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
    exchange_particles!(chunk_exchanger, particles_chunks, pia_chunks, cell_chunks, species)

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
function exchange_particles!(chunk_exchanger, particles_chunks, pia_chunks, cell_chunks, species)
    n_chunks = length(cell_chunks)
    @inbounds for i in 1:n_chunks-1
        for j in i+1:n_chunks
            # find how many particles need to be transferred from i to j
            # we find first index of particles in chunk i that belong to a cell
            # assigned to chunk j
            s_ij = 0
            s_ci_ij = 0 # index of the cell
            for cj in cell_chunks[j]
                if pia_chunks[i].indexer[cj, species].start1 > 0
                    s_ij = pia_chunks[i].indexer[cj, species].start1
                    s_ci_ij = cj
                    break
                end
            end

            # we find last index of particles in chunk i that belong to a cell
            # assigned to chunk j
            l_cj = length(cell_chunks[j])
            e_ij = -1
            e_ci_ij = 0 # index of the cell
            for cji in l_cj:-1:1
                cj = cell_chunks[j][cji]
                if pia_chunks[i].indexer[cj, species].end1 > 0
                    e_ij = pia_chunks[i].indexer[cj, species].end1
                    e_ci_ij = cj
                    break
                end
            end

            np_from_i_to_j = e_ij - s_ij + 1

            # now we do the same, but for particles in chunk j
            # that should be transferred to chunk i
            s_ji = 0
            s_ci_ji = 0 # index of the cell
            for ci in cell_chunks[i]
                if pia_chunks[j].indexer[ci, species].start1 > 0
                    s_ji = pia_chunks[j].indexer[ci, species].start1
                    s_ci_ji = ci
                    break
                end
            end

            l_ci = length(cell_chunks[i])
            e_ji = -1
            e_ci_ji = 0 # index of the cell
            for cij in l_ci:-1:1
                ci = cell_chunks[i][cij]
                if pia_chunks[j].indexer[ci, species].end1 > 0
                    e_ji = pia_chunks[j].indexer[ci, species].end1
                    e_ci_ji = ci
                    break
                end
            end

            np_from_j_to_i = e_ji - s_ji + 1

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

            offset_ij = 0
            offset_ji = 0

            s_ci_ij2 = s_ci_ij
            s_ci_ji2 = s_ci_ji

            # println("n_swap = $n_swap")
            # now update chunk_exchanger indexing
            if n_swap > 0
                for nsw in 1:n_swap
                    swap_particles!(particles_chunks[i][species], particles_chunks[j][species],
                                    s_ij+nsw-1, s_ji+nsw-1)
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
    end
end


"""
    sort_particles_after_exchange!(chunk_exchanger, gridsort, particles, pia, cell_chunk, species)

Restore indexing of a `ParticleVector` and the associated `ParticleIndexerArray`
after particles have been swapped and pushed between chunks.

# Positional arguments
* `chunk_exchanger`: the `ChunkExchanger` instance used to track post-swap and post-push indices
* `gridsort`: The `GridSortInPlace` associated with the chunk
* `particles_chunks`: the `ParticleVector` for which to restore the indexing
* `pia`: the `ParticleIndexerArray` instances associated with the chunk
* `cell_chunk`: list or range of cells belonging to the chunk
* `species`: the particle species being for which the indexing is being restored
"""
function sort_particles_after_exchange!(chunk_exchanger, gridsort, particles, pia, cell_chunk, species)
    @inbounds n_tot = pia.n_total[species] 
    @inbounds if n_tot > length(gridsort.sorted_indices)
        resize!(gridsort.sorted_indices, n_tot + DELTA_PARTICLES)
    end
    
    ci = 0
    n_tot = 0
    offset = 0
    @inbounds for cell in cell_chunk
        gridsort.cell_counts[cell] = pia.indexer[cell, species].n_group1

        s1 = pia.indexer[cell, species].start1
        e1 = pia.indexer[cell, species].end1

        for i in s1:e1
            ci += 1
            gridsort.sorted_indices[ci] = particles.index[i]
        end

        for chunk_id in 1:chunk_exchanger.n_chunks
            gridsort.cell_counts[cell] += chunk_exchanger.indexer[chunk_id,cell].n_group1
            gridsort.cell_counts[cell] += chunk_exchanger.indexer[chunk_id,cell].n_group2

            s1 = chunk_exchanger.indexer[chunk_id,cell].start1
            e1 = chunk_exchanger.indexer[chunk_id,cell].end1

            for i in s1:e1
                ci += 1
                gridsort.sorted_indices[ci] = particles.index[i]
            end

            s2 = chunk_exchanger.indexer[chunk_id,cell].start2
            e2 = chunk_exchanger.indexer[chunk_id,cell].end2

            for i in s2:e2
                ci += 1
                gridsort.sorted_indices[ci] = particles.index[i]
            end
        end
        n_tot += gridsort.cell_counts[cell]

        pia.indexer[cell,species].n_group1 = gridsort.cell_counts[cell]
        pia.indexer[cell,species].n_local = gridsort.cell_counts[cell]

        if pia.indexer[cell,species].n_group1 > 0
            pia.indexer[cell,species].start1 = offset + 1
            pia.indexer[cell,species].end1 = offset + gridsort.cell_counts[cell]
        else
            pia.indexer[cell,species].start1 = 0
            pia.indexer[cell,species].end1 = -1
        end

        offset += gridsort.cell_counts[cell]

        pia.indexer[cell,species].n_group2 = 0
        pia.indexer[cell,species].start2 = 0
        pia.indexer[cell,species].end2 = -1
    end

    @inbounds pia.n_total[species] = n_tot

    @inbounds for i in 1:n_tot
        particles.index[i] = gridsort.sorted_indices[i]
    end
end