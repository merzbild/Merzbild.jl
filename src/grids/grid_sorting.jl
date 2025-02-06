"""
    GridSortInPlace

Struct for in-place sorting of particles.

# Fields
* `cell_counts`: vector to store the number of particles in each cell + number of particles in all previous cells
* `sorted_indices`: vector to store sorted particle indices
* `non_contiguous_indices`: vector to store non-contiguous indices of the particles used in sorting in a contiguous array,
    in case the indices pointed to by a `ParticleIndexerArray` are not contiguous
"""
mutable struct GridSortInPlace
    cell_counts::Vector{Int64}  # n_cells + 1
    sorted_indices::Vector{Int64}
    non_contiguous_indices::Vector{Int64}
end

@doc """
    GridSortInPlace(grid, n_particles::Integer)

Create a `GridSortInPlace` instance given a grid and number of particles.

# Positional arguments
* `grid`: the grid on which to sort the particles
* `n_particles`: the (expected) number of particles in the simulation
    (to pre-allocate the `sorted_indices` vector) - it is recommended
    to set this to the maximum expected number of particles in the simulation to avoid resizing of arrays
    during a simulation
"""
GridSortInPlace(grid, n_particles::Integer) = GridSortInPlace(zeros(Int64, grid.n_cells + 1), zeros(Int64, n_particles),
                                                              zeros(Int64, n_particles))

"""
    sort_particles!(gridsort::GridSortInPlace, grid, particles, pia, species)

Sort particles on a grid using an in-place sorting algorithm. The `pia` instance is allowed to
have non-contiguous indices (arising for example from merging).

# Positional arguments
* `gridsort`: the `GridSortInPlace` structure
* `grid`: the grid (should have an `n_cells` field, and a `get_cell` function has to be defined for the grid type)
* `particles`: the `ParticleVector` of particles to be sorted
* `pia`: the `ParticleIndexerArray` instance
* `species`: the index of the species being sorted
"""
function sort_particles!(gridsort::GridSortInPlace, grid, particles, pia, species)
    if pia.n_total[species] > length(gridsort.sorted_indices)
        resize!(gridsort.sorted_indices, pia.n_total[species] + DELTA_PARTICLES)
        resize!(gridsort.non_contiguous_indices, pia.n_total[species] + DELTA_PARTICLES)
    end

    gridsort.cell_counts[:] .= 0

    if pia.contiguous[species]
        @inbounds @simd for i in 1:pia.n_total[species]
            newcell = get_cell(grid, particles[i].x)
            particles.cell[i] = newcell
            gridsort.cell_counts[newcell+1] += 1
        end
    else
        ci = 0
        for cell in 1:grid.n_cells
            if pia.indexer[cell,species].n_group1 > 0
                for nci in pia.indexer[cell,species].start1:pia.indexer[cell,species].end1
                    ci += 1
                    gridsort.non_contiguous_indices[ci] = nci
                end
            end

            if pia.indexer[cell,species].n_group2 > 0
                for nci in pia.indexer[cell,species].start2:pia.indexer[cell,species].end2
                    ci += 1
                    gridsort.non_contiguous_indices[ci] = nci
                end
            end
        end

        @inbounds @simd for i in 1:pia.n_total[species]
            nci = gridsort.non_contiguous_indices[i]
            newcell = get_cell(grid, particles[nci].x)
            particles.cell[i] = newcell
            gridsort.cell_counts[newcell+1] += 1
        end
    end

    @inbounds for cell in 1:grid.n_cells
        gridsort.cell_counts[cell+1] = gridsort.cell_counts[cell+1] + gridsort.cell_counts[cell]

        cell_start = gridsort.cell_counts[cell] + 1
        cell_np = gridsort.cell_counts[cell+1] - gridsort.cell_counts[cell]
        cell_end = gridsort.cell_counts[cell+1]

        pia.indexer[cell,species].start2 = 0
        pia.indexer[cell,species].end2 = -1
        pia.indexer[cell,species].n_group2 = 0

        if cell_np > 0
            pia.indexer[cell,species].start1 = cell_start
            pia.indexer[cell,species].end1 = cell_end
        else
            pia.indexer[cell,species].start1 = 0
            pia.indexer[cell,species].end1 = -1
        end
        pia.indexer[cell,species].n_group1 = cell_np
        pia.indexer[cell,species].n_local = cell_np
    end

    if pia.contiguous[species]
        @inbounds for i in pia.n_total[species]:-1:1
            curr_cell = particles.cell[i]
            gridsort.sorted_indices[gridsort.cell_counts[curr_cell+1]] = particles.index[i]

            gridsort.cell_counts[curr_cell+1] -= 1
        end
    else
        @inbounds for i in pia.n_total[species]:-1:1
            nci = gridsort.non_contiguous_indices[i]
            curr_cell = particles.cell[i]
            gridsort.sorted_indices[gridsort.cell_counts[curr_cell+1]] = particles.index[nci]

            gridsort.cell_counts[curr_cell+1] -= 1
        end
    end

    @inbounds @simd for i in 1:pia.n_total[species]
        particles.index[i] = gridsort.sorted_indices[i]
    end

    pia.contiguous[species] = true
end