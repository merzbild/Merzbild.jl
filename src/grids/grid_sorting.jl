"""
Struct for in-place sorting of particles
"""
mutable struct GridSortInPlace
    cell_counts::Vector{Int64}  # n_cells + 1
    sorted_indices::Vector{Int64}
end

"""
Create GridSortInPlace instance given a grid and number of particles
"""
GridSortInPlace(grid, n_particles::Integer) = GridSortInPlace(zeros(Int64, grid.n_cells + 1), zeros(Int64, n_particles))

"""
Sort particles on a grid
"""
function sort_particles!(gridsort::GridSortInPlace, grid, particles, pia, species)
    if pia.n_total[species] > length(gridsort.sorted_indices)
        resize!(gridsort.sorted_indices, pia.n_total[species] + DELTA_PARTICLES)
    end

    gridsort.cell_counts[:] .= 0

    @inbounds @simd for i in 1:pia.n_total[species]
        newcell = get_cell(grid, particles[i].x)
        particles.cell[i] = newcell
        gridsort.cell_counts[newcell+1] += 1
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

    @inbounds for i in pia.n_total[species]:-1:1
        curr_cell = particles.cell[i]
        gridsort.sorted_indices[gridsort.cell_counts[curr_cell+1]] = particles.index[i]

        gridsort.cell_counts[curr_cell+1] -= 1
    end

    @inbounds @simd for i in 1:pia.n_total[species]
        particles.index[i] = gridsort.sorted_indices[i]
    end
end