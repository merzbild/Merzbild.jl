"""
    GridSortInPlace

Struct for in-place sorting of particles.

# Fields
* `cell_counts`: vector to store the number of particles in each cell + number of particles in all previous cells
* `sorted_indices`: vector to store sorted particle indices
"""
mutable struct GridSortInPlace
    cell_counts::Vector{Int64}  # n_cells + 1
    sorted_indices::Vector{Int64}
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
GridSortInPlace(grid, n_particles::Integer) = GridSortInPlace(zeros(Int64, grid.n_cells + 1), zeros(Int64, n_particles))

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
    @inbounds n_tot = pia.n_total[species] 
    @inbounds if n_tot > length(gridsort.sorted_indices)
        resize!(gridsort.sorted_indices, n_tot + DELTA_PARTICLES)
    end

    gridsort.cell_counts[:] .= 0

    @inbounds if !pia.contiguous[species]
        squash_pia!(particles, pia, species)
    end

    @simd for i in 1:n_tot
        @inbounds newcell = get_cell(grid, particles[i].x)
        @inbounds particles.cell[i] = newcell
        @inbounds gridsort.cell_counts[newcell+1] += 1
    end

    for cell in 1:grid.n_cells
        @inbounds gridsort.cell_counts[cell+1] = gridsort.cell_counts[cell+1] + gridsort.cell_counts[cell]

        @inbounds cell_start = gridsort.cell_counts[cell] + 1
        @inbounds cell_np = gridsort.cell_counts[cell+1] - gridsort.cell_counts[cell]
        @inbounds cell_end = gridsort.cell_counts[cell+1]

        @inbounds pia.indexer[cell,species].start2 = 0
        @inbounds pia.indexer[cell,species].end2 = -1
        @inbounds pia.indexer[cell,species].n_group2 = 0

        if cell_np > 0
            @inbounds pia.indexer[cell,species].start1 = cell_start
            @inbounds pia.indexer[cell,species].end1 = cell_end
        else
            # this is done so that we can safely write for i in e1:s1 without worrying about accessing particles at index 0
            @inbounds pia.indexer[cell,species].start1 = 0
            @inbounds pia.indexer[cell,species].end1 = -1
        end
        @inbounds pia.indexer[cell,species].n_group1 = cell_np
        @inbounds pia.indexer[cell,species].n_local = cell_np
    end

    for i in n_tot:-1:1
        @inbounds curr_cell = particles.cell[i]
        @inbounds gridsort.sorted_indices[gridsort.cell_counts[curr_cell+1]] = particles.index[i]

        @inbounds gridsort.cell_counts[curr_cell+1] -= 1
    end

    @inbounds @simd for i in 1:n_tot
        @inbounds particles.index[i] = gridsort.sorted_indices[i]
    end

    @inbounds pia.contiguous[species] = true
end