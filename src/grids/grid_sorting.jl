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
    GridSortInPlace(n_cells::Integer, n_particles::Integer)

Create a `GridSortInPlace` instance given a number of grid cells and number of particles.

# Positional arguments
* `n_cells`: the number of grid cells
* `n_particles`: the (expected) number of particles in the simulation
    (to pre-allocate the `sorted_indices` vector) - it is recommended
    to set this to the maximum expected number of particles in the simulation to avoid resizing of arrays
    during a simulation
"""
GridSortInPlace(n_cells::Integer, n_particles::Integer) = GridSortInPlace(zeros(Int64, n_cells + 1), zeros(Int64, n_particles))

@doc """
    GridSortInPlace(grid::G, n_particles::Integer) where {G<:AbstractGrid}

Create a `GridSortInPlace` instance given a grid and number of particles.

# Positional arguments
* `grid`: the grid on which to sort the particles
* `n_particles`: the (expected) number of particles in the simulation
    (to pre-allocate the `sorted_indices` vector) - it is recommended
    to set this to the maximum expected number of particles in the simulation to avoid resizing of arrays
    during a simulation
"""
GridSortInPlace(grid::G, n_particles::Integer) where {G<:AbstractGrid} = GridSortInPlace(grid.n_cells, n_particles)

"""
    sort_particles!(gridsort::GridSortInPlace, grid, particles::ParticleVector{D}, pia, species) where D

Sort particles on a grid using an in-place sorting algorithm. The `pia` instance is allowed to
have non-contiguous indices (arising for example from merging). This function
assumes that at the start of the sorting, it is **not known** in which cell each particle is located,
and therefore the cell for each particle has to be determined (by calling `get_cell`).

# Positional arguments
* `gridsort`: the `GridSortInPlace` structure
* `grid`: the grid (should have an `n_cells` field, and a `get_cell` function has to be defined for the grid type)
* `particles`: the `ParticleVector` of particles to be sorted
* `pia`: the `ParticleIndexerArray` instance
* `species`: the index of the species being sorted
"""
function sort_particles!(gridsort::GridSortInPlace, grid, particles::ParticleVector{D}, pia, species) where D
    cell_counts = gridsort.cell_counts
    sorted_indices = gridsort.sorted_indices
    p_index = particles.index
    p_cell = particles.cell
    p_particles = particles.particles

    @inbounds n_tot = pia.n_total[species] 
    @inbounds if n_tot > length(sorted_indices)
        resize!(sorted_indices, n_tot + DELTA_PARTICLES)
    end

    fill!(cell_counts, 0)

    @inbounds if !pia.contiguous[species]
        squash_pia!(particles, pia, species)
    end

    @inbounds for i in 1:n_tot
        newcell = get_cell(grid, p_particles[p_index[i]].x)
        p_cell[i] = newcell
        cell_counts[newcell+1] += 1
    end

    n_cells = grid.n_cells
    @inbounds for cell in 1:n_cells
        cell_start = cell_counts[cell] + 1
        cell_np = cell_counts[cell+1]
        cell_counts[cell+1] = cell_counts[cell+1] + cell_counts[cell]
        cell_end = cell_counts[cell+1]

        indexer = pia.indexer[cell,species]

        indexer.start2 = 0
        indexer.end2 = -1
        indexer.n_group2 = 0

        if cell_np > 0
            indexer.start1 = cell_start
            indexer.end1 = cell_end
        else
            # this is done so that we can safely write for i in e1:s1 without worrying about accessing particles at index 0
            indexer.start1 = 0
            indexer.end1 = -1
        end
        indexer.n_group1 = cell_np
        indexer.n_local = cell_np
    end

    @inbounds for i in n_tot:-1:1
        curr_cell = p_cell[i]
        sorted_indices[cell_counts[curr_cell+1]] = p_index[i]

        cell_counts[curr_cell+1] -= 1
    end

    unsafe_copyto!(p_index, 1, sorted_indices, 1, n_tot)

    @inbounds pia.contiguous[species] = true
end

"""
    sort_particles!(gridsort::GridSortInPlace, particles, pia, species)

Sort particles on a grid using an in-place sorting algorithm. The `pia` instance is allowed to
have non-contiguous indices (arising for example from merging). This function
assumes that at the start of the sorting, it is **known** in which cell each particle is located.

# Positional arguments
* `gridsort`: the `GridSortInPlace` structure
* `particles`: the `ParticleVector` of particles to be sorted
* `pia`: the `ParticleIndexerArray` instance
* `species`: the index of the species being sorted
"""
function sort_particles!(gridsort::GridSortInPlace, particles, pia, species)
    cell_counts = gridsort.cell_counts
    sorted_indices = gridsort.sorted_indices
    p_index = particles.index
    p_cell = particles.cell

    n_cells = pia.n_cells
    @inbounds n_tot = pia.n_total[species] 
    @inbounds if n_tot > length(sorted_indices)
        resize!(sorted_indices, n_tot + DELTA_PARTICLES)
    end

    fill!(cell_counts, 0)

    @inbounds if !pia.contiguous[species]
        squash_pia!(particles, pia, species)
    end

    @inbounds for i in 1:n_tot
        cell_counts[p_cell[i]+1] += 1
    end

    @inbounds for cell in 1:n_cells
        cell_start = cell_counts[cell] + 1
        cell_np = cell_counts[cell+1]
        cell_counts[cell+1] = cell_counts[cell+1] + cell_counts[cell]
        cell_end = cell_counts[cell+1]

        indexer = pia.indexer[cell,species]

        indexer.start2 = 0
        indexer.end2 = -1
        indexer.n_group2 = 0

        if cell_np > 0
            indexer.start1 = cell_start
            indexer.end1 = cell_end
        else
            # this is done so that we can safely write for i in e1:s1 without worrying about accessing particles at index 0
            indexer.start1 = 0
            indexer.end1 = -1
        end
        indexer.n_group1 = cell_np
        indexer.n_local = cell_np
    end

    @inbounds for i in n_tot:-1:1
        curr_cell = p_cell[i]
        sorted_indices[cell_counts[curr_cell+1]] = p_index[i]

        cell_counts[curr_cell+1] -= 1
    end

    unsafe_copyto!(p_index, 1, sorted_indices, 1, n_tot)

    @inbounds pia.contiguous[species] = true
end