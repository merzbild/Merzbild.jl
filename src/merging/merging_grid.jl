using StaticArrays

"""
    GridCell

Struct for keeping track of merging-related quantities in a velocity grid cell.

# Fields
* `np`: number of particles in cell
* `w`: total computational weight of particles  in cell
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
mutable struct GridCell
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
end

"""
    GridN2Merge

Struct for merging using a grid in velocity space. Particles in each cell are merged down.
Particles outside of the grid are merged based on the octant they are in.
So for an N:2 merge one would expect at most `2*(Nx*Ny*Nz+8)` post-merge particles, where
`Nx`, `Ny`, `Nz` are the number of grid cells in each velocity direction.
The grid bounds in each direction can be computed using the mean thermal velocity of the particles and
the mean streaming velocity: `[v0-extent_multiplier*sqrt(2*k_B*T/m),v0+extent_multiplier*sqrt(2*k_B*T/m)]`.
Here `T` is the temperature of the species in question, `m` is the molecular mass, `v0` is the mean velocity
and `extent_multiplier` is a user-defined parameter (3.5 is a reasonable choice) defining the extent of the grid.

# Fields
* `Nx`: number of grid cells in x velocity direction
* `Ny`: number of grid cells in y velocity direction
* `Nz`: number of grid cells in z velocity direction
* `NyNz`: product of `Ny` and `Nz`
* `Ntotal`: total number of grid cells (equal to `Nx*Ny*Nz+8`, as we account for the external octants)
* `extent_multiplier`: the vector of factors by which to multiply the thermal velocity to determine the grid bounds
    in each velocity direction
* `extent_v_lower`: the lower bounds of the velocity grid in each velocity direction
* `extent_v_upper`: the upper bounds of the velocity grid in each velocity direction
* `extent_v_mid`: `(extent_v_lower + extent_v_upper)/2`
* `Δv`: the grid cell size in each velocity direction
* `Δv_inv`: the inverse grid cell size in each velocity direction
* `direction_vec`: used to store randomly sampled direction signs
* `cells`: vector of `GridCell` instances for each grid cell, as well as the external octants
"""
mutable struct GridN2Merge
    Nx::Int8
    Ny::Int8
    Nz::Int8
    NyNz::Int64
    Ntotal::Int64

    extent_multiplier::SVector{3,Float64}
    extent_v_lower::SVector{3,Float64}
    extent_v_upper::SVector{3,Float64}
    extent_v_mid::SVector{3,Float64}
    Δv::SVector{3,Float64}
    Δv_inv::SVector{3,Float64}
    direction_vec::SVector{3,Float64}

    cells::Vector{GridCell}

    @doc """
        GridN2Merge(Nx::Int, Ny::Int, Nz::Int, extent_multiplier::T) where T <: AbstractArray

    Create velocity grid-based merging.

    # Positional arguments
    * `Nx`: number of cells in vx direction
    * `Ny`: number of cells in vy direction
    * `Nz`: number of cells in vz direction
    * `extent_multiplier`: the vector of factors by which to multiply the thermal velocity to determine the grid bounds
    in each velocity direction
    """
    function GridN2Merge(Nx::Int, Ny::Int, Nz::Int, extent_multiplier::T) where T <: AbstractArray
        Ntotal = Nx * Ny * Nz + 8
        cells = Vector{GridCell}(undef, Nx * Ny * Nz + 8)

        for i in 1:Ntotal
            cells[i] = GridCell(0, 0.0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0, 0,
                                0.0, 0.0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0])
        end

        return new(Nx, Ny, Nz, Ny*Nz, Ntotal, extent_multiplier, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0],
                   [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], cells)
    end
end

"""
    GridN2Merge(N::Int, extent_multiplier::T) where T <: AbstractArray 

Create velocity grid-based merging with equal number of cells in each direction.

# Positional arguments
* `N`: number of cells in each velocity direction
* `extent_multiplier`: the vector of factors by which to multiply the thermal velocity to determine the grid bounds
in each velocity direction
"""
GridN2Merge(N::Int, extent_multiplier::T) where T <: AbstractArray = GridN2Merge(N, N, N, extent_multiplier)

"""
    GridN2Merge(Nx::Int, Ny::Int, Nz::Int, extent_multiplier::Float64)

Create velocity grid-based merging with equal multipliers in each direction.

# Positional arguments
* `Nx`: number of cells in vx direction
* `Ny`: number of cells in vy direction
* `Nz`: number of cells in vz direction
* `extent_multiplier`: the factor by which to multiply the thermal velocity to determine the grid bounds
in each velocity direction
"""
GridN2Merge(Nx::Int, Ny::Int, Nz::Int, extent_multiplier::Float64) = GridN2Merge(Nx, Ny, Nz, [extent_multiplier, extent_multiplier, extent_multiplier])

"""
    GridN2Merge(Nx::Int, Ny::Int, Nz::Int,
                extent_multiplier_x::Float64,
                extent_multiplier_y::Float64,
                extent_multiplier_z::Float64) 

Create velocity grid-based merging
* `Nx`: number of cells in vx direction
* `Ny`: number of cells in vy direction
* `Nz`: number of cells in vz direction
* `extent_multiplier_x`: the factor by which to multiply the thermal velocity to determine the grid bounds
in the x-velocity direction
* `extent_multiplier_y`: the factor by which to multiply the thermal velocity to determine the grid bounds
in the y-velocity direction
* `extent_multiplier_z`: the factor by which to multiply the thermal velocity to determine the grid bounds
in the z-velocity direction
"""
GridN2Merge(Nx::Int, Ny::Int, Nz::Int,
            extent_multiplier_x::Float64,
            extent_multiplier_y::Float64,
            extent_multiplier_z::Float64) = GridN2Merge(Nx, Ny, Nz, [extent_multiplier_x,
                                                                     extent_multiplier_y,
                                                                     extent_multiplier_z])

"""
    GridN2Merge(N::Int, extent_multiplier::Float64)

Create velocity grid-based merging with equal number of cells in each direction and equal multipliers
in each direction.

# Positional arguments
* `N`: number of cells in each velocity direction
* `extent_multiplier`: the factor by which to multiply the thermal velocity to determine the grid bounds
in each velocity direction
"""
GridN2Merge(N::Int, extent_multiplier::Float64) = GridN2Merge(N, N, N, extent_multiplier)

"""
    compute_velocity_extent!(merging_grid, cell, species, species_data, phys_props)

Compute extent of velocity grid based on temperature in the cell.

# Positional arguments:
* `merging_grid`: the grid merging (`GridN2Merge`) instance for which to compute the extent
* `cell`: the index of the cell in which the temperature was computed
* `species`: the index of the species for which the temperature was computed
* `species_data`: the array of `Species` data
* `phys_props`: the `PhysProps` instance containing the computed temperature
"""
function compute_velocity_extent!(merging_grid, cell, species, species_data, phys_props)
    dv = merging_grid.extent_multiplier .* sqrt.(2 * phys_props.T[cell, species] * k_B / species_data[species].mass)
    merging_grid.extent_v_lower = phys_props.v[:, cell, species] .- dv
    merging_grid.extent_v_upper = phys_props.v[:, cell, species] .+ dv
    merging_grid.extent_v_mid = phys_props.v[:, cell, species]
    merging_grid.Δv = SVector{3}(2 * dv[1] / merging_grid.Nx, 2 * dv[2] / merging_grid.Ny, 2 * dv[3] / merging_grid.Nz)
    merging_grid.Δv_inv = 1.0 ./ merging_grid.Δv
end

"""
    compute_grid_index(merging_grid, v)

Compute index of cell on the merging grid in which a velocity is located (the last 8 indices
correspond to the velocity octants outside the grid).

# Positional arguments:
* `merging_grid`: the grid merging (`GridN2Merge`) instance defining the velocity space grid
* `v`: the velocity vector

# Returns:
Index of cell on the merging grid.
"""
function compute_grid_index(merging_grid, v)
    outside_flag = false
    
    if (v[1] < merging_grid.extent_v_lower[1]) || (v[1] > merging_grid.extent_v_upper[1])
        outside_flag = true
    elseif (v[2] < merging_grid.extent_v_lower[2]) || (v[2] > merging_grid.extent_v_upper[2])
        outside_flag = true
    elseif (v[3] < merging_grid.extent_v_lower[3]) || (v[3] > merging_grid.extent_v_upper[3])
        outside_flag = true
    end
        
    if (!outside_flag)

        i_x = (v[1] - merging_grid.extent_v_lower[1]) * merging_grid.Δv_inv[1]
        i_y = (v[2] - merging_grid.extent_v_lower[2]) * merging_grid.Δv_inv[2]
        i_z = (v[3] - merging_grid.extent_v_lower[3]) * merging_grid.Δv_inv[3]

        index = floor(Int64, i_x) * merging_grid.NyNz + floor(Int64, i_y) * merging_grid.Nz + floor(Int64, i_z) + 1
    else
        index = merging_grid.Ntotal - 7

        if (v[1] > merging_grid.extent_v_mid[1])
            index += 1
        end
        if (v[2] > merging_grid.extent_v_mid[2])
            index += 2
        end
        if (v[3] > merging_grid.extent_v_mid[3])
            index += 4
        end
    end
  
    return index
end

"""
    clear_merging_grid!(merging_grid)

Resets all data for a merging grid instance.

# Positional arguments:
* `merging_grid`: the grid merging (`GridN2Merge`) instance defining the velocity space grid
"""
function clear_merging_grid!(merging_grid)
    for index in 1:merging_grid.Ntotal
        merging_grid.cells[index].w = 0.0
        merging_grid.cells[index].v_mean = SVector{3,Float64}(0.0, 0.0, 0.0)
        merging_grid.cells[index].v_std_sq = SVector{3,Float64}(0.0, 0.0, 0.0)
        merging_grid.cells[index].x_mean = SVector{3,Float64}(0.0, 0.0, 0.0)
        merging_grid.cells[index].x_std_sq = SVector{3,Float64}(0.0, 0.0, 0.0)
        merging_grid.cells[index].np = 0
        merging_grid.cells[index].particle_index1 = 0
        merging_grid.cells[index].particle_index2 = 0
    end
end

"""
    compute_grid!(merging_grid::GridN2Merge, particles, pia, cell, species)

Compute all the required cell properties for a grid-based merge.

# Positional arguments:
* `merging_grid`: the grid merging (`GridN2Merge`) instance defining the velocity space grid
* `particles`: the `ParticleVector` instance of the particles to be merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the cell index
* `species`: the species index
"""
function compute_grid!(merging_grid::GridN2Merge, particles, pia, cell, species)
    clear_merging_grid!(merging_grid)

    for i in pia.indexer[cell,species].start1:pia.indexer[cell,species].end1
        index = compute_grid_index(merging_grid, particles[i].v)

        merging_grid.cells[index].np += 1
        merging_grid.cells[index].w += particles[i].w
        merging_grid.cells[index].v_mean = merging_grid.cells[index].v_mean + particles[i].v * particles[i].w
        merging_grid.cells[index].x_mean = merging_grid.cells[index].x_mean + particles[i].x * particles[i].w

        if (merging_grid.cells[index].np == 1)
            merging_grid.cells[index].particle_index1 = i
        elseif (merging_grid.cells[index].np == 2)
            merging_grid.cells[index].particle_index2 = i
        end
    end

    if pia.indexer[cell,species].start2 > 0
        for i in pia.indexer[cell,species].start2:pia.indexer[cell,species].end2
            index = compute_grid_index(merging_grid, particles[i].v)

            merging_grid.cells[index].np += 1
            merging_grid.cells[index].w += particles[i].w
            merging_grid.cells[index].v_mean = merging_grid.cells[index].v_mean + particles[i].v * particles[i].w
            merging_grid.cells[index].x_mean = merging_grid.cells[index].x_mean + particles[i].x * particles[i].w

            if (merging_grid.cells[index].np == 1)
                merging_grid.cells[index].particle_index1 = i
            elseif (merging_grid.cells[index].np == 2)
                merging_grid.cells[index].particle_index2 = i
            end
        end
    end

    for index in 1:merging_grid.Ntotal
        if (merging_grid.cells[index].w > 0.0)
            merging_grid.cells[index].v_mean = merging_grid.cells[index].v_mean / merging_grid.cells[index].w
            merging_grid.cells[index].x_mean = merging_grid.cells[index].x_mean / merging_grid.cells[index].w
        end
    end

    for i in pia.indexer[cell,species].start1:pia.indexer[cell,species].end1
        index = compute_grid_index(merging_grid, particles[i].v)

        merging_grid.cells[index].v_std_sq = merging_grid.cells[index].v_std_sq + (particles[i].v - merging_grid.cells[index].v_mean).^2 * particles[i].w
        merging_grid.cells[index].x_std_sq = merging_grid.cells[index].x_std_sq + (particles[i].x - merging_grid.cells[index].x_mean).^2 * particles[i].w
    end

    if pia.indexer[cell,species].start2 > 0
        for i in pia.indexer[cell,species].start2:pia.indexer[cell,species].end2
            index = compute_grid_index(merging_grid, particles[i].v)

            merging_grid.cells[index].v_std_sq = merging_grid.cells[index].v_std_sq + (particles[i].v - merging_grid.cells[index].v_mean).^2 * particles[i].w
            merging_grid.cells[index].x_std_sq = merging_grid.cells[index].x_std_sq + (particles[i].x - merging_grid.cells[index].x_mean).^2 * particles[i].w
        end
    end

    for index in 1:merging_grid.Ntotal
        if (merging_grid.cells[index].w > 0.0)
            merging_grid.cells[index].v_std_sq = merging_grid.cells[index].v_std_sq / merging_grid.cells[index].w
            merging_grid.cells[index].x_std_sq = merging_grid.cells[index].x_std_sq / merging_grid.cells[index].w
        end
    end
end

"""
    compute_new_particles!(rng, merging_grid::GridN2Merge, particles, pia, cell, species)

Compute new particles based on the grid cell properties without checking particle locations.
So particles may end up outside of the domain.

# Positional arguments:
* `rng`: the random number generator instance
* `merging_grid`: the grid merging (`GridN2Merge`) instance defining the velocity space grid
* `particles`: the `ParticleVector` instance of the particles to be merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the cell index
* `species`: the species index
"""
function compute_new_particles!(rng, merging_grid::GridN2Merge, particles, pia, cell, species)
    # no limits on particle location, i.e. 0-D

    for index in 1:merging_grid.Ntotal
        if (merging_grid.cells[index].np > 2)
            merging_grid.cells[index].w1 = 0.5 * merging_grid.cells[index].w
            merging_grid.cells[index].w2 = merging_grid.cells[index].w1

            merging_grid.cells[index].v_std_sq = sqrt.(merging_grid.cells[index].v_std_sq)
            merging_grid.cells[index].x_std_sq = sqrt.(merging_grid.cells[index].x_std_sq)
            
            merging_grid.direction_vec = @SVector rand(rng, direction_signs, 3)
            merging_grid.cells[index].v1 = merging_grid.cells[index].v_mean + merging_grid.direction_vec .* merging_grid.cells[index].v_std_sq
            merging_grid.cells[index].v2 = merging_grid.cells[index].v_mean - merging_grid.direction_vec .* merging_grid.cells[index].v_std_sq

            merging_grid.direction_vec = @SVector rand(rng, direction_signs, 3)
            merging_grid.cells[index].x1 = merging_grid.cells[index].x_mean + merging_grid.direction_vec .* merging_grid.cells[index].x_std_sq
            merging_grid.cells[index].x2 = merging_grid.cells[index].x_mean - merging_grid.direction_vec .* merging_grid.cells[index].x_std_sq
        elseif (merging_grid.cells[index].np == 2)
            # get the particle indices we saved and just write data based on them
            i = merging_grid.cells[index].particle_index1
            merging_grid.cells[index].w1 = particles[i].w
            merging_grid.cells[index].v1 = particles[i].v
            merging_grid.cells[index].x1 = particles[i].x

            i = merging_grid.cells[index].particle_index2
            merging_grid.cells[index].w2 = particles[i].w
            merging_grid.cells[index].v2 = particles[i].v
            merging_grid.cells[index].x2 = particles[i].x
        elseif (merging_grid.cells[index].np == 1)
            # get the particle indices we saved and just write data based on them
            i = merging_grid.cells[index].particle_index1
            merging_grid.cells[index].w1 = particles[i].w
            merging_grid.cells[index].v1 = particles[i].v
            merging_grid.cells[index].x1 = particles[i].x
        end
    end

    curr_particle_index = 0
    for index in 1:merging_grid.Ntotal
        if (merging_grid.cells[index].np >= 2)
            i = map_cont_index(pia.indexer[cell,species], curr_particle_index)
            curr_particle_index += 1
            particles[i].w = merging_grid.cells[index].w1
            particles[i].v = merging_grid.cells[index].v1
            particles[i].x = merging_grid.cells[index].x1

            i = map_cont_index(pia.indexer[cell,species], curr_particle_index)
            curr_particle_index += 1
            particles[i].w = merging_grid.cells[index].w2
            particles[i].v = merging_grid.cells[index].v2
            particles[i].x = merging_grid.cells[index].x2
        elseif (merging_grid.cells[index].np == 1)
            i = map_cont_index(pia.indexer[cell,species], curr_particle_index)
            curr_particle_index += 1
            particles[i].w = merging_grid.cells[index].w1
            particles[i].v = merging_grid.cells[index].v1
            particles[i].x = merging_grid.cells[index].x1
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
    merge_grid_based!(rng, merging_grid, particles, pia, cell, species, species_data, phys_props)

Merge particles using a velocity grid-based merging approach. A Cartesian grid in velocity
space is used to group particles together (particles outside of the grid are group by velocity
octant), and in each cell/octant, particles are merged down to 2 particles.
The extent of the grid is based on the temperature for the species in question
in the physical grid cell being considered, as stored in the `phys_props` parameter.

# Positional arguments:
* `rng`: the random number generator instance
* `merging_grid`: the grid merging (`GridN2Merge`) instance defining the velocity space grid
* `particles`: the `ParticleVector` instance of the particles to be merged
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the cell index
* `species`: the species index
* `species_data`: the array of `Species` data
* `phys_props`: the `PhysProps` instance containing the computed temperature
"""
function merge_grid_based!(rng, merging_grid, particles, pia, cell, species, species_data, phys_props)
    # 0-D, no grid, particles in single cell
    compute_velocity_extent!(merging_grid, cell, species, species_data, phys_props)
    compute_grid!(merging_grid, particles, pia, cell, species)
    compute_new_particles!(rng, merging_grid, particles, pia, cell, species)
end