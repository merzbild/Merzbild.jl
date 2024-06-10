using StaticArrays

mutable struct MergingGrid
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
end

function create_merging_grid(Nx, Ny, Nz, extent_multiplier::T) where T <: AbstractArray
    Ntotal = Nx * Ny * Nz + 8
    cells = Vector{GridCell}(undef, Nx * Ny * Nz + 8)

    for i in 1:Ntotal
        cells[i] = GridCell(0, 0.0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0, 0,
                            0.0, 0.0, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0])
    end

    return MergingGrid(Nx, Ny, Nz, Ny*Nz, Ntotal, extent_multiplier, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], cells)
end

function create_merging_grid(N, extent_multiplier::T) where T <: AbstractArray
    return create_merging_grid(N, N, N, extent_multiplier)
end

function create_merging_grid(Nx, Ny, Nz, extent_multiplier::Float64)
    return create_merging_grid(Nx, Ny, Nz, [extent_multiplier, extent_multiplier, extent_multiplier])
end

function create_merging_grid(Nx, Ny, Nz, extent_multiplier_x::Float64, extent_multiplier_y::Float64, extent_multiplier_z::Float64)
    return create_merging_grid(Nx, Ny, Nz, [extent_multiplier_x, extent_multiplier_y, extent_multiplier_z])
end

function create_merging_grid(N, extent_multiplier::Float64)
    return create_merging_grid(N, N, N, extent_multiplier)
end

function compute_velocity_extent!(cell, species, merging_grid, phys_props, species_data)
    dv = merging_grid.extent_multiplier .* sqrt.(2 * phys_props.T[cell, species] * k_B / species_data[species].mass)
    merging_grid.extent_v_lower = phys_props.v[:, cell, species] .- dv
    merging_grid.extent_v_upper = phys_props.v[:, cell, species] .+ dv
    merging_grid.extent_v_mid = phys_props.v[:, cell, species]
    merging_grid.Δv = SVector{3}(2 * dv[1] / merging_grid.Nx, 2 * dv[2] / merging_grid.Ny, 2 * dv[3] / merging_grid.Nz)
    merging_grid.Δv_inv = 1.0 ./ merging_grid.Δv
end

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

function compute_grid!(cell, species, merging_grid::MergingGrid, particles, particle_indexer_array)
    clear_merging_grid!(merging_grid)

    for i in particle_indexer_array.indexer[cell,species].start1:particle_indexer_array.indexer[cell,species].end1
        index = compute_grid_index(merging_grid, particles[species][i].v)

        merging_grid.cells[index].np += 1
        merging_grid.cells[index].w += particles[species][i].w
        merging_grid.cells[index].v_mean = merging_grid.cells[index].v_mean + particles[species][i].v * particles[species][i].w
        merging_grid.cells[index].x_mean = merging_grid.cells[index].x_mean + particles[species][i].x * particles[species][i].w

        if (merging_grid.cells[index].np == 1)
            merging_grid.cells[index].particle_index1 = i
        elseif (merging_grid.cells[index].np == 2)
            merging_grid.cells[index].particle_index2 = i
        end
    end

    if particle_indexer_array.indexer[cell,species].start2 > 0
        for i in particle_indexer_array.indexer[cell,species].start2:particle_indexer_array.indexer[cell,species].end2
            index = compute_grid_index(merging_grid, particles[species][i].v)

            merging_grid.cells[index].np += 1
            merging_grid.cells[index].w += particles[species][i].w
            merging_grid.cells[index].v_mean = merging_grid.cells[index].v_mean + particles[species][i].v * particles[species][i].w
            merging_grid.cells[index].x_mean = merging_grid.cells[index].x_mean + particles[species][i].x * particles[species][i].w

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

    for i in particle_indexer_array.indexer[cell,species].start1:particle_indexer_array.indexer[cell,species].end1
        index = compute_grid_index(merging_grid, particles[species][i].v)

        merging_grid.cells[index].v_std_sq = merging_grid.cells[index].v_std_sq + (particles[species][i].v - merging_grid.cells[index].v_mean).^2 * particles[species][i].w
        merging_grid.cells[index].x_std_sq = merging_grid.cells[index].x_std_sq + (particles[species][i].x - merging_grid.cells[index].x_mean).^2 * particles[species][i].w
    end

    if particle_indexer_array.indexer[cell,species].start2 > 0
        for i in particle_indexer_array.indexer[cell,species].start2:particle_indexer_array.indexer[cell,species].end2
            index = compute_grid_index(merging_grid, particles[species][i].v)

            merging_grid.cells[index].v_std_sq = merging_grid.cells[index].v_std_sq + (particles[species][i].v - merging_grid.cells[index].v_mean).^2 * particles[species][i].w
            merging_grid.cells[index].x_std_sq = merging_grid.cells[index].x_std_sq + (particles[species][i].x - merging_grid.cells[index].x_mean).^2 * particles[species][i].w
        end
    end

    for index in 1:merging_grid.Ntotal
        if (merging_grid.cells[index].w > 0.0)
            merging_grid.cells[index].v_std_sq = merging_grid.cells[index].v_std_sq / merging_grid.cells[index].w
            merging_grid.cells[index].x_std_sq = merging_grid.cells[index].x_std_sq / merging_grid.cells[index].w
        end
    end
end

function compute_new_particles!(cell, species, merging_grid::MergingGrid, particles, particle_indexer_array)
    # no limits on particle location, i.e. 0-D

    for index in 1:merging_grid.Ntotal
        if (merging_grid.cells[index].np > 2)
            merging_grid.cells[index].w1 = 0.5 * merging_grid.cells[index].w
            merging_grid.cells[index].w2 = merging_grid.cells[index].w1

            merging_grid.cells[index].v_std_sq = sqrt.(merging_grid.cells[index].v_std_sq)
            merging_grid.cells[index].x_std_sq = sqrt.(merging_grid.cells[index].x_std_sq)
            
            merging_grid.direction_vec = @SVector rand(direction_signs, 3)
            merging_grid.cells[index].v1 = merging_grid.cells[index].v_mean + merging_grid.direction_vec .* merging_grid.cells[index].v_std_sq
            merging_grid.cells[index].v2 = merging_grid.cells[index].v_mean - merging_grid.direction_vec .* merging_grid.cells[index].v_std_sq

            merging_grid.direction_vec = @SVector rand(direction_signs, 3)
            merging_grid.cells[index].x1 = merging_grid.cells[index].x_mean + merging_grid.direction_vec .* merging_grid.cells[index].x_std_sq
            merging_grid.cells[index].x2 = merging_grid.cells[index].x_mean - merging_grid.direction_vec .* merging_grid.cells[index].x_std_sq
        elseif (merging_grid.cells[index].np == 2)
            # get the particle indices we saved and just write data based on them
            i = merging_grid.cells[index].particle_index1
            merging_grid.cells[index].w1 = particles[species][i].w
            merging_grid.cells[index].v1 = particles[species][i].v
            merging_grid.cells[index].x1 = particles[species][i].x

            i = merging_grid.cells[index].particle_index2
            merging_grid.cells[index].w2 = particles[species][i].w
            merging_grid.cells[index].v2 = particles[species][i].v
            merging_grid.cells[index].x2 = particles[species][i].x
        elseif (merging_grid.cells[index].np == 1)
            # get the particle indices we saved and just write data based on them
            i = merging_grid.cells[index].particle_index1
            merging_grid.cells[index].w1 = particles[species][i].w
            merging_grid.cells[index].v1 = particles[species][i].v
            merging_grid.cells[index].x1 = particles[species][i].x
        end
    end

    curr_particle_index = 0
    for index in 1:merging_grid.Ntotal
        if (merging_grid.cells[index].np >= 2)
            i = map_cont_index(particle_indexer_array.indexer[cell,species], curr_particle_index)
            curr_particle_index += 1
            particles[species][i].w = merging_grid.cells[index].w1
            particles[species][i].v = merging_grid.cells[index].v1
            particles[species][i].x = merging_grid.cells[index].x1

            i = map_cont_index(particle_indexer_array.indexer[cell,species], curr_particle_index)
            curr_particle_index += 1
            particles[species][i].w = merging_grid.cells[index].w2
            particles[species][i].v = merging_grid.cells[index].v2
            particles[species][i].x = merging_grid.cells[index].x2
        elseif (merging_grid.cells[index].np == 1)
            i = map_cont_index(particle_indexer_array.indexer[cell,species], curr_particle_index)
            curr_particle_index += 1
            particles[species][i].w = merging_grid.cells[index].w1
            particles[species][i].v = merging_grid.cells[index].v1
            particles[species][i].x = merging_grid.cells[index].x1
        end
    end

    update_particle_indexer_new_lower_count(species, cell, particle_indexer_array, curr_particle_index)
end

function merge_grid_based!(cell, species, merging_grid, phys_props, species_data, particles, particle_indexer_array)
    # 0-D, no grid, particles in single cell
    compute_velocity_extent!(cell, species, merging_grid, phys_props, species_data)
    compute_grid!(cell, species, merging_grid, particles, particle_indexer_array)
    compute_new_particles!(cell, species, merging_grid, particles, particle_indexer_array)
end