using StaticArrays
using TOML

"""
    Particle

A structure to store information about a single particle.

# Fields
* `w`: the computational weight of the particle
* `v`: the 3-dimensional velocity vector of the particle
* `x`: the 3-dimensional position of the particle
"""
mutable struct Particle
    w::Float64
    v::SVector{3,Float64}
    x::SVector{3,Float64}
end

"""
    Species

A structure to store information about a chemical species.

# Fields
* `name`: the name of the species
* `mass`: the molecular mass of the species
* `charge`: the charge of the species in terms of elementary charge (i.e. 1, -1, etc.)
* `charge_div_mass`: the charge of the species divided by its mass, C/kg
"""
struct Species
    name::String
    mass::Float64
    charge::Float64  # in terms of elemenary charge
    charge_div_mass::Float64  # q/m, C/kg
end

"""
    ParticleIndexer

The structure used to index particles of a given species in a given cell.
It is assumed that the particle indices are contiguous; they may be split
across two groups of contiguous indices.

# Fields
* `n_local`: the number of particles in the cell
* `start1`: the first index in the first group of particle indices
* `end1`: the last index in the first group of particle indices
* `n_group1`: the number of particles in the first group (`n_group1 = end1 - start1 + 1`)
* `start2`: the first index in the second group of particle indices, if no particles are present
    in the group, it should be <= 0
* `end2`: the last index in the second group of particle indices
* `n_group2`: the number of particles in the second group (`n_group2 = end2 - start2 + 1`,
    unless no particles are present in the second group)
"""
mutable struct ParticleIndexer
    n_local::Int64

    start1::Int64  
    end1::Int64
    n_group1::Int64  # = end1 - start1 + 1

    start2::Int64
    end2::Int64
    n_group2::Int64  # = end2 - start2 + 1
end

"""
    ParticleIndexer(n_particles)

Create a ParticleIndexer given a number of particles. All
particles are in group 1, with indices starting from 1.

# Positional arguments
* `n_particles`: the number of particles
"""
ParticleIndexer(n_particles) = ParticleIndexer(n_particles, 1, n_particles, n_particles, 0, 0, 0)

"""
    ParticleIndexer()

Create an empty ParticleIndexer. `end1` and `end2` are set to -1, the other fields are set to 0.
"""
ParticleIndexer() = ParticleIndexer(0, 0, -1, 0, 0, -1, 0)

"""
    ParticleIndexerArray

A structure to store an array of `ParticleIndexer` instances
for each species in each cell. It also stores information about the whether the ParticleIndexer instances for each species are "contiguous".
A set of ParticleIndexer instances (for a specific species) are "contiguous" if the following conditions are fulfilled
(assuming all cells contain particles and in all cells have both `n_group1>0` and `n_group2>0`, for a more
detailed explanation, one is referred to [the documentation on contiguous indexing](@ref "Contiguous indexing")):
* `pia.indexer[cell,species].end1 + 1 == pia.indexer[cell+1,species].start1`, `cell = 1,...,n_cells - 1`
* `pia.indexer[n_cells,species].end1 + 1 == pia.indexer[1,species].start2`
* `pia.indexer[cell,species].end2 + 1 == pia.indexer[cell+1,species].start2`, `cell = 1,...,n_cells - 1`

# Fields
* `indexer`: the array of size `(n_cells, n_species)` (number of grid cells * number of species in the simulation)
  storing the `ParticleIndexer` instances
* `n_total`: vector of length `n_species` storing the total number of particles of each species
* `contiguous`: vector of length `n_species` storing a boolean flag whether the ParticleIndexer instances for a species are "contiguous"
"""
mutable struct ParticleIndexerArray
    indexer::Array{ParticleIndexer,2}  # cells x species
    n_total::Vector{Int64}  # per-species
    contiguous::Vector{Bool}  # per-species

    @doc """
        ParticleIndexerArray(indexer_arr::Array{ParticleIndexer,2}, n_total)

    Creates a `ParticleIndexerArray` from a 2-D array of `ParticleIndexer` instances.

    # Positional arguments
    * `indexer_arr`: the 2-D array of `ParticleIndexer` instances
    * `n_total`: the vector of the total number of particles of each species
    """
    function ParticleIndexerArray(indexer_arr::Array{ParticleIndexer,2}, n_total)
        return new(indexer_arr, n_total, [true for i in 1:length(n_total)])
    end

    @doc """
        ParticleIndexerArray(n_cells::Integer, n_species::Integer)

    Create an empty `ParticleIndexerArray` given the number of cells and species

    # Positional arguments
    * `n_cells`: the number of grid cells
    * `n_species`: the number of species
    """
    function ParticleIndexerArray(n_cells::Integer, n_species::Integer)  # most generic version
        pia_indexer = Array{ParticleIndexer, 2}(undef, (n_cells, n_species))

        for j in 1:n_species
            for i in 1:n_cells
                @inbounds pia_indexer[i, j] = ParticleIndexer()
            end
        end
        return new(pia_indexer, [0 for i in 1:n_species], [true for i in 1:n_species])
    end
end

"""
    ParticleIndexerArray(n_particles::Integer)

Create a single-species/single-cell `ParticleIndexerArray`.

# Positional arguments
* `n_particles`: the number of particles (integer number)
"""
ParticleIndexerArray(n_particles::Integer) = ParticleIndexerArray(hcat(ParticleIndexer(n_particles)), [n_particles])


"""
    ParticleIndexerArray(n_particles::T) where T<:AbstractVector

Create a multi-species/single-cell `ParticleIndexerArray`.

# Positional arguments
* `n_particles`: the number of particles of each species (vector-like)
"""
ParticleIndexerArray(n_particles::T) where T<:AbstractVector = ParticleIndexerArray(reshape([ParticleIndexer(np) for np in n_particles], 1, :),
                                                                                    copy(n_particles))
"""
    ParticleIndexerArray(grid, species_data::Array{Species}) where T<:AbstractVector

Create an empty multi-species/multi-cell `ParticleIndexerArray`.

# Positional arguments
* `grid`: the simulation grid
* `species_data`: array of `Species` data for all of the species in the simulation
"""
ParticleIndexerArray(grid, species_data::Array{Species}) = ParticleIndexerArray(grid.n_cells, length(species_data))

"""
    ParticleVector

The structure used to store particles, sort and keep track of particle indices, and keep track of unused particles.
The lengths of the `particles`, `index`, `cell`, and `buffer` vectors are all the same (and stay the same during
resizing of a `ParticleVector` instance). Only the first `nbuffer` elements of the `buffer` vector store
indices of the actually unused particles.

Accessing `ParticleVector[i]` will return a `Particle`, with the actual particle returned being
`ParticleVector.particles[ParticleVector.index[i]]`.

# Fields
* `particles`: the vector of particles of a single species
* `index`: the vector of indices of the particles (these are sorted in grid sorting, not the particles themselves)
* `cell`: the vector storing information in which cell a particle is located (used in grid sorting routines)
* `buffer`: a last-in-first-out (LIFO) queue keeping track of pre-allocated but unused particles
* `nbuffer`: the number of elements in the buffer
"""
mutable struct ParticleVector
    particles::Vector{Particle}
    index::Vector{Int64}
    cell::Vector{Int64}
    buffer::Vector{Int64}  # LIFO queue to keep track which particles we can write to
    nbuffer::Int64  # length of buffer
end

"""
    ParticleVector(np)

Create an empty `ParticleVector` instance of length `np` (all vectors will have length `np`).

# Positional arguments
* `np`: the length of the `ParticleVector` instance to create
"""
ParticleVector(np) = ParticleVector(Vector{Particle}(undef, np), Vector{Int64}(1:np), zeros(Int64, np),
                                    Vector{Int64}(np:-1:1), np)

"""
    Base.getindex(pv::ParticleVector, i)

Returns the underlying particle in a `ParticleVector` instance with index `i`.

Is usually called as `ParticleVector[i]`.

# Positional arguments
* `pv`: `ParticleVector` instance
* `i`: the index of the particle to be selected
"""
@inline function Base.getindex(pv::ParticleVector, i)
    @inbounds return pv.particles[pv.index[i]]
end

"""
    Base.setindex!(pv::ParticleVector, p::Particle, i::Integer)

Set the underlying particle in a `ParticleVector` instance with index `i` to a new particle.

Is usually called as `ParticleVector[i] = p`.

# Positional arguments
* `pv`: `ParticleVector` instance
* `p`: the `Particle` instance to write
* `i`: the index of the particle to be written to
"""
@inline function Base.setindex!(pv::ParticleVector, p::Particle, i::Integer)
    @inbounds pv.particles[pv.index[i]] = p
end

"""
    Base.length(pv::ParticleVector)

Returns the length of a `ParticleVector` instance.

Is usually called as `length(ParticleVector)`.

# Positional arguments
* `pv`: `ParticleVector` instance
"""
@inline function Base.length(pv::ParticleVector)
    return length(pv.particles)
end

"""
    Base.resize!(pv::ParticleVector, n::Integer)

Resize a `ParticleVector` instance.

Is usually called as `resize!(ParticleVector, n)`.

# Positional arguments
* `pv`: `ParticleVector` instance
* `n`: the new length of the `ParticleVector` instance (i.e. the length of all the vector fields of the instance)
"""
function Base.resize!(pv::ParticleVector, n::Integer)
    old_len = length(pv.particles)
    resize!(pv.particles, n)
    resize!(pv.index, n)
    resize!(pv.cell, n)
    resize!(pv.buffer, n)
    n_diff = n - old_len

    # fill with new indices
    @inbounds pv.index[old_len + 1:n] = old_len + 1:n

    # [1,2,3,(4,5,6,7)] -> [4,5,6,7,(1,2,3)]
    # move older particles to tail of buffer so that they get used up first
    @inbounds pv.buffer[n_diff + 1:n] .= pv.buffer[1:old_len]

    # write new indices to buffer in reverse order (first particles in the array to be used up first)
    @inbounds pv.buffer[1:n_diff] = n:-1:old_len+1
    @inbounds pv.nbuffer += n_diff
end

"""
    update_particle_buffer_new_particle!(pv::ParticleVector, position)

Update the buffer in a `ParticleVector` instance when a new particle is created. This writes the index of the new particle
(the last index stored in the active part of the buffer) to the `index` vector at position `position`, and reduces
the length of the active part of the buffer by 1.

# Positional arguments
* `pv`: `ParticleVector` instance
* `position`: the position in the `index` vector to which to write the index of the new particle
"""
@inline function update_particle_buffer_new_particle!(pv::ParticleVector, position)
    # position is where we will be writing to
    @inbounds pv.index[position] = pv.buffer[pv.nbuffer]
    pv.nbuffer -= 1
end

"""
    update_particle_buffer_new_particle!(pv::ParticleVector, pia, species)

Update the buffer in a `ParticleVector` instance when a new particle is created at the end of the particle array, and reduces
the length of the active part of the buffer by 1.
This assumes that the `pia` structure has already an updated particle count (that accounts for the newly created particle),
as the index of the new particle taken from the buffer is written to `pv.index.[pia.n_total[species]]`.

# Positional arguments
* `pv`: `ParticleVector` instance
* `pia`: the `ParticleIndexerArray` instance
* `species`: the index of the species of which a new particle is created
"""
@inline function update_particle_buffer_new_particle!(pv::ParticleVector, pia, species)
    update_particle_buffer_new_particle!(pv, pia.n_total[species])
end

"""
    update_particle_buffer_new_particle!(pv::Vector{Particle}, pia, species)

Dummy function in case `Vector{Particle}` is used and not a `ParticleVector`, just to make the simplest 0-D examples work.
"""
@inline function update_particle_buffer_new_particle!(pv::Vector{Particle}, pia, species)
    # dummy function, might remove it at some point
    nothing
end

"""
    update_particle_buffer_new_particle!(pv::Vector{Particle}, position)

Dummy function in case `Vector{Particle}` is used and not a `ParticleVector`, just to make the simplest 0-D examples work.
"""
@inline function update_particle_buffer_new_particle!(pv::Vector{Particle}, position)
    # dummy function, might remove it at some point
    nothing
end

"""
    map_cont_index(particle_indexer, i)

Maps a continuous index in the range `[0, n_local-1]` to an index in the particle array given
a `ParticleIndexer` instance describing how particle indices are split across 2 groups.

# Positional arguments
* `particle_indexer`: the `ParticleIndexer` instance
* `i`: the index to map
"""
@inline function map_cont_index(particle_indexer, i)
    return i < particle_indexer.n_group1 ? i + particle_indexer.start1 : (i - particle_indexer.n_group1) + particle_indexer.start2
end

"""
    map_cont_index(pia, cell, species, i)

Maps a continuous index in the range `[0, n_local-1]` to an index in the particle array
for particles of a specific species in a specific cell.

# Positional arguments
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the cell in which the particles are located
* `species`: the index of the particles' species
* `i`: the index to map
"""
@inline function map_cont_index(pia, cell, species, i)
    @inbounds return map_cont_index(pia.indexer[cell, species], i)
end

"""
    update_particle_indexer_new_lower_count!(pia, cell, species, new_lower_count)

Update a `ParticleIndexerArray` instance when the particle count of a given species in a given cell is reduced.

# Positional arguments
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the cell in which the particles are located
* `species`: the index of the particles' species
* `new_lower_count`: the new number of particles of the given species in the given cell
"""
function update_particle_indexer_new_lower_count!(pia, cell, species, new_lower_count)
    @inbounds diff = pia.indexer[cell, species].n_local - new_lower_count
    @inbounds pia.indexer[cell, species].n_local = new_lower_count

    @inbounds pia.n_total[species] -= diff

    @inbounds if (new_lower_count > pia.indexer[cell, species].n_group1)
        @inbounds pia.indexer[cell, species].end2 -= diff
        @inbounds pia.indexer[cell, species].n_group2 -= diff
    else
        @inbounds diff -= pia.indexer[cell, species].n_group2
        @inbounds pia.indexer[cell, species].start2 = 0
        @inbounds pia.indexer[cell, species].end2 = 0
        @inbounds pia.indexer[cell, species].n_group2 = 0

        @inbounds pia.indexer[cell, species].end1 -= diff
        @inbounds pia.indexer[cell, species].n_group1 -= diff
    end
end

"""
    update_particle_indexer_new_particle!(pia, cell, species)

Update a `ParticleIndexerArray` instance when a particle of a given species in a given cell is created.
This places the particle index in the 2-nd group of particle indices in the `ParticleIndexer` instance.

# Positional arguments
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the cell in which the particle is created
* `species`: the index of the species of which the particle is created
"""
@inline function update_particle_indexer_new_particle!(pia, cell, species)
    @inbounds pia.n_total[species] += 1
    @inbounds pia.indexer[cell, species].n_local += 1
    @inbounds pia.indexer[cell, species].n_group2 += 1

    @inbounds pia.indexer[cell, species].start2 = pia.indexer[cell, species].start2 > 0 ? pia.indexer[cell, species].start2 : pia.n_total[species]
    @inbounds pia.indexer[cell, species].end2 = pia.n_total[species]
end

"""
    delete_particle!(pv::ParticleVector, pia, cell, species, i)

Delete particle with index i of species `species` in cell `cell`
and update the particle indexers and buffers accordingly. This changes the ordering of the non-deleted particles in the cell.

# Positional arguments
* `pv`: `ParticleVector` instance
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the cell in which the particle is deleted
* `species`: the index of the species of which the particle is deleted
* `i`: the index of the particle to delete
"""
@inline function delete_particle!(pv::ParticleVector, pia, cell, species, i)
    # check in which group we are in
    if (pia.indexer[cell, species].n_group2 > 0) && (i >= pia.indexer[cell, species].start2) && (i <= pia.indexer[cell, species].end2)
        @inbounds last_index_group2 = pv.index[pia.indexer[cell, species].end2]
        @inbounds pv.index[pia.indexer[cell, species].end2] = pv.index[i]
        @inbounds pv.index[i] = last_index_group2
        delete_particle_end_group2!(pv, pia, cell, species)
    else
        @inbounds last_index_group1 = pv.index[pia.indexer[cell, species].end1]
        @inbounds pv.index[pia.indexer[cell, species].end1] = pv.index[i]
        @inbounds pv.index[i] = last_index_group1
        delete_particle_end_group1!(pv, pia, cell, species)
    end
end

"""
    delete_particle_end!(pv::ParticleVector, pia, cell, species)

Delete the last particle of species `species` in cell `cell`: if particles are present in the 2nd
group of the indices stored in the `ParticleIndexer` instance, it will delete the last particle in that group;
otherwise it will delete the last particle in the 1st group of particles pointed to by the `ParticleIndexer` instance.
If no particles are present in the cell, the function does nothing. This does not set the value of the `contiguous` field of
`pia` to `false` even if `pia` becomes discontinuous; this has to be done outside of this function.

# Positional arguments
* `pv`: `ParticleVector` instance
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the cell in which the particle is deleted
* `species`: the index of the species of which the particle is deleted
"""
@inline function delete_particle_end!(pv::ParticleVector, pia, cell, species)
    if pia.indexer[cell, species].n_group2 > 0
        delete_particle_end_group2!(pv, pia, cell, species)
    else
        delete_particle_end_group1!(pv, pia, cell, species)
    end
end

"""
    delete_particle_end_group1!(pv::ParticleVector, pia, cell, species)

Delete particle with index `pia.indexer[cell, species].end1`` of species `species` in cell `cell`
and update the particle indexers and buffers accordingly (i.e. delete the last particle in the 1st group of particles
of a given species in a given cell). This also sets the weight of the deleted particle to 0.
If no particles are present in the 1st group of particles, the function does nothing. This does not set the value of the `contiguous` field of
`pia` to `false` even if `pia` becomes discontinuous; this has to be done outside of this function.

# Positional arguments
* `pv`: `ParticleVector` instance
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the cell in which the particle is deleted
* `species`: the index of the species of which the particle is deleted
"""
@inline function delete_particle_end_group1!(pv::ParticleVector, pia, cell, species)
    @inbounds if pia.indexer[cell, species].n_group1 == 0
        return
    end

    @inbounds index_of_deleted = pia.indexer[cell, species].end1

    # set weight to 0
    @inbounds pv[index_of_deleted].w = 0.0

    @inbounds pia.indexer[cell, species].n_local -= 1
    @inbounds pia.indexer[cell, species].end1 -= 1
    @inbounds pia.indexer[cell, species].n_group1 -= 1
    @inbounds pia.n_total[species] -= 1

    # deleted last particle from group1
    @inbounds if pia.indexer[cell, species].end1 < pia.indexer[cell, species].start1
        @inbounds pia.indexer[cell, species].start1 = 0
        @inbounds pia.indexer[cell, species].end1 = -1
    end

    # add the deleted particle to the buffer
    pv.nbuffer += 1
    @inbounds pv.buffer[pv.nbuffer] = pv.index[index_of_deleted]
end

"""
    delete_particle_end_group2!(pv::ParticleVector, pia, cell, species)

Delete particle with index `pia.indexer[cell, species].end2`` of species `species` in cell `cell`
and update the particle indexers and buffers accordingly (i.e. delete the last particle in the 2nd group of the particles
of a given species in a given cell). This also sets the weight of the deleted particle to 0.
If no particles are present in the 2nd group of particles, the function does nothing. This does not set the value of the `contiguous` field of
`pia` to `false` even if `pia` becomes discontinuous; this has to be done outside of this function.

# Positional arguments
* `pv`: `ParticleVector` instance
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the cell in which the particle is deleted
* `species`: the index of the species of which the particle is deleted
"""
@inline function delete_particle_end_group2!(pv::ParticleVector, pia, cell, species)
    @inbounds if pia.indexer[cell, species].n_group2 == 0
        return
    end

    @inbounds index_of_deleted = pia.indexer[cell, species].end2

    # set weight to 0
    @inbounds pv[index_of_deleted].w = 0.0

    @inbounds pia.indexer[cell, species].n_local -= 1
    @inbounds pia.indexer[cell, species].end2 -= 1
    @inbounds pia.indexer[cell, species].n_group2 -= 1
    @inbounds pia.n_total[species] -= 1

    # deleted last particle from group2
    @inbounds if pia.indexer[cell, species].end2 < pia.indexer[cell, species].start2
        @inbounds pia.indexer[cell, species].start2 = 0
        @inbounds pia.indexer[cell, species].end2 = 0
    end

    # add the deleted particle to the buffer
    @inbounds pv.nbuffer += 1
    @inbounds pv.buffer[pv.nbuffer] = pv.index[index_of_deleted]
end

"""
    load_species_data(species_filename, species_names)

Load a vector of species data (mass, charge, etc.) from a TOML file.

# Positional arguments
* `species_filename`: the path to the TOML file containing the data
* `species_names`: a list of the names of the species for which to load the data

# Returns
Vector of `Species` filled with data loaded from the file. 
"""
function load_species_data(species_filename, species_names)
    species_data = TOML.parsefile(species_filename)

    res = Vector{Species}()

    for species_name in species_names
        push!(res, Species(species_name, species_data[species_name]["mass"], species_data[species_name]["charge"],
                           q_e * species_data[species_name]["charge"] / species_data[species_name]["mass"]))
    end

    return res
end


"""
    load_species_data(species_filename, species_name::String)

Load a vector of species data (mass, charge, etc.) from a TOML file for a single species.

# Positional arguments
* `species_filename`: the path to the TOML file containing the data
* `species_name`: the name of the species for which to load the data

# Returns
Vector of `Species` filled with data loaded from the file.
"""
function load_species_data(species_filename, species_name::String)
    return load_species_data(species_filename, [species_name])
end

"""
    squash_pia!(pia, species)

Restore the continuity of indices in a `ParticleVector and associated
`ParticleIndexerArray` instance for a specific species.
If for this species the instance has `contiguous == true`, nothing will be done.

# Positional arguments
* `pv`: the `ParticleVector`
* `pia`: the `ParticleIndexerArray` instance
* `species`: the index of the species for which to restore continuity of indices
"""
function squash_pia!(pv, pia, species)
    if pia.contiguous[species]
        return
    else
        n_cells = size(pia.indexer)[1]
        if n_cells == 1
            @inbounds if pia.indexer[1, species].n_group2 > 0
                @inbounds e1 = pia.indexer[1, species].end1 > 0 ? pia.indexer[1, species].end1 : 0
                @inbounds offset = pia.indexer[1, species].start2 - (e1 + 1)

                if offset > 0
                    @inbounds pia.indexer[1, species].start2 -= offset
                    @inbounds pia.indexer[1, species].end2 -= offset

                    @inbounds s2 = pia.indexer[1, species].start2
                    @inbounds e2 = pia.indexer[1, species].end2
                    for j in s2:e2
                        @inbounds pv.index[j] = pv.index[j+offset]
                    end
                end
            end
        else
            last_end = pia.indexer[1, species].end1 > 0 ? pia.indexer[1, species].end1 : 0
            for i in 1:n_cells-1
                @inbounds offset = pia.indexer[i+1, species].start1 - (last_end + 1)
                if offset > 0
                    @inbounds pia.indexer[i+1, species].start1 -= offset
                    @inbounds pia.indexer[i+1, species].end1 -= offset

                    @inbounds s1 = pia.indexer[i+1, species].start1
                    @inbounds e1 = pia.indexer[i+1, species].end1
                    for j in s1:e1
                        @inbounds pv.index[j] = pv.index[j+offset]
                    end
                end
                @inbounds last_end = pia.indexer[i+1, species].end1 > 0 ? pia.indexer[i+1, species].end1 : last_end
            end

            for i in 1:n_cells
                @inbounds if pia.indexer[i, species].n_group2 > 0
                    offset = pia.indexer[i, species].start2 - (last_end + 1)
                    if offset > 0
                        @inbounds pia.indexer[i, species].start2 -= offset
                        @inbounds pia.indexer[i, species].end2 -= offset
                        
                        @inbounds s2 = pia.indexer[i, species].start2
                        @inbounds e2 = pia.indexer[i, species].end2
                        @inbounds for j in s2:e2
                            @inbounds pv.index[j] = pv.index[j+offset]
                        end
                    end
                    @inbounds last_end = pia.indexer[i, species].end2 > 0 ? pia.indexer[i, species].end2 : last_end
                end
            end
        end
        @inbounds pia.contiguous[species] = true
    end
end

"""
    squash_pia!(pia)

Restore the continuity of indices in a list of `ParticleVector`s and the associated
`ParticleIndexerArray` instance for all species.
If for a specific species the instance has `contiguous == true`, nothing will be done.

# Positional arguments
* `particles`: the list of `ParticleVector`s for all species in the flow
* `pia`: the `ParticleIndexerArray` instance
"""
function squash_pia!(particles, pia)
    for species in 1:length(pia.contiguous)
        @inbounds if !pia.contiguous[species]
            @inbounds squash_pia!(particles[species], pia, species)
        end
    end
end

"""
    update_buffer_index_new_particle!(pv, pia, cell, species)

Update a `ParticleIndexerArray` and the buffer in a `ParticleVector` instance
when a particle of a given species in a given cell is created. The particle index is added
to the second group of particles pointed to by the `ParticleIndexer`.
See the documentation of [`update_particle_indexer_new_particle!`](@ref update_particle_indexer_new_particle!)
and [`update_particle_buffer_new_particle!`](@ref update_particle_buffer_new_particle!)

# Positional arguments
* `pv`: `ParticleVector` instance
* `pia`: the `ParticleIndexerArray` instance
* `cell`: the index of the cell in which the particle is created
* `species`: the index of the species of which the particle is created
"""
@inline function update_buffer_index_new_particle!(pv, pia, cell, species)
    update_particle_indexer_new_particle!(pia, cell, species)
    update_particle_buffer_new_particle!(pv, pia, species)
end

"""
    add_particle!(pv, position, w, v, x)

Create a new particle in a `ParticleVector` instance at position `position`.
The `ParticleIndexer`/`ParticleIndexerArray` instances should be updated
independently. See [`update_particle_buffer_new_particle!`](@ref) for more information
regarding how the buffer of the `ParticleVector` is updated. This should not be used
to update an existing particle. Particles at positions before `position` should exist
in the `ParticleVector` array.

# Positional arguments
* `pv`: `ParticleVector` instance
* `position`: the position in the `index` vector to which to write the index of the new particle
* `w`: the computational weight of the particle to create
* `v`: the velocity of the particle to create
"""
@inline function add_particle!(pv, position, w, v, x)
    update_particle_buffer_new_particle!(pv, position)
    @inbounds pv[position] = Particle(w, v, x)
end


"""
    pretty_print_pia(pia)

Display a `ParticleIndexerArray` instance by showing the starting/ending indices of the groups over all cells
for a specific species.

# Positional arguments
* `pia`: the `ParticleIndexerArray` instance
* `species`: the index of the species for which the indices are displayed
"""
function pretty_print_pia(pia, species)
    n_cells = size(pia.indexer)[1]
    println("Total: $(pia.n_total[species])")
    for cell in 1:n_cells
        if pia.indexer[cell,species].n_group1 > 0
            println("Cell $cell: [$(pia.indexer[cell,species].start1), $(pia.indexer[cell,species].end1)]")
        end
    end
    for cell in 1:n_cells
        if pia.indexer[cell,species].n_group2 > 0
            println("Cell $cell: [$(pia.indexer[cell,species].start2), $(pia.indexer[cell,species].end2)]")
        end
    end
end