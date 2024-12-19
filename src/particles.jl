using StaticArrays
using TOML

"""
The Particle struct
"""
mutable struct Particle
    w::Float64
    v::SVector{3,Float64}
    x::SVector{3,Float64}
end

"""
The Species struct
"""
struct Species
    name::String
    mass::Float64
    charge::Float64  # in terms of elemenary charge
    charge_div_mass::Float64  # q/m, C/kg
end

"""
The ParticleIndexer struct
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
Create a ParticleIndexer given a number of particles
"""
ParticleIndexer(n_particles) = ParticleIndexer(n_particles, 1, n_particles, n_particles, 0, 0, 0)

"""
Create an empty ParticleIndexer
"""
ParticleIndexer() = ParticleIndexer(0, 0, -1, 0, 0, -1, 0)

"""
The ParticleIndexerArray struct
"""
mutable struct ParticleIndexerArray
    indexer::Array{ParticleIndexer,2}  # cells x species
    n_total::Vector{Int64}  # per-species

    """
    Create a ParticleIndexerArray from a 2-D array of ParticleIndexers
    """
    function ParticleIndexerArray(indexer_arr, n_total)
        return new(indexer_arr, n_total)
    end

    """
    Create a ParticleIndexerArray given the number of cells and species
    """
    function ParticleIndexerArray(n_cells::Int, n_species::Int)  # most generic version
        pia_indexer = Array{ParticleIndexer, 2}(undef, (n_cells, n_species))

        for j in 1:n_species
            for i in 1:n_cells
                pia_indexer[i, j] = ParticleIndexer()
            end
        end
        return new(pia_indexer, [0 for i in 1:n_species])
    end
end

"""
Create a single-species/single-cell ParticleIndexerArray
"""
ParticleIndexerArray(n_particles::Int) = ParticleIndexerArray(hcat(ParticleIndexer(n_particles)), [n_particles])

"""
Create a multi-species/single-cell ParticleIndexerArray
"""
ParticleIndexerArray(n_particles::T) where T<:AbstractVector = ParticleIndexerArray(reshape([ParticleIndexer(np) for np in n_particles], 1, :),
                                                                                    copy(n_particles))
"""
Create a multi-species/multi-cell ParticleIndexerArrays
"""
ParticleIndexerArray(grid, species_data::Array{Species}) = ParticleIndexerArray(grid.n_cells, length(species_data))

"""
The ParticleVector struct
"""
mutable struct ParticleVector
    particles::Vector{Particle}
    index::Vector{Int64}
    cell::Vector{Int64}
    buffer::Vector{Int64}  # LIFO queue to keep track which particles we can write to
    nbuffer::Int64  # length of buffer
end

"""
Create a ParticleVector of length np
"""
ParticleVector(np) = ParticleVector(Vector{Particle}(undef, np), Vector{Int64}(1:np), zeros(Int64, np),
                                    Vector{Int64}(np:-1:1), np)

"""
Get underlying particle in ParticleVector with index i
"""
function Base.getindex(pv::ParticleVector, i)
    return pv.particles[pv.index[i]]
end

"""
Set underlying particle in ParticleVector with index i to a new particle
"""
function Base.setindex!(pv::ParticleVector, p::Particle, i::Int)
    pv.particles[pv.index[i]] = p
end

"""
Get length of a `ParticleVector`` instance
"""
function Base.length(pv::ParticleVector)
    return length(pv.particles)
end

"""
Resize a `ParticleVector`` instance
"""
function Base.resize!(pv::ParticleVector, n::Int)
    old_len = length(pv.particles)
    resize!(pv.particles, n)
    resize!(pv.index, n)
    resize!(pv.cell, n)
    resize!(pv.buffer, n)
    n_diff = n - old_len

    # [1,2,3,(4,5,6,7)] -> [4,5,6,7,(1,2,3)]
    # move older particles to tail of buffer so that they get used up first
    pv.buffer[n_diff + 1:n] .= pv.buffer[1:old_len]

    # write new indices to buffer
    pv.buffer[1:n_diff] = old_len+1:n
    pv.nbuffer += n_diff
end

"""
Update the buffer of unused `Particle`s in a `ParticleVector` instance
"""
function update_particle_buffer_new_particle(pv::ParticleVector, position)
    # position is where we will be writing to
    pv.index[position] = pv.buffer[pv.nbuffer]
    pv.nbuffer -= 1
end

"""
Update the buffer of unused `Particle`s in a `ParticleVector` instance
"""
function update_particle_buffer_new_particle(pv::ParticleVector, pia, species)
    update_particle_buffer_new_particle(pv, pia.n_total[species])
end

"""
Dummy function in case `Vector{Particle}` is used and not a `ParticleVector``
"""
function update_particle_buffer_new_particle(pv::Vector{Particle}, pia, species)
    # dummy function, might remove it at some point
    nothing
end

"""
Dummy function in case `Vector{Particle}` is used and not a `ParticleVector``
"""
function update_particle_buffer_new_particle(pv::Vector{Particle}, position)
    # dummy function, might remove it at some point
    nothing
end

"""
map a continuous index in `[0, n_local-1]`
to an index in the particle array given particle_indexer struct describing the split
"""
function map_cont_index(particle_indexer, i)
    return i < particle_indexer.n_group1 ? i + particle_indexer.start1 : (i - particle_indexer.n_group1) + particle_indexer.start2
end

"""
map a continuous index in `[0, n_local-1]`
to an index in the particle array given particle_indexer struct describing the split
"""
function map_cont_index(pia, cell, species, i)
    return map_cont_index(pia.indexer[cell, species], i)
end

"""
Update particle indexer when we reduced the particle count
"""
function update_particle_indexer_new_lower_count(pia, cell, species, new_lower_count)
    diff = pia.indexer[cell, species].n_local - new_lower_count
    pia.indexer[cell, species].n_local = new_lower_count

    pia.n_total[species] -= diff

    if (new_lower_count > pia.indexer[cell, species].n_group1)
        pia.indexer[cell, species].end2 -= diff
        pia.indexer[cell, species].n_group2 -= diff
    else
        diff -= pia.indexer[cell, species].n_group2
        pia.indexer[cell, species].start2 = 0
        pia.indexer[cell, species].end2 = 0
        pia.indexer[cell, species].n_group2 = 0

        pia.indexer[cell, species].end1 -= diff
        pia.indexer[cell, species].n_group1 -= diff
    end
end

"""
Update particle indexer when we add a new particle
"""
function update_particle_indexer_new_particle(pia, cell, species)
    pia.n_total[species] += 1
    pia.indexer[cell, species].n_local += 1
    pia.indexer[cell, species].n_group2 += 1

    pia.indexer[cell, species].start2 = pia.indexer[cell, species].start2 > 0 ? pia.indexer[cell, species].start2 : pia.n_total[species]
    pia.indexer[cell, species].end2 = pia.n_total[species]
end

"""
Load a vector of species data
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
Load a vector of species data for a single species
"""
function load_species_data(species_filename, species_name::String)
    return load_species_data(species_filename, [species_name])
end