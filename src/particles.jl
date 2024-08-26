using StaticArrays
using TOML

mutable struct Particle
    w::Float64
    v::SVector{3,Float64}
    x::SVector{3,Float64}
end

struct Species
    name::String
    mass::Float64
    charge::Float64  # in terms of elemenary charge
    charge_div_mass::Float64  # q/m, C/kg
end

# species x cells
mutable struct ParticleIndexer
    n_local::Int64

    start1::Int64  
    end1::Int64
    n_group1::Int64  # = end1 - start1 + 1

    start2::Int64
    end2::Int64
    n_group2::Int64  # = end2 - start2 + 1
end

ParticleIndexer(n_particles) = ParticleIndexer(n_particles, 1, n_particles, n_particles, 0, 0, 0)

ParticleIndexer() = ParticleIndexer(0, 0, -1, 0, 0, -1, 0)

mutable struct ParticleIndexerArray
    indexer::Array{ParticleIndexer,2}  # cells x species
    n_total::Vector{Int64}  # per-species

    function ParticleIndexerArray(indexer_arr, n_total)
        return new(indexer_arr, n_total)
    end

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

ParticleIndexerArray(n_particles::Int64) = ParticleIndexerArray(hcat(ParticleIndexer(n_particles)), [n_particles])

ParticleIndexerArray(n_particles::T) where T<:AbstractVector = ParticleIndexerArray(reshape([ParticleIndexer(np) for np in n_particles], 1, :),
                                                                                    copy(n_particles))

ParticleIndexerArray(grid, species_data::Array{Species}) = ParticleIndexerArray(grid.n_cells, length(species_data))

mutable struct ParticleVector
    particles::Vector{Particle}
    index::Vector{Int64}
    cell::Vector{Int64}
end

ParticleVector(np) = ParticleVector(Vector{Particle}(undef, np), Vector{Int64}(1:np), zeros(Int64, np))

function Base.getindex(pv::ParticleVector, i)
    return pv.particles[pv.index[i]]
end

function Base.setindex!(pv::ParticleVector, p::Particle, i::Int64)
    pv.particles[pv.index[i]] = p
end

function Base.length(pv::ParticleVector)
    return length(pv.particles)
end

function map_cont_index(particle_indexer, i)
    # map a continuous index in [0, n_local-1]
    # to a index in the particle array given particle_indexer struct describing the split
    return i < particle_indexer.n_group1 ? i + particle_indexer.start1 : (i - particle_indexer.n_group1) + particle_indexer.start2
end

function map_cont_index(pia, cell, species, i)
    # map a continuous index in [0, n_local-1]
    # to a index in the particle array given particle_indexer struct describing the split
    return map_cont_index(pia.indexer[cell, species], i)
end

function update_particle_indexer_new_lower_count(pia, cell, species, new_lower_count)
    # update particle indexer when we reduced the particle count
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

function update_particle_indexer_new_particle(pia, cell, species)
    # update particle indexer when we add a new particle
    pia.n_total[species] += 1
    pia.indexer[cell, species].n_local += 1
    pia.indexer[cell, species].n_group2 += 1

    pia.indexer[cell, species].start2 = pia.indexer[cell, species].start2 > 0 ? pia.indexer[cell, species].start2 : pia.n_total[species]
    pia.indexer[cell, species].end2 = pia.n_total[species]
end

function load_species_data(species_filename, species_names)
    species_data = TOML.parsefile(species_filename)

    res = Vector{Species}()

    for species_name in species_names
        push!(res, Species(species_name, species_data[species_name]["mass"], species_data[species_name]["charge"],
                           q_e * species_data[species_name]["charge"] / species_data[species_name]["mass"]))
    end

    return res
end

function load_species_data(species_filename, species_name::String)
    return load_species_data(species_filename, [species_name])
end