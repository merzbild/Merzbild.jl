using StaticArrays
using TOML

mutable struct Particle
    w::Float64
    v::SVector{3,Float64}
    x::SVector{3,Float64}
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

mutable struct ParticleIndexerArray
    indexer::Array{ParticleIndexer,2}  # cells x species
    n_total::Vector{Int64}  # per-species
end

struct Species
    name::String
    mass::Float64
end

function map_cont_index(particle_indexer, i)
    # map a continuous index in [0, n_local-1]
    # to a index in the particle array given particle_indexer struct describing the split
    return i < particle_indexer.n_group1 ? i + particle_indexer.start1 : (i - particle_indexer.n_group1) + particle_indexer.start2
end

function map_cont_index(species, cell, pia, i)
    # map a continuous index in [0, n_local-1]
    # to a index in the particle array given particle_indexer struct describing the split
    return map_cont_index(species, cell, pia.indexer[cell, species], i)
end

# function update_particle_indexer_new_lower_count(particle_indexer, new_lower_count)
#     # update particle indexer when we reduced the particle count
#     particle_indexer.n_total += 1
#     particle_indexer.n_group2 += 1

#     particle_indexer.start2 = particle_indexer.start2 > 0 ? particle_indexer.start2 : particle_indexer.n_total 
#     particle_indexer.end2 = particle_indexer.n_total
# end

function update_particle_indexer_new_particle(species, cell, pia)
    # update particle indexer when we add a new particle
    pia.n_total[species] += 1
    pia.indexer[cell, species].n_local += 1
    pia.indexer[cell, species].n_group2 += 1

    pia.indexer[cell, species].start2 = pia.indexer[cell, species].start2 > 0 ? pia.indexer[cell, species].start2 : pia.n_total[species]
    pia.indexer[cell, species].end2 = pia.n_total[species]
end

function load_species_list(species_filename, species_names)
    species_data = TOML.parsefile(species_filename)

    species_list = Vector{Species}()

    for species_name in species_names
        push!(species_list, Species(species_name, species_data[species_name]["mass"]))
    end

    return species_list
end

function load_species_list(species_filename, species_name::String)
    return load_species_list(species_filename, [species_name])
end

function create_particle_indexer(n_particles)
    return ParticleIndexer(n_particles, 1, n_particles, n_particles, 0, 0, 0)
end

function create_particle_indexer_array(n_particles::Int64)  # 1 cell 1 species
    return ParticleIndexerArray(hcat(ParticleIndexer(n_particles, 1, n_particles, n_particles, 0, 0, 0)), [n_particles])
end


function create_particle_indexer_array(n_particles::T) where T<:AbstractVector  # 1 cell multi-species
    return ParticleIndexerArray(reshape([ParticleIndexer(np, 1, np, np, 0, 0, 0) for np in n_particles], 1, :), copy(n_particles))
end