module Particles
using StaticArrays

export Particle, ParticleIndexer, Species

mutable struct Particle
    w::Float64
    v::MVector{3,Float64}
    x::MVector{3,Float64}
end


# species x cells
mutable struct ParticleIndexer
    n_total::Int64

    start1::Int64
    end1::Int64
    n_group1::Int64  # = end1 - start1 + 1

    start2::Int64
    end2::Int64
end


struct Species
    name::String
    mass::Float64
end
end