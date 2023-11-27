module Merging

include("constants.jl")
include("particles.jl")
include("distributions_and_sampling.jl")
include("physical_props.jl")
include("collisions.jl")

export load_species_list, Particle, sample_particles_equal_weight!, create_particle_indexer
export create_props, compute_props!, ParticleIndexer

end # module merging
