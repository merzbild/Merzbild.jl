module Merzbild

include("constants.jl")
include("particles.jl")
include("distributions_and_sampling.jl")
include("physical_props.jl")
include("collisions.jl")
include("io.jl")
include("merging/merging.jl")

const OCTREE_DEFAULT_BUFFER_SIZE::Int32 = 8192
const DELTA_PARTICLES::Int32 = 256

export sample_maxwellian_on_grid!, sample_on_grid!, bkw, maxwellian
export load_species_list, Particle, sample_particles_equal_weight!, create_particle_indexer, create_particle_indexer_array
export Species, Interaction
export create_props, compute_props!, compute_props_sorted_without_moments!
export ParticleIndexer, PhysProps, CollisionFactors, CollisionData
export create_netcdf_phys_props
export write_netcdf_phys_props
export load_interaction_data
export create_collision_factors, create_collision_data
export estimate_sigma_g_w_max, estimate_sigma_g_w_max!
export ntc!
export k_B
export create_merging_grid, merge_grid_based!
export OctreeBinMidSplit, OctreeBinMeanSplit, OctreeBinMedianSplit
export OctreeInitBinMinMaxVel, OctreeInitBinMinMaxVelSym, OctreeInitBinC
export create_merging_octree, merge_octree_N2_based!

end # module merzbild
