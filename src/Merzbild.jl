module Merzbild

include("constants.jl")
include("utils.jl")
include("particles.jl")
include("distributions_and_sampling.jl")
include("physical_props.jl")
include("collisions/collisions.jl")
include("io.jl")
include("merging/merging.jl")
include("pic/pic.jl")
include("grids/grids.jl")
include("convection/convection.jl")

const OCTREE_DEFAULT_BUFFER_SIZE::Int32 = 8192
const DELTA_PARTICLES::Int32 = 256

export sample_maxwellian_on_grid!, sample_on_grid!, bkw, maxwellian
export load_species_list, Particle, sample_particles_equal_weight!, create_particle_indexer, create_particle_indexer_array
export Species, Interaction
export create_props, compute_props!, compute_props_sorted_without_moments!
export clear_props!, avg_props!
export ParticleIndexer, PhysProps, CollisionFactors, CollisionData
export create_netcdf_phys_props, create_IO_skip_list
export write_netcdf_phys_props
export close_netcdf
export load_interaction_data, load_interaction_data_with_dummy
export create_collision_factors, create_collision_data
export estimate_sigma_g_w_max, estimate_sigma_g_w_max!
export ntc!
export k_B
export create_merging_grid, merge_grid_based!
export OctreeBinMidSplit, OctreeBinMeanSplit, OctreeBinMedianSplit
export OctreeInitBinMinMaxVel, OctreeInitBinMinMaxVelSym, OctreeInitBinC
export create_merging_octree, merge_octree_N2_based!
export create_nnls_merging, compute_multi_index_moments, merge_nnls_based!
export create_nnls_merging_rate_conserving, merge_nnls_based_rate_preserving!
export load_electron_neutral_interactions, create_computed_crosssections, DataMissingException
export ScatteringIsotropic, ScatteringOkhrimovskyy
export ElectronEnergySplitEqual, ElectronEnergySplitZeroE
export accelerate_constant_field_x!
export estimate_sigma_g_w_max_ntc_n_e!, ntc_n_e!, ntc_n_e_es!
export create_particle_vector, create_grid1D_uniform
export create_grid_sort_inplace, sort_particles!
export create_1D_boundaries, convect_particles!

end # module merzbild
