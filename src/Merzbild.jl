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
export load_species_data, Particle, sample_particles_equal_weight!
export Species, Interaction
export compute_props!, compute_props_sorted!, compute_props_with_total_moments!
export clear_props!, avg_props!
export ParticleIndexer, ParticleIndexerArray, PhysProps, CollisionFactors, CollisionData, CollisionDataFP
export NCDataHolder, IOSkipList
export write_netcdf_phys_props
export close_netcdf
export load_interaction_data, load_interaction_data_with_dummy, load_species_and_interaction_data
export create_collision_factors_array
export estimate_sigma_g_w_max, estimate_sigma_g_w_max!
export ntc!, fp!
export k_B
export GridN2Merge, merge_grid_based!
export OctreeN2Merge, OctreeBinMidSplit, OctreeBinMeanSplit, OctreeBinMedianSplit
export OctreeInitBinMinMaxVel, OctreeInitBinMinMaxVelSym, OctreeInitBinC
export OctreeBinBoundsInherit, OctreeBinBoundsRecompute
export merge_octree_N2_based!
export NNLSMerge, compute_multi_index_moments, merge_nnls_based!
export merge_nnls_based_rate_preserving!
export load_electron_neutral_interactions, create_computed_crosssections, DataMissingException
export ScatteringIsotropic, ScatteringOkhrimovskyy
export ElectronEnergySplitEqual, ElectronEnergySplitZeroE
export accelerate_constant_field_x!
export estimate_sigma_g_w_max_ntc_n_e!, ntc_n_e!, ntc_n_e_es!
export ParticleVector
export Grid1DUniform, write_grid
export GridSortInPlace, sort_particles!
export MaxwellWallBC, MaxwellWalls1D, convect_particles!

end # module merzbild
