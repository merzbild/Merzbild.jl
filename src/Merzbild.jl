module Merzbild

using MuladdMacro

include("abstract_types.jl")
include("constants.jl")
include("utils.jl")
include("particles.jl")
include("distributions_and_sampling.jl")
include("physical_props.jl")
include("flux_props.jl")
include("collisions/collisions.jl")
include("grids/grids.jl")
include("merging/merging.jl")
include("pic/pic.jl")
include("surface_props.jl")
include("io.jl")
include("convection/convection.jl")
include("parallel.jl")

const OCTREE_DEFAULT_BUFFER_SIZE::Int32 = 8192
const DELTA_PARTICLES::Int32 = 256

export sample_maxwellian_on_grid!, sample_on_grid!, bkw, maxwellian
export load_species_data, Particle, sample_particles_equal_weight!
export Species, Interaction
export compute_props!, compute_props_sorted!, compute_props_with_total_moments!
export clear_props!, avg_props!
export ParticleIndexer, ParticleIndexerArray, PhysProps, CollisionFactors, CollisionData, CollisionDataFP
export SurfProps, reduce_surf_props!
export squash_pia!
export NCDataHolder, NCDataHolderSurf, IOSkipList, IOSkipListSurf, NCDataHolderFlux, IOSkipListFlux
export write_netcdf
export close_netcdf
export load_interaction_data, load_interaction_data_with_dummy, load_species_and_interaction_data
export create_collision_factors_array
export estimate_sigma_g_w_max, estimate_sigma_g_w_max!
export ntc!, fp_linear!
export k_B
export GridN2Merge, merge_grid_based!
export OctreeN2Merge, OctreeBinMidSplit, OctreeBinMeanSplit, OctreeBinMedianSplit
export OctreeInitBinMinMaxVel, OctreeInitBinMinMaxVelSym, OctreeInitBinC
export OctreeBinBoundsInherit, OctreeBinBoundsRecompute
export merge_octree_N2_based!
export merge_roulette!
export NNLSMerge, compute_multi_index_moments, merge_nnls_based!
export merge_nnls_based_rate_preserving!
export ElectronNeutralInteractions, ComputedCrossSections
export load_electron_neutral_interactions, create_computed_crosssections, DataMissingException
export ScatteringIsotropic, ScatteringOkhrimovskyy
export ElectronEnergySplitEqual, ElectronEnergySplitZeroE
export CSExtendZero, CSExtendConstant
export accelerate_constant_field_x!
export estimate_sigma_g_w_max_ntc_n_e!, ntc_n_e!, ntc_n_e_es!
export ParticleVector
export AbstractGrid, Grid1DUniform, write_grid
export GridSortInPlace, sort_particles!
export MaxwellWallBC, MaxwellWalls1D, convect_particles!
export pretty_print_pia
export ChunkExchanger, exchange_particles!, reset!, sort_particles_after_exchange!
export count_disordered_particles, check_pia_is_correct, check_unique_index
export FluxProps, compute_flux_props!, compute_flux_props_sorted!

end # module merzbild
