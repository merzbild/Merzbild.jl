module Merzbild

using MuladdMacro
using StaticArrays

include("abstract_types.jl")
include("constants.jl")
include("utils.jl")
include("particles.jl")
include("distributions_and_sampling.jl")
include("properties/physical_props.jl")
include("properties/moment_props.jl")
include("properties/flux_props.jl")
include("properties/collisional_props.jl")
include("collisions/collisions.jl")
include("grids/grids.jl")
include("merging/merging.jl")
include("pic/pic.jl")
include("properties/surface_props.jl")
include("io.jl")
include("io_moments.jl")
include("convection/convection.jl")
include("parallel/parallel.jl")

const OCTREE_DEFAULT_BUFFER_SIZE::Int64 = 8192
const DELTA_PARTICLES::Int64 = 256

"""
    MERZBILD_DATA_PATH

Absolute path to the `data` directory bundled with Merzbild.jl, holding the built-in species
(`particles.toml`) and interaction (`vhs.toml`, `pseudo_maxwell.toml`) data files. Use it to load
bundled data independently of the current working directory, e.g.
`load_species_data(joinpath(MERZBILD_DATA_PATH, "particles.toml"), "Ar")`.
"""
const MERZBILD_DATA_PATH = joinpath(pkgdir(Merzbild), "data")
export MERZBILD_DATA_PATH

"""
    MERZBILD_SIMULATIONS_PATH

Absolute path to the `simulations` directory bundled with Merzbild.jl, holding various simulation examples.
Use it to run bundled simulations independently of the current working directory, e.g.
`include(joinpath(MERZBILD_SIMULATIONS_PATH, "0D/BKW/bkw.jl"))`.
"""
const MERZBILD_SIMULATIONS_PATH = joinpath(pkgdir(Merzbild), "simulations")
export MERZBILD_SIMULATIONS_PATH

"""
    MERZBILD_SCRIPTS_PATH

Absolute path to the `scripts` directory bundled with Merzbild.jl, holding various Python scripts for
data post-processing and plotting.
"""
const MERZBILD_SCRIPTS_PATH = joinpath(pkgdir(Merzbild), "scripts")
export MERZBILD_SCRIPTS_PATH

"""
Relative tolerance below which a post-merge variance is treated as having collapsed to zero,
see [`variance_scaling`](@ref).
"""
const variance_scaling_rel_tol = 1e-10

export sample_maxwellian_on_grid!, sample_on_grid!, bkw, maxwellian
export load_species_data, Particle, sample_particles_equal_weight!
export sample_particles_phase_box_weighted!
export Species, Interaction
export compute_props!, compute_props_sorted!
export compute_moment_scaling!, compute_moments!
export clear_props!, avg_props!
export ParticleIndexer, ParticleIndexerArray
export PhysProps
export CollisionFactors, CollisionFactorsSWPM
export CollisionData, CollisionDataFP
export SurfProps, reduce_surf_props!
export squash_pia!
export NCDataHolder, NCDataHolderSurf, IOSkipList, IOSkipListSurf, NCDataHolderFlux, IOSkipListFlux, NCDataHolderMoments
export write_netcdf
export close_netcdf
export load_interaction_data, load_interaction_data_with_dummy, load_species_and_interaction_data
export create_collision_factors_array
export create_collision_factors_swpm_array
export estimate_sigma_g_w_max, estimate_sigma_g_w_max!, estimate_sigma_g_max!
export ntc!, fp_linear!, swpm!, ntc_equal_weight!
export k_B
export GridN2Merge, merge_grid_based!
export OctreeMerge, OctreeBinMidSplit, OctreeBinMeanSplit, OctreeBinMedianSplit
export OctreeInitBinMinMaxVel, OctreeInitBinMinMaxVelSym, OctreeInitBinC
export OctreeBinBoundsInherit, OctreeBinBoundsRecompute
export merge_octree!
export merge_roulette!
export NNLSMerge, compute_multi_index_moments, merge_nnls_based!
export merge_nnls_based_rate_preserving!
export ElectronNeutralInteractions, ComputedCrossSections
export create_computed_crosssections, DataMissingException
export ScatteringIsotropic, ScatteringOkhrimovskyy
export ElectronEnergySplitEqual, ElectronEnergySplitZeroE
export CSExtendZero, CSExtendConstant
export accelerate_constant_field_x!
export estimate_sigma_g_w_max_ntc_n_e!, ntc_n_e!, ntc_n_e_es!
export ParticleVector
export AbstractGrid, Grid1DUniform, write_grid
export GridSortInPlace, sort_particles!
export FullyDiffuseBC1D, MaxwellWallBC1D, FullySpecularBC1D, convect_particles!, convect_particles_and_compute_cell!
export convect_particles_periodic!, convect_particles_and_compute_cell_periodic!
export pretty_print_pia
export ChunkExchanger, exchange_particles!, reset!, sort_particles_after_exchange!
export count_disordered_particles, check_pia_is_correct, check_unique_index, check_unique_buffer
export FluxProps, compute_flux_props!, compute_flux_props_sorted!
export mean_free_path, mean_collision_frequency
export generate_1_factorization
export restore_particle_ordering!

end # module merzbild
