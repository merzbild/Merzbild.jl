# Merzbild.jl public API reference

## Particles
```@docs
Particle
ParticleVector
ParticleVector(np)
Base.getindex(pv::ParticleVector{D}, i) where D
Base.setindex!(pv::ParticleVector{D}, p::Particle{D}, i::Integer) where D
Base.length(pv::ParticleVector{D}) where D
Base.resize!(pv::ParticleVector{D}, n::Integer) where D
```

## Particle indexing
```@docs
ParticleIndexer
ParticleIndexer()
ParticleIndexer(n_particles)
ParticleIndexerArray
ParticleIndexerArray(indexer_arr::Array{ParticleIndexer,2}, n_total)
ParticleIndexerArray(n_cells::Integer, n_species::Integer)
ParticleIndexerArray(n_particles::Integer) 
ParticleIndexerArray(n_particles::T) where T<:AbstractVector
ParticleIndexerArray(grid, species_data::Array{Species}) 
squash_pia!
restore_particle_ordering!
count_disordered_particles
check_unique_index
check_pia_is_correct
pretty_print_pia
check_unique_buffer
```

## Loading species and interaction data
```@docs
MERZBILD_DATA_PATH
Species
Interaction
Interaction(m1::Float64, m2::Float64, vhs_d::Float64, vhs_o::Float64, vhs_Tref::Float64)
load_species_data
load_interaction_data
load_interaction_data_with_dummy
load_species_and_interaction_data
```

## Sampling
```@docs
maxwellian
bkw
sample_on_grid!
sample_maxwellian_on_grid!
sample_particles_equal_weight!
sample_particles_phase_box_weighted!
```

## Computing grid and surface macroscopic properties
```@docs
PhysProps
PhysProps(n_cells, n_species; ndens_not_Np=false)
PhysProps(pia::ParticleIndexerArray; ndens_not_Np=false)
SurfProps
SurfProps(n_elements, n_species, areas, normals)
SurfProps(pia, grid::Grid1DUniform)
FluxProps
FluxProps(n_cells, n_species)
FluxProps(pia)
ElectrostaticFieldProps
ElectrostaticFieldProps(n_nodes::Integer)
ElectrostaticFieldProps(grid::Grid1DUniform)
clear_charge_density!
compute_props!
compute_props_sorted!
compute_flux_props!
compute_flux_props_sorted!
avg_props!
clear_props!
compute_moment_scaling!
compute_moments!
```

## Collisional properties
```@docs
mean_free_path
mean_collision_frequency
debye_length
plasma_frequency
```

## Collision computations
```@docs
CollisionData
CollisionData()
CollisionFactors
CollisionFactors()
CollisionFactorsSWPM
CollisionFactorsSWPM()
CollisionDataFP
CollisionDataFP()
CollisionDataFP(n_particles_in_cell)
create_collision_factors_array
create_collision_factors_swpm_array
create_computed_crosssections
estimate_sigma_g_w_max
estimate_sigma_g_w_max!
estimate_sigma_g_w_max_ntc_n_e!
estimate_sigma_g_max!
ntc!
ntc_equal_weight!
ntc_n_e!
ntc_n_e_es!
swpm!
```

## Fokker-Planck computations
```@docs
fp_linear!
```

## Electron-neutral interactions
```@docs
ElectronNeutralInteractions
ElectronNeutralInteractions(species_data, filename, databases, scattering_laws, energy_splits)
ComputedCrossSections
Merzbild.ElectronEnergySplit
Merzbild.ScatteringLaw
Merzbild.CSExtend
```

## Merging

### Grid merging
```@docs
GridN2Merge{D}
GridN2Merge{D}(Nx::Int, Ny::Int, Nz::Int, extent_multiplier::T) where {D, T <: AbstractArray}
GridN2Merge{D}(N::Int, extent_multiplier::T) where {D, T <: AbstractArray}
GridN2Merge{D}(Nx::Int, Ny::Int, Nz::Int, extent_multiplier::Float64) where D
GridN2Merge{D}(Nx::Int, Ny::Int, Nz::Int, extent_multiplier_x::Float64, extent_multiplier_y::Float64, extent_multiplier_z::Float64) where D
GridN2Merge{D}(N::Int, extent_multiplier::Float64) where D
GridN2Merge(N::Int, extent_multiplier::T) where T <: AbstractArray 
GridN2Merge(Nx::Int, Ny::Int, Nz::Int, extent_multiplier::Float64)
GridN2Merge(Nx::Int, Ny::Int, Nz::Int, extent_multiplier_x::Float64, extent_multiplier_y::Float64, extent_multiplier_z::Float64)
merge_grid_based!
GridN2Merge(N::Int, extent_multiplier::Float64)
```

### NNLS merging
```@docs
NNLSMerge{D}
NNLSMerge{D}(multi_index_moments, init_np; rate_preserving=false, multi_index_moments_pos=[], matrix_ncol_nprealloc=0) where D
NNLSMerge(multi_index_moments, init_np; rate_preserving=false, multi_index_moments_pos=[], matrix_ncol_nprealloc=0)
compute_multi_index_moments
merge_nnls_based!
merge_nnls_based_rate_preserving!
```

### Octree merging
```@docs
Merzbild.OctreeBinSplit
Merzbild.OctreeInitBin
Merzbild.OctreeBinBounds
OctreeMerge
OctreeMerge{D,M}(split::Merzbild.OctreeBinSplit; init_bin_bounds=OctreeInitBinMinMaxVel, bin_bounds_compute=OctreeBinBoundsInherit, max_Nbins=4096, max_depth=10) where {D,M}
OctreeMerge(split::Merzbild.OctreeBinSplit; init_bin_bounds=OctreeInitBinMinMaxVel, bin_bounds_compute=OctreeBinBoundsInherit, max_Nbins=4096, max_depth=10) 
merge_octree!
```

### Roulette merging
```@docs
merge_roulette!
```

## Grids and particle sorting
```@docs
AbstractGrid
Grid1DUniform
Grid1DUniform(L, nx; wall_offset=1e-12)
GridSortInPlace
GridSortInPlace(n_cells::Integer, n_particles::Integer)
GridSortInPlace(grid::G, n_particles::Integer) where {G<:AbstractGrid}
sort_particles!
```

## Particle movement
```@docs
convect_particles!
convect_particles_periodic!
convect_particles_and_compute_cell!
convect_particles_and_compute_cell_periodic!
```

## Particle-surface interactions
```@docs
MaxwellWallBC1D
MaxwellWallBC1D(species, species_data,T::Float64, v, accommodation::Float64)
FullyDiffuseBC1D
FullyDiffuseBC1D(species, species_data, T::Float64, v)
FullySpecularBC1D
```

## I/O
```@docs
write_grid
IOSkipList
IOSkipList(list_of_variables_to_skip)
IOSkipList()
IOSkipListSurf
IOSkipListSurf(list_of_variables_to_skip)
IOSkipListSurf()
IOSkipListFlux
IOSkipListFlux(list_of_variables_to_skip)
IOSkipListFlux()
NCDataHolder
NCDataHolder(nc_filename, names_skip_list, species_data, phys_props; global_attributes=Dict{Any,Any}(), mode=NC_64BIT_OFFSET)
NCDataHolder(nc_filename, species_data, phys_props; global_attributes=Dict{Any,Any}(), mode=NC_64BIT_OFFSET)
NCDataHolderSurf
NCDataHolderSurf(nc_filename, names_skip_list, species_data, surf_props; global_attributes=Dict{Any,Any}(), mode=NC_64BIT_OFFSET)
NCDataHolderSurf(nc_filename, species_data, surf_props; global_attributes=Dict{Any,Any}(), mode=NC_64BIT_OFFSET)
NCDataHolderFlux
NCDataHolderFlux(nc_filename, names_skip_list, species_data, flux_props; global_attributes=Dict{Any,Any}(), mode=NC_64BIT_OFFSET)
NCDataHolderFlux(nc_filename, species_data, flux_props; global_attributes=Dict{Any,Any}(), mode=NC_64BIT_OFFSET)
NCDataHolderMoments
NCDataHolderMoments(nc_filename, species_data, n_cells, n_species, moment_powers; global_attributes=Dict{Any,Any}(), mode=NC_64BIT_OFFSET)
write_netcdf
close_netcdf
```

## Parallel computations
```@docs
ChunkExchanger
ChunkExchanger(chunks, n_cells) 
update_occupancy_bounds!
exchange_particles!
sort_particles_after_exchange!
reset!
reduce_surf_props!
reduce_field_props!
generate_1_factorization
LoadBalancerCellQ
LoadBalancerCellQ(n_cells, n_chunks)
update_lb_cellq!
rebalance_lb!
reset_lb!
```

## Particle-in-Cell
```@docs
accelerate_constant_field_x!
DirichletFieldBC1D
NeumannFieldBC1D
PeriodicFieldBC1D
PoissonSolver1DUniform
PoissonSolver1DUniform(grid::Grid1DUniform, bc_left::Merzbild.AbstractFieldBC1D, bc_right::Merzbild.AbstractFieldBC1D)
deposit_charge!
normalize_charge_density!
solve_poisson!
accelerate_electric_field_x!
```

## Constants
```@docs
k_B
eps_0
```

## Misc
```@docs
DataMissingException
MERZBILD_SIMULATIONS_PATH
MERZBILD_SCRIPTS_PATH
```