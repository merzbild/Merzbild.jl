# Merzbild.jl public API reference

## Particles
```@docs
Particle
ParticleVector
ParticleVector(np)
Base.getindex(pv::ParticleVector, i)
Base.setindex!(pv::ParticleVector, p::Particle, i::Integer)
Base.length(pv::ParticleVector)
Base.resize!(pv::ParticleVector, n::Integer)
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
count_disordered_particles
check_unique_index
check_pia_is_correct
pretty_print_pia
```

## Loading species and interaction data
```@docs
Species
Interaction
load_species_data
load_interaction_data
load_interaction_data_with_dummy
load_species_and_interaction_data
load_electron_neutral_interactions
```

## Sampling
```@docs
maxwellian
bkw
sample_on_grid!
sample_maxwellian_on_grid!
sample_particles_equal_weight!
```

## Computing grid and surface macroscopic properties
```@docs
PhysProps
PhysProps(n_cells, n_species, moments_list; Tref=300.0)
PhysProps(pia, moments_list; Tref=300.0)
PhysProps(pia)
SurfProps
SurfProps(n_elements, n_species, areas, normals)
SurfProps(pia, grid::Grid1DUniform)
FluxProps
FluxProps(n_cells, n_species)
FluxProps(pia)
compute_props!
compute_props_with_total_moments!
compute_props_sorted!
compute_flux_props!
compute_flux_props_sorted!
avg_props!
clear_props!
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
create_collision_factors_array(n_species)
create_collision_factors_array(n_species, n_cells)
create_collision_factors_array(pia::ParticleIndexerArray)
create_collision_factors_array(pia, interactions, species_data, T::Real, Fnum::Real; mult_factor=1.0)
create_collision_factors_array(pia, interactions, species_data, T_list, Fnum::Real; mult_factor=1.0)
create_collision_factors_swpm_array(n_species)
create_collision_factors_swpm_array(n_species, n_cells)
create_collision_factors_swpm_array(pia::ParticleIndexerArray)
create_collision_factors_swpm_array(pia, interactions, species_data, T::Real; mult_factor=1.0)
create_collision_factors_swpm_array(pia, interactions, species_data, T_list; mult_factor=1.0)
create_computed_crosssections
estimate_sigma_g_w_max
estimate_sigma_g_w_max!
estimate_sigma_g_w_max_ntc_n_e!
estimate_sigma_g_max!
ntc!(rng, collision_factors, collision_data, interaction, particles, pia, cell, species, Δt, V)
ntc!(rng, collision_factors, collision_data, interaction,
              particles_1, particles_2, pia,
              cell, species1, species2, Δt, V)
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
ComputedCrossSections
Merzbild.ElectronEnergySplit
Merzbild.ScatteringLaw
Merzbild.CSExtend
```

## Merging

### Grid merging
```@docs
GridN2Merge
GridN2Merge(Nx::Int, Ny::Int, Nz::Int, extent_multiplier::T) where T <: AbstractArray
GridN2Merge(N::Int, extent_multiplier::T) where T <: AbstractArray
GridN2Merge(Nx::Int, Ny::Int, Nz::Int, extent_multiplier::Float64)
GridN2Merge(Nx::Int, Ny::Int, Nz::Int, extent_multiplier_x::Float64, extent_multiplier_y::Float64, extent_multiplier_z::Float64)
GridN2Merge(N::Int, extent_multiplier::Float64)
merge_grid_based!
```

### NNLS merging
```@docs
NNLSMerge
NNLSMerge(multi_index_moments, init_np; rate_preserving=false)
compute_multi_index_moments
merge_nnls_based!
merge_nnls_based_rate_preserving!
```

### Octree merging
```@docs
Merzbild.OctreeBinSplit
Merzbild.OctreeInitBin
Merzbild.OctreeBinBounds
OctreeN2Merge
OctreeN2Merge(split::Merzbild.OctreeBinSplit;
              init_bin_bounds=OctreeInitBinMinMaxVel,
              bin_bounds_compute=OctreeBinBoundsInherit,
              max_Nbins=4096, max_depth=10)
merge_octree_N2_based!
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
```

## Particle-surface interactions
```@docs
MaxwellWallBC
MaxwellWalls1D
MaxwellWalls1D(species_data, T_l::Float64, T_r::Float64, vy_l::Float64, vy_r::Float64, accomodation_l::Float64,          
    accomodation_r::Float64)
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
NCDataHolder(nc_filename, names_skip_list, species_data, phys_props; global_attributes=Dict{Any,Any}())
NCDataHolder(nc_filename, species_data, phys_props; global_attributes=Dict{Any,Any}())
NCDataHolderSurf
NCDataHolderSurf(nc_filename, names_skip_list, species_data, surf_props; global_attributes=Dict{Any,Any}())
NCDataHolderSurf(nc_filename, species_data, surf_props; global_attributes=Dict{Any,Any}())
NCDataHolderFlux
NCDataHolderFlux(nc_filename, names_skip_list, species_data, flux_props; global_attributes=Dict{Any,Any}())
NCDataHolderFlux(nc_filename, species_data, flux_props; global_attributes=Dict{Any,Any}())
write_netcdf(ds, phys_props::PhysProps, timestep; sync_freq=0)
write_netcdf(ds, surf_props::SurfProps, timestep; sync_freq=0)
write_netcdf(ds, flux_props::FluxProps, timestep; sync_freq=0)
write_netcdf(nc_filename, pv::ParticleVector, pia, species, species_data; global_attributes=Dict{Any,Any}())
write_netcdf(nc_filename, particles::Vector{ParticleVector}, pia, species_data; global_attributes=Dict{Any,Any}())
close_netcdf
```

## Parallel computations
```@docs
ChunkExchanger
ChunkExchanger(chunks, n_cells) 
exchange_particles!
sort_particles_after_exchange!
reset!
reduce_surf_props!
```

## Particle-in-Cell
```@docs
accelerate_constant_field_x!
```

## Constants
```@docs
k_B
```

## Misc
```@docs
DataMissingException
```