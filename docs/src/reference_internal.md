# Merzbild.jl internal API reference

The public API functions are exported by the module

## Particle indexing
```@docs
Merzbild.map_cont_index
Merzbild.update_particle_indexer_new_lower_count!
Merzbild.update_particle_indexer_new_particle!
Base.getindex
Base.setindex!
Base.length
Base.resize!
Merzbild.add_particle!
Merzbild.update_particle_buffer_new_particle!(pv::ParticleVector, position)
Merzbild.update_particle_buffer_new_particle!(pv::ParticleVector, pia, species)
Merzbild.update_particle_buffer_new_particle!(pv::Vector{Particle}, pia, species)
Merzbild.update_particle_buffer_new_particle!(pv::Vector{Particle}, position)
Merzbild.update_buffer_index_new_particle!
Merzbild.delete_particle!
Merzbild.delete_particle_end!
Merzbild.delete_particle_end_group1!
Merzbild.delete_particle_end_group2!
```

## Loading species and interaction data
```@docs
Merzbild.compute_mu_ref
```

## Sampling
```@docs
Merzbild.UnitDVGrid
Merzbild.DVGrid
Merzbild.VDF
Merzbild.sample_bkw!
Merzbild.evaluate_distribution_on_grid!
Merzbild.sample_maxwellian!
```

## Collision computations
```@docs
Merzbild.compute_n_coll_single_species
Merzbild.compute_n_coll_two_species
Merzbild.compute_vhs_factor
Merzbild.compute_com!
Merzbild.compute_g!
Merzbild.compute_g_new_ionization!
Merzbild.scatter_vhs!
Merzbild.scatter_electron_vhs!
Merzbild.scatter_ionization_electrons!
Merzbild.sigma_vhs
Merzbild.compute_tabulated_cs_constant_continuation
Merzbild.compute_tabulated_cs_zero_continuation
Merzbild.compute_cross_sections_only!
Merzbild.compute_cross_sections!
Merzbild.get_cs_total
Merzbild.get_cs_elastic
Merzbild.get_cs_ionization
Merzbild.get_ionization_threshold
Merzbild.get_electron_energy_split
```

## Electron-neutral interactions
```@docs
Merzbild.TabulatedCSData
Merzbild.ElasticScattering
Merzbild.ExcitationSink
Merzbild.Ionization
Merzbild.find_species_in_db
Merzbild.load_elastic_data
Merzbild.load_ionization_data
```

## Merging
```@docs
Merzbild.GridCell
Merzbild.compute_velocity_extent!
Merzbild.compute_grid_index
Merzbild.clear_merging_grid!
Merzbild.compute_grid!
Merzbild.compute_new_particles!
Merzbild.vx_sign
Merzbild.vy_sign
Merzbild.vz_sign
Merzbild.check_speed_bound
Merzbild.base_multi_index_moments
Merzbild.compute_w_total_v0!
Merzbild.ccm
Merzbild.compute_lhs_and_rhs!
Merzbild.compute_lhs_and_rhs_rate_preserving!
Merzbild.compute_lhs_particles_additional!
Merzbild.compute_lhs_particles_additional_rate_preserving!
Merzbild.scale_lhs_rhs!
Merzbild.scale_lhs_rhs_rate_preserving!
Merzbild.compute_post_merge_particles_nnls!
Merzbild.OctreeCell
Merzbild.OctreeFullCell
Merzbild.fill_bins
Merzbild.fill_full_bins
Merzbild.clear_octree!
Merzbild.resize_octree_buffers!
Merzbild.compute_octant
Merzbild.bin_bounds_inherit!
Merzbild.bin_bounds_recompute!
Merzbild.compute_v_mean!
Merzbild.compute_v_median!
Merzbild.get_new_bin_id
Merzbild.split_bin!
Merzbild.compute_bin_props!
Merzbild.get_bin_post_merge_np
Merzbild.init_octree!
Merzbild.compute_octree!
```

## Grids
```@docs
Merzbild.Cell1D
Merzbild.Cell1D(xlo, xhi, V)
Merzbild.get_cell
```

## Particle movement
```@docs
Merzbild.convect_single_particle!
```

## Particle-surface interactions
```@docs
Merzbild.specular_reflection_x!
Merzbild.diffuse_reflection_x!
Merzbild.reflect_particle_x!
Merzbild.update_surface_incident!
Merzbild.update_surface_reflected!
Merzbild.surface_props_scale!
```

## I/O
```@docs
Merzbild.AbstractNCDataHolder
```

## Constants
```@docs
Merzbild.c_light
Merzbild.eV
Merzbild.eV_J
Merzbild.eV_J_inv
Merzbild.twopi
Merzbild.e_mass_div_electron_volt
Merzbild.direction_signs
Merzbild.q_e
```

## Misc
```@docs
Merzbild.compute_thermal_velocity
Merzbild.binary_search
Merzbild.linear_interpolation
Merzbild.compute_mixed_moment
```

## NNLS

```@docs
Merzbild.solve!
Merzbild.construct_householder!
Merzbild.fastview
Merzbild.solve_triangular_system!
Merzbild.UnsafeVectorView
Merzbild.orthogonal_rotmat
Merzbild.apply_householder!
```
