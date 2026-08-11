# Merzbild.jl internal API reference

The public API functions are exported by the module

## Particle indexing
```@docs
Merzbild.map_cont_index
Merzbild.update_particle_indexer_new_particle!
Base.getindex
Base.setindex!
Base.length
Base.resize!
Merzbild.add_particle!
Merzbild.update_particle_buffer_new_particle!
Merzbild.update_buffer_index_new_particle!
Merzbild.delete_particle!
Merzbild.delete_particle_end!
Merzbild.delete_particle_end_group1!
Merzbild.delete_particle_end_group2!
Merzbild.delete_batch_end!
Merzbild.delete_batch_end_group1!
Merzbild.delete_batch_end_group2!
Merzbild.find_index_last_after_group1_delete!
Merzbild.find_index_last_after_group2_delete!
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
Merzbild.norm3
Merzbild.compute_n_coll_single_species
Merzbild.compute_n_coll_two_species
Merzbild.collide_2particles_vhs!
Merzbild.collide_2particles_vhs_equal_weight!
Merzbild.compute_vhs_factor
Merzbild.compute_com!
Merzbild.compute_g!
Merzbild.compute_g_new_ionization!
Merzbild.scatter_vhs!
Merzbild.scatter_electron_vhs!
Merzbild.scatter_ionization_electrons_and_ion!
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

## Fokker-Planck computations
```@docs
Merzbild.sample_normal_rands!
Merzbild.scale_norm_rands!
Merzbild.compute_relaxation_time
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
Merzbild.write_back_to_particles!
Merzbild.variance_scaling_rel_tol
Merzbild.variance_scaling
Merzbild.vx_sign
Merzbild.vy_sign
Merzbild.vz_sign
Merzbild.base_multi_index_moments
Merzbild.compute_w_total_v0!
Merzbild.fill_powers!
Merzbild.ccm_vel
Merzbild.ccm_pos
Merzbild.compute_lhs_and_rhs!
Merzbild.compute_lhs_and_rhs_rate_preserving!
Merzbild.scale_lhs_rhs_variance!
Merzbild.scale_lhs_rhs_vref!
Merzbild.scale_lhs_rhs_spatial_variance!
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
Merzbild.convect_single_particle_periodic!
Merzbild.set_x
```

## Particle-surface interactions
```@docs
Merzbild.AbstractBC
Merzbild.specular_reflection_x!
Merzbild.diffuse_reflection_x!
Merzbild.update_surface_incident!
Merzbild.update_surface_reflected!
Merzbild.surface_props_scale!
Merzbild.apply_bc!
Merzbild.apply_bc_dispatched!
```

## I/O
```@docs
Merzbild.AbstractNCDataHolder
```

## Parallel computations
```@docs
Merzbild.swap_particles_true_index!
Merzbild.swap_particles!
Merzbild.update_swap_indexing!
Merzbild.push_particles!
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
Merzbild.create_position_vector_D
Merzbild.compute_thermal_velocity
Merzbild.binary_search
Merzbild.linear_interpolation
Merzbild.compute_mixed_moment
Merzbild.scale_columns!
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
Merzbild.apply_householder_sweep!
```
