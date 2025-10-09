# Coding conventions

## Naming conventions

1. A variable of type `ParticleIndexer` struct is to be called `particle_indexer` (since `pi` could mean a particle index),
but a variable of type `ParticleIndexerArray` is to be called `pia`.

## Argument order
In general, arguments follow from least function/process-specific to most process-specific. Macroscopic parameters come last.

Rule-of-thumb: mutable structs are passed species-specific, immutable/constant quantities are passed as a whole along with species indices, and
the function figures out internally what to use.
If a function operates only particles of a specific species, only the vector of those particles should be passed in.
Collision data should also be the specific pair used in a function.

Species data and interaction data should be passed in full (the whole vectors), not the specific
pair used in a function.

Product species follow reactant species (same applies to particle arrays).

1. `rng` is to be the first argument in a function signature

2. This is then followed by any function-specific structs.
For example, a merging setting struct for particle merging, or collision structs
for collisions: `collision_factors`, `collision_data`, `interaction`,
followed by any more specific data (e.g. `n_e_interactions`, `n_e_cs`).

3. If particle arrays are passed to the function, they come next

4. Then comes the particle indexing struct `pia` 

5. Then come the cell and species indices: `cells` are to be followed by `species`

6. Then comes the species data list

7. Then the `phys_props` struct

8. Then the `surf_props` struct

9. Then the `flux_props` struct

10. Then any function-specific non-struct parameters

For collisions with possible reactant products (e.g. ionizing collisions), first everything reactant-related,
then everything product-related.

## Argument order examples

1. Single-species NTC routine: `ntc!(rng, collision_factors, collision_data, interaction, particles, pia, cell, species, Δt, V)`

2. Multi-species NTC routine: `ntc!(rng, collision_factors, collision_data, interaction, particles_1, particles_2, pia, cell, species_1, species_2, Δt, V)`

3. Electron-neutral ionizing collision routine with even splitting: `ntc_n_e_es!(rng, collision_factors, collision_data, interaction n_e_interactions, n_e_cs, particles_n, particles_e, particles_ion, pia, cell, species_n, species_e, species_ion, Δt, V)`

4. Getting electron-neutral cross-section data: `get_cs_ionization(electron_neutral_interactions, computed_cs, neutral_species_index)`
5. Octree-based N:2 merging: `merge_octree_N2_based!(octree, particles, pia, cell, species, target_np)`

6. Computation of physical properties: `compute_props!(particles, pia, species_data, phys_props)`

# Checklist before contributing

- [ ] Do all tests pass?

- [ ] Does the new feature have corresponding tests added?

- [ ] Is it documented?

- [ ] Is the documentation added to the correct section in the `docs/src/reference.md` file?

# Versioning
The package uses semantic versioning. The MINOR version number is bumped for breaking changes; the PATCH version number is bumped when new features and/or bug fixes are added.
