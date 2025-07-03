# Particle buffers and contiguous indexing

This section discusses some more advanced concepts and associated issues concerning particle indexing in scenarios where
particle counts may become lower during the course of a simulation.

## Particle buffer
As discussed in [the overview of the basic building blocks of Merzbild.jl](@ref "Overview of basic building blocks"),
a `ParticleVector` instance has a `particles` vector, an `index` array to index the particles, and `buffer` array to keep track of
which particles in the `particles` vector are unused (as a LIFO queue).

So to create a new particle in the simulation (assuming the `ParticleVector` instance `pv` has enough unused pre-allocated particles,
otherwise one needs to call `resize!` first), one needs to do the following:

1. Update the `pia` structure to keep track for which species in which cell and group the particle is being created
2. Grab the index of the appropriate element in the buffer: `i = pv.buffer[pv.nbuffer]`
    and decrease the number of elements in the buffer (`pv.nbuffer -= 1`)
3. Writes this index `i` to the last position in the `index` array (the particle count in `pia` has already been updated):
    `pv.index[pia.n_total[species]] = i`. 

A utility function [`update_buffer_index_new_particle!`](@ref Merzbild.update_buffer_index_new_particle!) is available which takes care of
all of these steps.

**Summary**: one needs to be aware of the buffer when creating new particles, to avoid any issues with over-writing existing particles.

## Contiguous indexing
One can iterate through all the particles in a simulation of a specific species by going through all the cells and then
going through all the groups of particles:
```julia
for cell in 1:grid.n_cells
    s = pia.indexer[cell, species].start1
    e = pia.indexer[cell, species].end1
    
    for i in s:e
        do_something!(particles[i])
    end

    if pia.indexer[cell, species].n_group2 > 0
        s = pia.indexer[cell, species].start2
        e = pia.indexer[cell, species].end2
    
        for i in s:e
            do_something!(particles[i])
        end
    end
end
```

However, it is more efficient to iterate in the following manner:
```julia
for i in 1:pia.n_total[species]
    do_something!(particles[i]) 
end
```

This however requires the indexing to be **contiguous**, which means the following should hold:
1. `pia.indexer[cell,species].end1 + 1 == pia.indexer[cell+1,species].start1`, `cell = 1,...,n_cells - 1`
2. `pia.indexer[n_cells,species].end1 + 1 == pia.indexer[1,species].start2`
3. `pia.indexer[cell,species].end2 + 1 == pia.indexer[cell+1,species].start2`, `cell = 1,...,n_cells - 1`.
So this basically correspond to the indexing as defined by the indexer having no "holes". For this purposes, `pia` has
the boolean `contiguous` property, which can be checked to decided how to iterate over particles.

It is expected that a particle sorting routine restores continuity of indices. There is also a utility function [`squash_pia!`](@ref squash_pia!)
which restores continuity of indices by moving around the indices in a `ParticleVector` instance, as well as the starts and ends of groups
in the `pia` particle indexing structure.

**Summary**: If iterating over particles in more than one cell, one should not assume that the particle indices are contiguous, and the value
of `pia.contiguous` should be checked first. Particle sorting or a call to `squash_pia!` restore continuity of indices.

## Deleting particles
So how do these holes in the indexing actually appear?

Deletion happens at the end of the group, if not, the ordering of the particles is changed by changing around the indices so
that the particle to be deleted moves to the end of the group and then is deleted.

A function [`delete_particle!`](@ref Merzbild.delete_particle!) is provided which does exactly this.
A more specialized function [`delete_particle_end!`](@ref Merzbild.delete_particle_end!) deletes the last particle
in the group.

So a hole might appear in the indexing if one is doing particle merging. For example, let's say we have 20 particles in two cells with 10 particles per cell, and the indexing looks like this: `pia.indexer[1,1].start1 = 1`, `pia.indexer[1,1].end1 = 10`,
`pia.indexer[2,1].start1 = 11`,  `pia.indexer[2,1].end1 = 20`. If we merge particles in cell 1 down to 2 particles,
the `buffer` in the `ParticleVector` instance `pv` will be updated, and `pia.indexer[1,1].end1` will be set to 2, but the particles
`pv[3:10]` can't be really accessed or used, unless we either 1) sort the particles 2) call [`squash_pia!`](@ref squash_pia!).
Some computations currently used the number of particles of a certain species to find where to place a newly created particle (for example,
in the variable-weight NTC collision functions); without restoring continuity of the indexing, this will lead to erroneous results,
as now the index of the last particle in the simulation `pv[20]` is no longer the same as the total number of particles in the simulation (which is 12 after the merge).

So, for a multi-dimensional simulation with variable-weight DSMC, the correct collide-merge procedure might take on the following form:
```julia

for cell in 1:grid.n_cells
    ntc!(rng, collision_factors[1, 1, cell],
         collision_data, interaction_data, particles[1], pia, cell, 1, Δt, grid.cells[cell].V)

    if pia.indexer[cell,1].n_local > merge_threshold
        # we need to merge
        merge_octree_N2_based!(rng, oc, particles[1], pia, cell, 1, merge_target, grid)
        # restore continuity of pia
        squash_pia!(particles, pia)
    end
end
```

So if no merging took place, then no particle deletion happened, and one doesn't need to call `squash_pia!`. Of course, `squash_pia!` checks
the value of `pia.contiguous` and does nothing if that value is set to `true`.

**Summary**: one needs to restore continuity of particle indexing if particles are deleted and created in a simulation, otherwise this
might lead to erroneous results.