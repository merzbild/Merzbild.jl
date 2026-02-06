# Debugging Merzbild

Debugging a particle code is no mean feat. The unit tests in Merzbild try to cover all of the functionality,
however, occasional bugs might still be present. Indexing bugs are probably the most annoying ones, and have been encountered
several times during development, since the structures used for indexing are non-trivial, and things like merging and
multi-threading do a lot of indexing work.

Therefore, several utility functions have been developed to help track down indexing-related bugs.

[`pretty_print_pia`](@ref) prints out the particle indexing for a given species over all cells, producing output like
```
Total: 9
Cell 1: group1: [1, 3] group2: [4, 5]
Cell 3: group2: [6, 8]
Cell 4: group1: [9, 9]
```
This can be helpful in checking that indexing doesn't overlap or in general looking at what happens to specific groups of particles in
certain cells.

One can also use [`check_pia_is_correct`](@ref) to verify that the start of the indexing in a cell follows the end of indexing in the
previous cell, without gaps. It also checks that local cell particle counts sum up to the total cell count.

[`check_unique_index`](@ref) verifies that the underlying indexing of a `ParticleVector` instance doesn't contain duplicates,
i.e. that `pv[i]` points to a different particle than `pv[j]` if `i!=j`.

Finally, [`check_unique_buffer`](@ref) checks whether the indexes stored in the buffer (that are used to instantiate new particles),
specifically, stored in the active part of the buffer, are all unique.
