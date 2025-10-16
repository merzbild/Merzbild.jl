# Multithreaded simulations

## Overview
Multithreaded simulations are possible via use of standard Julia threads.
The approach is relatively straightforward:
1. The domain is decomposed into `n_chunks` "cell chunks", which are lists of cells.
    Within each chunk the cell indices should be continuous.
    For example, a 4-cell domain may be decomposed into `[[1,2],[3,4]]` (2 chunks),
    `[[1], [2,3], [4]` (3 chunks), `[[1], [2], [3], [4]]` (4 chunks), etc.
2. Several variables are instantiated for each chunk:
    - `n_chunks` arrays of `ParticleVector`s (each hold `n_species` `ParticleVector`s)
    - `n_chunks` RNGs
    - `n_chunks` `ParticleIndexerArray` instances
    - `n_chunks` `GridSortInPlace` instances
    - `n_chunks+1` `SurfProps` instances (one additional instance is required for the reduce operation)
3. Collisions, convection, sorting are all performed per-chunk, completely independently,
    so multithreading can be used
4. After particles have been sorted, some might need to be moved between chunks,
    if they ended up in cells assigned to another chunk.
    This is done in two steps: first, one calls [`exchange_particles!`](@ref) to move the particle
    data between chunks. This is done in serial mode to avoid race conditions.
    Then, [`sort_particles_after_exchange!`](@ref) is called to reset indexing without having to
    completely resort all the new particles. This can be done using multithreading.
5. Physical grid properties are computed using multithreading, as they can easily be computed
    only for cells assigned to the chunk, thus avoiding any race conditions.
    In case surface properties were computed during particle movement, a
    reduce operation is needed to sum the per-chunk values. This is done by
    a serial call to [`reduce_surf_props!`](@ref).

The trade-off of the approach is that even in case no new particles are created in the simulation,
particle data will be moved around, as when after a movement step a particle ends up in a cell
corresponding to a different chunk, it needs to be moved to the `ParticleVector` for that chunk.
In addition, some duplication of data structures (i.e. multiple `GridSortInPlace` instances)
is required. A special struct [`ChunkExchanger`](@ref) is used to facilitate exchange of particles between chunks.

However, the code logic is more straightforward, as most of the data is independent and no race conditions
can occur. Some of this functionality will be re-used to make MPI simulations possible.

## Example: multithreaded Couette flow simulation

Below is an example of a fixed-weight DSMC multithreaded Couette flow simulation. To run in multithreaded mode, one should start
julia specifying the number of threads: `julia --threads NTHREADS`. The file can also be found under `simulations/1D/couette_multithreaded.jl`.
A variable-weight example can be found under `simulations/1D/couette_multithreaded_varweight_octree.jl`.

The [ChunkSplitters.jl](https://github.com/JuliaFolds2/ChunkSplitters.jl)
library is used to perform domain decomposition, by splitting the range of cell indices `1:nx` into independent chunks.
By default, the number of chunks is set equal to the number of threads; it can be set to a multiple of the number of threads
by setting a value of the `chunk_count_multiplier` parameter. `preallocation_margin_multiplier` allocates additional unused
particles in the per-chunk `ParticleVector` instances, since otherwise the transfer of particles might lead to frequent
calls to `resize!` at the start of the simulation as the solution approaches steady state and the average number of particles in a chunk
changes significantly. A value of `1.0` means no additional particle storage is allocated.

The surface properties are collected into `surf_props_reduced` via a call to `reduce_surf_props`.
Before the start of the time loop, the sampling procedure is multithreaded via the `@threads` macro.
The physical properties are also computed in multithreaded mode.

Inside the time loop, collisions, convection, and sorting are performed inside a `@threads` block. The `chunk_exchanger` data
is also cleared in this multithreaded loop to prepare it for the movement of particles between chunks.
Once the block finishes, the particles are moved between chunks via a serial call to [`sort_particles_after_exchange!`](@ref).
Finally, the indexing is reset, and physical grid properties are computed, again inside a `@threads` block.
The reduction operation for the surface properties, as well as averaging of grid and surface properties, is performed serially
at the end of the timestep.

```julia
using Merzbild
using Random
using TimerOutputs
using Base.Threads
using ChunkSplitters

function run(seed, T_wall, v_wall, L, ndens, nx, ppc, Δt, n_timesteps, avg_start; chunk_count_multiplier=1,
             preallocation_margin_multiplier=1.0)
    reset_timer!()

    n_threads = Threads.nthreads()
    n_chunks = n_threads * chunk_count_multiplier
    println("Running on $n_threads threads, will split cells into $n_chunks chunks")

    rng_chunks = [Xoshiro(seed + i) for i in 0:n_chunks-1]

    # load particle and interaction data
    particles_data_path = joinpath("data", "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")
    interaction_data_path = joinpath("data", "vhs.toml")
    interaction_data::Array{Interaction, 2} = load_interaction_data(interaction_data_path, species_data)

    # create our grid and BCs
    grid = Grid1DUniform(L, nx)
    boundaries = MaxwellWalls1D(species_data, T_wall, T_wall, -v_wall, v_wall, 1.0, 1.0)

    # split cell indices into chunks
    cell_indices = Vector(1:nx)
    cell_chunks = chunks(cell_indices; n=n_chunks)

    # init per-chunk particle vectors, particle indexers, grid particle sorters
    n_particles_chunks = [floor(Int64, ppc * length(cell_chunk) * preallocation_margin_multiplier) for cell_chunk in cell_chunks]
    particles_chunks = [[ParticleVector(n_particles)] for n_particles in n_particles_chunks]
    pia_chunks = [ParticleIndexerArray(grid.n_cells, 1) for cell_chunk in cell_chunks]
    gridsorter_chunks = [GridSortInPlace(grid, n_particles) for n_particles in n_particles_chunks]

    # this is used for moving particles between chunks after they have been sorted into grid cells
    chunk_exchanger = ChunkExchanger(cell_chunks, nx)

    # sample particles
    # Fnum * ppc = Np in cell = ndens * V_cell
    Fnum = grid.cells[1].V * ndens / ppc

    # sample particles per-chunk
    @timeit "sampling" @threads for (chunk_id, cell_chunk) in enumerate(cell_chunks)
        @inbounds sample_particles_equal_weight!(rng_chunks[chunk_id], grid, particles_chunks[chunk_id][1],
                                                    pia_chunks[chunk_id],
                                                    1, species_data, ndens, T_wall, Fnum, cell_chunk)
    end
     
    # create collision structs
    collision_data = [CollisionData() for cell_chunk in cell_chunks]
    
    # create struct for computation of physical properties, sizes of pia are the same
    phys_props = PhysProps(pia_chunks[1])

    # create second struct for averaging of physical properties, sizes of pia are the same
    phys_props_avg = PhysProps(pia_chunks[1])

    # create struct for computation of surface properties, need a SurfProps instance per chunk
    surf_props_chunks = [SurfProps(pia_chunks[1], grid) for cell_chunk in cell_chunks]

    # we sum up all the surf props here
    surf_props_reduced = SurfProps(pia_chunks[1], grid)

    # create second struct for averaging of physical properties
    surf_props_avg = SurfProps(pia_chunks[1], grid)

    # create struct for netCDF for time-averaged physical properties I/O
    ds_avg = NCDataHolder("scratch/data/avg_mt_couette_$(L)_$(nx)_$(v_wall)_$(T_wall)_$(ppc)_after$(avg_start).nc",
                          species_data, phys_props)

    # create struct for netCDF for time-averaged surface properties I/O
    ds_surf_avg = NCDataHolderSurf("scratch/data/avg_mt_couette_$(L)_$(nx)_$(v_wall)_$(T_wall)_$(ppc)_surf_after$(avg_start).nc",
                                   species_data, surf_props_avg)

    # create and estimate collision factors
    collision_factors = [create_collision_factors_array(pia, interaction_data, species_data, T_wall, Fnum)
                         for pia in pia_chunks]

    # compute data at t=0
    @timeit "props compute" @threads for (chunk_id, cell_chunk) in enumerate(cell_chunks)
        compute_props_sorted!(particles_chunks[chunk_id], pia_chunks[chunk_id], species_data, phys_props, cell_chunk)
    end

    n_avg = n_timesteps - avg_start + 1

    for t in 1:n_timesteps

        if t % 500 == 0
            println(t)
        end
        
        # collide, convect, sort particles
        @timeit "collide+convect+sort" @threads for chunk_id in 1:n_chunks
            for cell in cell_chunks[chunk_id]
                @inbounds ntc!(rng_chunks[chunk_id], collision_factors[chunk_id][1, 1, cell],
                               collision_data[chunk_id], interaction_data, particles_chunks[chunk_id][1],
                               pia_chunks[chunk_id], cell, 1, Δt, grid.cells[cell].V)
            end

            if (t >= avg_start)
                @inbounds convect_particles!(rng_chunks[chunk_id], grid, boundaries,
                                    particles_chunks[chunk_id][1], pia_chunks[chunk_id],
                                    1, species_data, surf_props_chunks[chunk_id], Δt)
            else
                # we don't need to compute surface properties before we start averaging
                @inbounds convect_particles!(rng_chunks[chunk_id], grid, boundaries,
                                    particles_chunks[chunk_id][1], pia_chunks[chunk_id],
                                    1, species_data, Δt)
            end

            # need to clear the data in the chunk exchanger
            @inbounds reset!(chunk_exchanger, chunk_id)

            # sort particles
            @inbounds sort_particles!(gridsorter_chunks[chunk_id], grid, particles_chunks[chunk_id][1], pia_chunks[chunk_id], 1)
        end

        # move particles between chunks
        @timeit "exchange" exchange_particles!(chunk_exchanger, particles_chunks, pia_chunks, cell_chunks, 1)

        # reset indexing, compute physical properties if needed
        @timeit "re-sort + compute props" @threads for chunk_id in 1:n_chunks
            sort_particles_after_exchange!(chunk_exchanger, gridsorter_chunks[chunk_id],
                                           particles_chunks[chunk_id][1], pia_chunks[chunk_id],
                                           cell_chunks[chunk_id], 1)
            if (t >= avg_start)
                @inbounds compute_props_sorted!(particles_chunks[chunk_id], pia_chunks[chunk_id],
                                                species_data, phys_props, cell_chunks[chunk_id])
            end
        end

        # reduce surface properties, average grid and surface properties
        if (t >= avg_start)
            @timeit "avg physprops" avg_props!(phys_props_avg, phys_props, n_avg)
            @timeit "reduce surf props" reduce_surf_props!(surf_props_reduced, surf_props_chunks)
            @timeit "avg surfprops" avg_props!(surf_props_avg, surf_props_reduced, n_avg)
        end
    end

    @timeit "I/O" write_netcdf(ds_avg, phys_props_avg, n_timesteps)
    @timeit "I/O" write_netcdf(ds_surf_avg, surf_props_avg, n_timesteps)

    close_netcdf(ds_avg)
    close_netcdf(ds_surf_avg)

    print_timer()
end

run(1234, 300.0, 500.0, 5e-4, 5e22, 500, 500, 2.59e-9, 50000, 14000; chunk_count_multiplier=1, preallocation_margin_multiplier=1.5)
```