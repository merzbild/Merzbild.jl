# include("../../src/Merzbild.jl")

using Merzbild
using Random
using TimerOutputs
using Base.Threads
using ChunkSplitters

function run(seed, T_wall, v_wall, L, ndens, nx, ppc_sampled, merge_threshold, merge_target, Δt, output_freq, n_timesteps, avg_start;
    chunk_count_multiplier=1, preallocation_margin_multiplier=1.0, parallel_exchange=true, final_debug=false)
    reset_timer!()
    main_to = TimerOutput()

    n_threads = Threads.nthreads()
    n_chunks = n_threads * chunk_count_multiplier
    println("Running on $n_threads threads, will split cells into $n_chunks chunks")

    rng_chunks = [Xoshiro(seed + i) for i in 0:n_chunks-1]

    # load particle and interaction data
    particles_data_path = joinpath(MERZBILD_DATA_PATH, "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")
    interaction_data_path = joinpath(MERZBILD_DATA_PATH, "vhs.toml")
    interaction_data::Array{Interaction, 2} = load_interaction_data(interaction_data_path, species_data)

    # create our grid and BCs
    grid = Grid1DUniform(L, nx)
    bc_list = (FullyDiffuseBC1D(1, species_data, T_wall, [0.0, -v_wall, 0.0]),
               FullyDiffuseBC1D(1, species_data, T_wall, [0.0, v_wall, 0.0]))

    # split cell indices into chunks
    cell_indices = Vector(1:nx)
    cell_chunks = chunks(cell_indices; n=n_chunks)

    # init per-chunk particle vectors, particle indexers, grid particle sorters
    n_particles_chunks = [floor(Int64, ppc_sampled * length(cell_chunk) * preallocation_margin_multiplier) for cell_chunk in cell_chunks]
    particles_chunks = [[ParticleVector{1}(n_particles)] for n_particles in n_particles_chunks]
    pia_chunks = [ParticleIndexerArray(grid.n_cells, 1) for cell_chunk in cell_chunks]
    gridsorter_chunks = [GridSortInPlace(grid, n_particles) for n_particles in n_particles_chunks]

    # this is used for moving particles between chunks after they have been sorted into grid cells
    chunk_exchanger = ChunkExchanger(cell_chunks, nx)

    # sample particles
    # Fnum * ppc_sampled = Np in cell = ndens * V_cell
    Fnum = grid.cells[1].V * ndens / ppc_sampled

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

    # create struct for netCDF for physical properties I/O
    ds = NCDataHolder("scratch/data/mt$(n_threads)_$(n_chunks)ch_couette_$(L)_$(nx)_$(v_wall)_$(T_wall)_$(ppc_sampled)_$(merge_threshold)_$(merge_target)_octree_.nc",
                      species_data, phys_props)

    # create struct for netCDF for time-averaged physical properties I/O
    ds_avg = NCDataHolder("scratch/data/avg_mt$(n_threads)_$(n_chunks)ch_couette_$(L)_$(nx)_$(v_wall)_$(T_wall)_$(ppc_sampled)_$(merge_threshold)_$(merge_target)_octree_after$(avg_start).nc",
                          species_data, phys_props)

    # create struct for netCDF for time-averaged surface properties I/O
    ds_surf_avg = NCDataHolderSurf("scratch/data/avg_mt$(n_threads)_$(n_chunks)ch_couette_$(L)_$(nx)_$(v_wall)_$(T_wall)_$(ppc_sampled)_$(merge_threshold)_$(merge_target)_octree_surf_after$(avg_start).nc",
                                   species_data, surf_props_avg)

    # create merging structs
    oc_chunks = [OctreeMerge{1,2}(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel, max_Nbins=6000) for cell_chunks in cell_chunks]

    # merge and compute data at t=0
    @threads for chunk_id in 1:n_chunks
        for cell in cell_chunks[chunk_id]
            merge_octree!(rng_chunks[chunk_id], oc_chunks[chunk_id], particles_chunks[chunk_id][1], pia_chunks[chunk_id], cell, 1, merge_target, grid)
        end
        squash_pia!(particles_chunks[chunk_id], pia_chunks[chunk_id])
        compute_props_sorted!(particles_chunks[chunk_id], pia_chunks[chunk_id], species_data, phys_props, cell_chunks[chunk_id])
    end

    # create and estimate collision factors
    # we correct Fnum estimate since we performed merging
    collision_factors = [create_collision_factors_array(pia, interaction_data, species_data, T_wall, Fnum * ppc_sampled / merge_target)
                         for pia in pia_chunks]

    # compute data at t=0
    @threads for (chunk_id, cell_chunk) in enumerate(cell_chunks)
        compute_props_sorted!(particles_chunks[chunk_id], pia_chunks[chunk_id], species_data, phys_props, cell_chunk)
    end

    index_inv_map = [zeros(Int64, floor(Int64, merge_threshold * length(cell_chunk) * preallocation_margin_multiplier)) for cell_chunk in cell_chunks]

    n_avg = n_timesteps - avg_start + 1

    chunk_timers = [TimerOutput() for _ in 1:n_chunks]

    exchange_list = generate_1_factorization(n_chunks)
    println("exchange_list = $(exchange_list)")

    @timeit main_to "main loop" @inbounds for t in 1:n_timesteps
        if t % 1000 == 0
            println(t)
        end
        
        # collide, convect, sort particles
        @threads for chunk_id in 1:n_chunks
            @inbounds rng_local = rng_chunks[chunk_id]
            @inbounds cells_local = cell_chunks[chunk_id]
            @inbounds collision_factors_local = collision_factors[chunk_id]
            @inbounds coll_data = collision_data[chunk_id]
            @inbounds pia_local = pia_chunks[chunk_id]
            @inbounds particles_local = particles_chunks[chunk_id]
            @inbounds local_timer = chunk_timers[chunk_id]
            
            for cell in cells_local
                @timeit local_timer "collide (t)" ntc!(rng_local, collision_factors_local[1, 1, cell],
                               coll_data, interaction_data, particles_local[1],
                               pia_local, cell, 1, Δt, grid.cells[cell].V)

                if pia_local.indexer[cell,1].n_local > merge_threshold
                    @timeit local_timer "merge (t)" merge_octree!(rng_local, oc_chunks[chunk_id], particles_local[1], pia_local, cell, 1, merge_target, grid)
                end
            end

            @timeit local_timer "squash_pia (t)" squash_pia!(particles_local, pia_local)

            if (t >= avg_start)
                @timeit local_timer "convect + surface compute (t)" @inbounds convect_particles!(rng_local, grid, bc_list,
                                    particles_local[1], pia_local,
                                    1, species_data, surf_props_chunks[chunk_id], Δt)
            else
                # we don't need to compute surface properties before we start averaging
                @timeit local_timer "convect (t)" @inbounds convect_particles!(rng_local, grid, bc_list,
                                    particles_local[1], pia_local,
                                    1, species_data, Δt)
            end

            # need to clear the data in the chunk exchanger
            reset!(chunk_exchanger, chunk_id)

            # sort particles
            @timeit local_timer "sort (t)" @inbounds sort_particles!(gridsorter_chunks[chunk_id], grid, particles_chunks[chunk_id][1], pia_chunks[chunk_id], 1)

            if t%10 == 0
                @timeit local_timer "restore ordering (t)" restore_particle_ordering!(particles_chunks[chunk_id][1], index_inv_map[chunk_id])
            end
        end

        # move particles between chunks
        if parallel_exchange
            @timeit main_to "exchange (t)" for exchange_partial in exchange_list
                @threads for exchange_pair in exchange_partial
                    exchange_particles!(chunk_exchanger, particles_chunks, pia_chunks, cell_chunks, 
                                        1, exchange_pair[1], exchange_pair[2])
                end
            end
        else
            @timeit main_to "exchange" exchange_particles!(chunk_exchanger, particles_chunks, pia_chunks, cell_chunks, 1)
        end

        # reset indexing, compute physical properties if needed
        @threads for chunk_id in 1:n_chunks
            @inbounds cells_local = cell_chunks[chunk_id]
            @inbounds particles_local = particles_chunks[chunk_id]
            @inbounds pia_local = pia_chunks[chunk_id]
            @inbounds local_timer = chunk_timers[chunk_id]

            @timeit local_timer "sort post-exchange (t)" sort_particles_after_exchange!(chunk_exchanger, gridsorter_chunks[chunk_id],
                                           particles_local[1], pia_local,
                                           cells_local, 1)
            if (t >= avg_start)
                @timeit local_timer "props compute (t)" compute_props_sorted!(particles_local, pia_local,
                                                species_data, phys_props, cells_local)
            elseif (t % output_freq == 0)
                @timeit local_timer "props compute (t)" compute_props_sorted!(particles_local, pia_local,
                                                species_data, phys_props, cells_local)
            end
        end

        if (t % output_freq == 0)
            @timeit main_to "I/O" write_netcdf(ds, phys_props, t)
        end

        # reduce surface properties, average grid and surface properties
        if (t >= avg_start)
            @timeit main_to "avg physprops" avg_props!(phys_props_avg, phys_props, n_avg)
            @timeit main_to "reduce surf props" reduce_surf_props!(surf_props_reduced, surf_props_chunks)
            @timeit main_to "avg surfprops" avg_props!(surf_props_avg, surf_props_reduced, n_avg)
        end
    end

    @timeit main_to "I/O final" write_netcdf(ds_avg, phys_props_avg, n_timesteps)
    @timeit main_to "I/O final" write_netcdf(ds_surf_avg, surf_props_avg, n_timesteps)

    close_netcdf(ds)
    close_netcdf(ds_avg)
    close_netcdf(ds_surf_avg)

    if final_debug
        println("ndens = $(sum(phys_props.n))")
        println("Np(total) = $(sum(phys_props.np))")
        println([check_pia_is_correct(pia_chunks[chunk_id], 1) for chunk_id in 1:n_chunks])
        println([check_unique_buffer(particles_chunks[chunk_id][1]) for chunk_id in 1:n_chunks])
        println([check_unique_index(particles_chunks[chunk_id][1], pia_chunks[chunk_id], 1) for chunk_id in 1:n_chunks])
    end

    for t in chunk_timers
        merge!(main_to, t, tree_point=["main loop"])
    end

    # we now fix ncall, timing, allocation counts for the threaded timers by hand
    # (t) denotes multi-threaded part
    timers_to_average = ["collide (t)", "merge (t)", "squash_pia (t)", "sort (t)", "restore ordering (t)", "convect (t)",
                         "convect + surface compute (t)", "sort post-exchange (t)", "props compute (t)"]

    for timer_name in timers_to_average
        try
            accd = main_to["main loop"].inner_timers[timer_name].accumulated_data
            accd.ncalls = round(Int64, accd.ncalls/n_chunks)
            accd.time = round(Int64, accd.time/n_chunks)
            accd.allocs = round(Int64, accd.allocs/n_chunks)
        catch
            nothing
        end
    end

    print_timer(main_to)
end

const n_t = 50000
run(1234, 300.0, 500.0, 5e-4, 5e22, 500, 500, 130, 100, 2.59e-9, 1000, n_t, 14000;
    chunk_count_multiplier=1, preallocation_margin_multiplier=1.5, parallel_exchange=true, final_debug=true)