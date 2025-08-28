include("../../src/Merzbild.jl")

using ..Merzbild
using Random
using TimerOutputs
using Base.Threads
using ChunkSplitters

function run(seed, T_wall, v_wall, L, ndens, nx, ppc, Δt, output_freq, n_timesteps, avg_start; chunk_count_multiplier=1,
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

    # create struct for netCDF for physical properties I/O
    ds = NCDataHolder("scratch/data/mt$(n_threads)_$(n_chunks)ch_couette_$(L)_$(nx)_$(v_wall)_$(T_wall)_$(ppc).nc",
                      species_data, phys_props)

    # create struct for netCDF for time-averaged physical properties I/O
    ds_avg = NCDataHolder("scratch/data/avg_mt$(n_threads)_$(n_chunks)ch_couette_$(L)_$(nx)_$(v_wall)_$(T_wall)_$(ppc)_after$(avg_start).nc",
                          species_data, phys_props)

    # create struct for netCDF for time-averaged surface properties I/O
    ds_surf_avg = NCDataHolderSurf("scratch/data/avg_mt$(n_threads)_$(n_chunks)ch_couette_$(L)_$(nx)_$(v_wall)_$(T_wall)_$(ppc)_surf_after$(avg_start).nc",
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
        if t % 1000 == 0
            println(t)
        end

        # if t % 10000 == 0
        #     println(t, ", ", sum([pia_chunks[chunk_id].n_total[1] for chunk_id in 1:n_chunks]))
        #     println([pia_chunks[chunk_id].n_total[1] for chunk_id in 1:n_chunks])
        #     println([check_pia_is_correct(pia_chunks[chunk_id], 1) for chunk_id in 1:n_chunks])
        #     println([maximum([pia_chunks[chunk_id].indexer[cell, 1].n_local for cell in cell_chunks[chunk_id]])
        #              for chunk_id in 1:n_chunks])
        #     println()
        # end
        
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
            elseif (t % output_freq == 0)
                @inbounds compute_props_sorted!(particles_chunks[chunk_id], pia_chunks[chunk_id],
                                                species_data, phys_props, cell_chunks[chunk_id])
            end
        end

        if (t % output_freq == 0)
            @timeit "I/O" write_netcdf(ds, phys_props, t)
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

    close_netcdf(ds)
    close_netcdf(ds_avg)
    close_netcdf(ds_surf_avg)

    print_timer()
end

run(1234, 300.0, 500.0, 5e-4, 5e22, 2000, 250, 2.59e-9, 1000, 50000, 14000; chunk_count_multiplier=1, preallocation_margin_multiplier=1.5)