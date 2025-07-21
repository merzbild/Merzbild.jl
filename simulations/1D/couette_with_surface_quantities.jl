include("../../src/Merzbild.jl")

using ..Merzbild
using Random
using TimerOutputs

function run(seed, T_wall, v_wall, L, ndens, nx, ppc, Δt, output_freq, n_timesteps, avg_start)
    reset_timer!()

    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    # load particle and interaction data
    particles_data_path = joinpath("data", "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")
    interaction_data_path = joinpath("data", "vhs.toml")
    interaction_data::Array{Interaction, 2} = load_interaction_data(interaction_data_path, species_data)

    # create our grid and BCs
    grid = Grid1DUniform(L, nx)
    boundaries = MaxwellWalls1D(species_data, T_wall, T_wall, -v_wall, v_wall, 1.0, 1.0)

    # init particle vector, particle indexer, grid particle sorter
    n_particles = ppc * nx
    particles = [ParticleVector(n_particles)]
    pia = ParticleIndexerArray(grid.n_cells, 1)
    gridsorter = GridSortInPlace(grid, n_particles)

    # sample particles
    # Fnum * ppc = Np in cell = ndens * V_cell
    Fnum = grid.cells[1].V * ndens / ppc
    @timeit "sampling" sample_particles_equal_weight!(rng, grid, particles[1], pia, 1,
                                                      species_data, ndens, T_wall, Fnum)

    # create collision structs
    collision_data = CollisionData()
    
    # create struct for computation of physical properties
    phys_props = PhysProps(pia)

    # create second struct for averaging of physical properties
    phys_props_avg = PhysProps(pia)

    # create struct for computation of surface properties
    surf_props = SurfProps(pia, grid)

    # create second struct for averaging of physical properties
    surf_props_avg = SurfProps(pia, grid)

    # create struct for netCDF output
    ds = NCDataHolder("scratch/data/couette_$(L)_$(nx)_$(v_wall)_$(T_wall)_$(ppc).nc", species_data, phys_props)

    # create struct for second netCDF, this one is for time-averaged 
    ds_avg = NCDataHolder("scratch/data/avg_couette_$(L)_$(nx)_$(v_wall)_$(T_wall)_$(ppc)_after$(avg_start).nc",
                          species_data, phys_props)

    # create struct for time-averaged surface properties I/O
    ds_surf_avg = NCDataHolderSurf("scratch/data/avg_couette_$(L)_$(nx)_$(v_wall)_$(T_wall)_$(ppc)_surf_after$(avg_start).nc",
                                   species_data, surf_props_avg)

    # create and estime collision factors
    collision_factors = create_collision_factors_array(pia, interaction_data, species_data, T_wall, Fnum)

    # compute and write data at t=0
    compute_props!(particles, pia, species_data, phys_props)
    write_netcdf_phys_props(ds, phys_props, 0)
    write_grid("scratch/data/couette_$(L)_$(nx)_grid.nc", grid)

    n_avg = n_timesteps - avg_start + 1

    for t in 1:n_timesteps
        if t % 1000 == 0
            println(t)
        end

        # collide particles
        for cell in 1:grid.n_cells
            @timeit "collide" ntc!(rng, collision_factors[1, 1, cell],
                                   collision_data, interaction_data, particles[1], pia, cell, 1, Δt, grid.cells[cell].V)
        end

        # convect particles
        if (t < avg_start)
            @timeit "convect" convect_particles!(rng, grid, boundaries, particles[1], pia, 1, species_data, Δt)
        else
            @timeit "convect + surface compute" convect_particles!(rng, grid, boundaries, particles[1], pia, 1, species_data, surf_props, Δt)
            @timeit "avg surfprops" avg_props!(surf_props_avg, surf_props, n_avg)
        end

        # sort particles
        @timeit "sort" sort_particles!(gridsorter, grid, particles[1], pia, 1)

        # compute props and do I/O

        if (t < avg_start)
            if (t % output_freq == 0)
                @timeit "props compute" compute_props_sorted!(particles, pia, species_data, phys_props)
            end
        else
            @timeit "props compute" compute_props_sorted!(particles, pia, species_data, phys_props)
            @timeit "avg physprops" avg_props!(phys_props_avg, phys_props, n_avg)
        end


        if (t % output_freq == 0)
            @timeit "I/O" write_netcdf_phys_props(ds, phys_props, t)
        end
    end

    @timeit "I/O" write_netcdf_phys_props(ds_avg, phys_props_avg, n_timesteps)
    @timeit "I/O" write_netcdf_surf_props(ds_surf_avg, surf_props_avg, n_timesteps)

    close_netcdf(ds)
    close_netcdf(ds_avg)
    close_netcdf(ds_surf_avg)

    print_timer()
end

run(1234, 300.0, 500.0, 5e-4, 5e22, 50, 1000, 2.59e-9, 1000, 30000, 14000)