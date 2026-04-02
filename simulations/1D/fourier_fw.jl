# include("../../src/Merzbild.jl")

using Merzbild
using Random
using TimerOutputs

"""
    run(seed, T_bg0, T_wall1, T_wall2, v_wall, L, p0, nx,
        ppc_sampled,
        Δt, n_timesteps, avg_start,
        name)

Run a Fourier flow simulation with fixed-weight particles. Output is time-averaged and written to file
with name `avg_[name]_ntc_[L]_[nx]_[v_wall]_[T_wall1]_[T_wall2]_after[avg_start].nc` (fluxes
and surface properties following similar naming scheme) in `scratch/data/`.

Positional arguments:
* `seed`: random seed value
* `T_bg0`: initial temperature of gas
* `T_wall1`: temperature of left wall
* `T_wall2`: temperature of right wall
* `v_wall`: left wall will have y-velocity of `-v_wall`, right wall will have y-velocity of `v_wall` (set to 0 for Fourier
simulation)
* `L`: length of domain in meters
* `p0`: initial gas pressure in Pa
* `nx`: number of grid cells
* `ppc_sampled`: number of particles sampled in each cell
* `Δt`: timestep size
* `n_timesteps`: total number of timesteps to run for
* `avg_start`: timestep after which averaging begins
* `name`: name to append to start of output files
"""
function run(seed, T_bg0, T_wall1, T_wall2, v_wall, L, p0, nx,
             ppc_sampled,
             Δt, n_timesteps, avg_start,
             name)
    reset_timer!()

    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    # load particle and interaction data
    particles_data_path = joinpath("data", "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")
    interaction_data_path = joinpath("data", "vhs.toml")
    interaction_data::Array{Interaction, 2} = load_interaction_data(interaction_data_path, species_data)


    ndens = p0 / (k_B * T_bg0)
    println("n = $ndens")
    nu = mean_collision_frequency(interaction_data, species_data, 1, ndens, 0.5*(T_wall1 + T_wall2))
    lam = mean_free_path(interaction_data, 1, ndens, 0.5*(T_wall1 + T_wall2))
    println("1/ν = $(1/nu)")
    println("n timesteps before avg: $(100 * (1/nu) / Δt)")
    println("λ = $lam; Kn = $(lam / L); λ/dx = $(lam / (L / nx))")
    vr = sqrt(2 * k_B * T_bg0 / species_data[1].mass)

    println("CFL time: $(L / nx / vr)")
    # return
    # create our grid and BCs
    grid = Grid1DUniform(L, nx)
    boundaries = MaxwellWalls1D(species_data, T_wall1, T_wall2, -v_wall, v_wall, 1.0, 1.0)

    # init particle vector, particle indexer, grid particle sorter
    n_particles = ppc_sampled * nx
    particles = [ParticleVector(n_particles)]
    pia = ParticleIndexerArray(grid.n_cells, 1)
    gridsorter = GridSortInPlace(grid, n_particles)

    # sample particles
    # Fnum * ppc_sampled = Np in cell = ndens * V_cell
    Fnum = grid.cells[1].V * ndens / ppc_sampled
    @timeit "sampling" sample_particles_equal_weight!(rng, grid, particles[1], pia, 1,
                                                      species_data, ndens, T_bg0, Fnum)

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

    flux_props = FluxProps(pia)
    flux_props_avg = FluxProps(pia)

    # create struct for netCDF output
    # ds = NCDataHolder("scratch/data/couette_swpm_$(L)_$(nx)_$(v_wall)_$(T_wall)_$(merge_threshold)_$(merge_target).nc", species_data, phys_props)

    # create struct for second netCDF, this one is for time-averaged 
    ds_avg = NCDataHolder("scratch/data/avg_$(name)_ntc_$(L)_$(nx)_$(v_wall)_$(T_wall1)_$(T_wall2)_after$(avg_start).nc",
                          species_data, phys_props)

    # create struct for time-averaged surface properties I/O
    ds_surf_avg = NCDataHolderSurf("scratch/data/avg_$(name)_ntc_$(L)_$(nx)_$(v_wall)_$(T_wall1)_$(T_wall2)_surf_after$(avg_start).nc",
                                   species_data, surf_props_avg)

    # create struct for time-averaged fluxal properties I/O
    ds_flux_avg = NCDataHolderFlux("scratch/data/avg_$(name)_ntc_$(L)_$(nx)_$(v_wall)_$(T_wall1)_$(T_wall2)_flux_after$(avg_start).nc",
                                   species_data, flux_props_avg)

    # create and estimate collision factors
    # because we do a merge at the start need to account for changed average Fnum
    collision_factors = create_collision_factors_array(pia, interaction_data, species_data, T_bg0, Fnum)


    write_grid("scratch/data/couette_$(L)_$(nx)_grid.nc", grid)

    n_avg = n_timesteps - avg_start + 1

    println(grid.cells[1].V, " ", grid.cells[2].V)
    for t in 1:n_timesteps
        if t % 1000 == 0
            println(t)
        end

        for cell in 1:grid.n_cells
            @timeit "collide" ntc!(rng, collision_factors[1, 1, cell],
                                   collision_data, interaction_data, particles[1], pia, cell, 1, Δt, grid.cells[cell].V)
        end

        # convect particles
        if (t < avg_start)
            @timeit "convect" convect_particles!(rng, grid, boundaries, particles[1], pia, 1, species_data, Δt)
        else
            @timeit "convect + surface compute" convect_particles!(rng, grid, boundaries, particles[1], pia, 1, species_data, surf_props, Δt)
            avg_props!(surf_props_avg, surf_props, n_avg)
        end

        # sort particles
        @timeit "sort" sort_particles!(gridsorter, grid, particles[1], pia, 1)

        if (t >= avg_start)
            @timeit "props compute" compute_props_sorted!(particles, pia, species_data, phys_props)
            @timeit "flux props compute" compute_flux_props_sorted!(particles, pia, species_data, phys_props, flux_props, grid)
            avg_props!(phys_props_avg, phys_props, n_avg)
            avg_props!(flux_props_avg, flux_props, n_avg)
        end
    end

    @timeit "I/O" write_netcdf(ds_avg, phys_props_avg, n_timesteps)
    @timeit "I/O" write_netcdf(ds_surf_avg, surf_props_avg, n_timesteps)
    @timeit "I/O" write_netcdf(ds_flux_avg, flux_props_avg, n_timesteps)

    close_netcdf(ds_avg)
    close_netcdf(ds_surf_avg)
    close_netcdf(ds_flux_avg)

    print_timer()
end

const n_t = 3500000
run(1234, 350.0, 300.0, 600.0, 0.0, 5e-3, 1.5e-3 * 101325.0, 1000, 500, 4e-10, n_t, 500_000, "Fourier_FW")
