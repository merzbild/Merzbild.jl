include("../../src/Merzbild.jl")

using ..Merzbild
using Random
using TimerOutputs

function run(seed, T_bg0, T_wall1, T_wall2, v_wall, L, p0, nx,
             ppc_sampled, merge_threshold, merge_target, Δt, n_timesteps, avg_start,
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

    # create struct for second netCDF, this one is for time-averaged 
    ds_avg = NCDataHolder("scratch/data/avg_$(name)_ntc_$(L)_$(nx)_$(v_wall)_$(T_wall1)_$(T_wall2)_$(merge_threshold)_$(merge_target)_after$(avg_start).nc",
                          species_data, phys_props)

    # create struct for time-averaged surface properties I/O
    ds_surf_avg = NCDataHolderSurf("scratch/data/avg_$(name)_ntc_$(L)_$(nx)_$(v_wall)_$(T_wall1)_$(T_wall2)_$(merge_threshold)_$(merge_target)_surf_after$(avg_start).nc",
                                   species_data, surf_props_avg)

    # create struct for time-averaged fluxal properties I/O
    ds_flux_avg = NCDataHolderFlux("scratch/data/avg_$(name)_ntc_$(L)_$(nx)_$(v_wall)_$(T_wall1)_$(T_wall2)_$(merge_threshold)_$(merge_target)_flux_after$(avg_start).nc",
                                   species_data, flux_props_avg)

    # create and estimate collision factors
    # because we do a merge at the start need to account for changed average Fnum
    collision_factors = create_collision_factors_array(pia, interaction_data, species_data, T_bg0, Fnum * (merge_threshold / merge_target))

    oc = OctreeN2Merge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel, max_Nbins=6000)

    for cell in 1:grid.n_cells
        if pia.indexer[cell,1].n_local > merge_threshold
            @timeit "merge (t=0)" merge_octree_N2_based!(rng, oc, particles[1], pia, cell, 1, merge_target, grid)
        end
    end
    @timeit "squash (t=0)" squash_pia!(particles, pia)

    write_grid("scratch/data/couette_$(L)_$(nx)_grid.nc", grid)

    n_avg = n_timesteps - avg_start + 1

    for t in 1:n_timesteps
        if t % 1000 == 0
            println(t)
        end

        ncolls = 0
        for cell in 1:grid.n_cells
            @timeit "collide" ntc!(rng, collision_factors[1, 1, cell],
                                   collision_data, interaction_data, particles[1], pia, cell, 1, Δt, grid.cells[cell].V)

            ncolls += collision_factors[1, 1, cell].n_coll_performed
            if pia.indexer[cell,1].n_local > merge_threshold
                @timeit "merge" merge_octree_N2_based!(rng, oc, particles[1], pia, cell, 1, merge_target, grid)
                @timeit "squash" squash_pia!(particles, pia)
            end
        end
        # println("ncolls = $(ncolls)")

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

        # if (t % 5000 == 0)
        #     compute_props_sorted!(particles, pia, species_data, phys_props)
        #     # println("n = ", sum(phys_props.n))
        #     println("T = ", sum(phys_props.T) / 5000.0)
        # end

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
# const params = [[45, 38], [69, 58], [105, 88], [146, 122], [200, 166], [266, 220]] - iterate over all N:M merges
const params = [[45, 38]]

for param in params
    run(1234, 350.0, 300.0, 600.0, 0.0, 5e-3, 1.5e-3 * 101325.0, 1000, 500, param[1], param[2], 4e-10, n_t, 500_000, "Fourier_octree")
end
