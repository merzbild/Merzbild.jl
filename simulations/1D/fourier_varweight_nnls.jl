# include("../../src/Merzbild.jl")

using Merzbild
using Random
using TimerOutputs

"""
    run(seed, T_bg0, T_wall1, T_wall2, v_wall, L, p0, nx,
        ppc_sampled, merge_threshold, merge_target, n_vel_up_total,
        n_moms_pos,
        Δt, n_timesteps, avg_start,
        name)

Run a Fourier flow simulation with variable-weight particles and NNLS merging. Output is time-averaged and written to file
with name `avg_[name]_ntc_[L]_[nx]_[v_wall]_[T_wall1]_[T_wall2]_[merge_threshold]_[merge_target]_after[avg_start].nc` (fluxes
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
* `merge_threshold`: threshold number of particles in a cell after which they are merged
* `merge_target`: target number of particles to merge down to in case octree merging is used as backup merging
* `n_vel_up_total`: all mixed velocity moments up to this order are conserved during merging
* `n_moms_pos`: all spatial moments in the x-direction up to this order are conserved during merging
* `Δt`: timestep size
* `n_timesteps`: total number of timesteps to run for
* `avg_start`: timestep after which averaging begins
* `name`: name to append to start of output files
"""
function run(seed, T_bg0, T_wall1, T_wall2, v_wall, L, p0, nx,
             ppc_sampled, merge_threshold, merge_target, n_vel_up_total,
             n_moms_pos,
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


    mim = []
    n_moms = n_vel_up_total
    for i in 1:n_moms
        append!(mim, compute_multi_index_moments(i))
    end

    # conserve center of mass and variance in x direction
    pos_moments = [[1,0,0],[2,0,0]]
    for i in 3:n_moms_pos
        append!(pos_moments, [[i,0,0]])
    end

    println(pos_moments)

    init_np = merge_threshold
    matrix_ncol_nprealloc = 25

    # the main NNLS
    @timeit "NNLSinit" mnnls = NNLSMerge(mim, init_np; multi_index_moments_pos=pos_moments, matrix_ncol_nprealloc=matrix_ncol_nprealloc)


    mim_backup = []
    n_moms_backup = n_vel_up_total - 1
    for i in 1:n_moms_backup
        append!(mim_backup, compute_multi_index_moments(i))
    end

    # the backup NNLS that conserves moments of total order 1 lower than the main one
    @timeit "NNLSinit backup" mnnls_backup = NNLSMerge(mim_backup, init_np; multi_index_moments_pos=pos_moments, matrix_ncol_nprealloc=matrix_ncol_nprealloc)

    for cell in 1:grid.n_cells
        if pia.indexer[cell,1].n_local > merge_threshold
            
            if pia.indexer[cell,1].n_local > merge_threshold
                @timeit "merge NNLS" nnls_success_flag = merge_nnls_based!(rng, mnnls, particles[1], pia, cell, 1;
                                                                           scaling=:variance, iteration_mult=5,
                                                                           centered_at_mean=false, v_multipliers=[], w_threshold=1e-12)
    
                if nnls_success_flag == -1
                    @timeit "merge NNLS backup" nnls_success_flag = merge_nnls_based!(rng, mnnls, particles[1], pia, cell, 1;
                                                                                      scaling=:variance, iteration_mult=5,
                                                                                      centered_at_mean=false, v_multipliers=[], w_threshold=1e-12)
                end
    
                if nnls_success_flag == -1
                    @timeit "merge octree" merge_octree_N2_based!(rng, oc, particles[1], pia, cell, 1, merge_target, grid)
                end
                @timeit "squash" squash_pia!(particles, pia)
            end
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
                @timeit "merge NNLS" nnls_success_flag = merge_nnls_based!(rng, mnnls, particles[1], pia, cell, 1;
                                                                           scaling=:variance, iteration_mult=5,
                                                                           centered_at_mean=false, v_multipliers=[], w_threshold=1e-12)
    
                if nnls_success_flag == -1
                    @timeit "merge NNLS backup" nnls_success_flag = merge_nnls_based!(rng, mnnls, particles[1], pia, cell, 1;
                                                                                      scaling=:variance, iteration_mult=5,
                                                                                      centered_at_mean=false, v_multipliers=[], w_threshold=1e-12)
                end
    
                if nnls_success_flag == -1
                    @timeit "merge octree" merge_octree_N2_based!(rng, oc, particles[1], pia, cell, 1, merge_target, grid)
                end
                @timeit "squash" squash_pia!(particles, pia)
            end
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


# Each elements of params is a list with 3 elements: number of velocity moments conserved, threshold number of particles,
# target number of particles for backup octree merging (in case octree fails)
# if NNLS with param[1] velocity moments fails, a backup NNLS is called with param[1]-1 moments, if that fails - octree is called
const params = [[4, 45, 38]]
# const params = [[4, 45, 38], [5, 69, 58], [6, 105, 88], [7, 146, 122], [8, 200, 166]] - iterate over all merging parameters
# uncomment line above to run over the parameter sets used for "Moment-preserving particle merging via non-negative least squares"

const x_moms = 2  # spatial moments of up to this order are conserved

const n_t = 3500000

for param in params
    println(x_moms, ", ", param)
    run(1234, 350.0, 300.0, 600.0, 0.0, 5e-3, 1.5e-3 * 101325.0, 1000, 500,
        param[2], param[3], param[1], x_moms, 4e-10, n_t, 500_000, "Fourier_NNLS")
end