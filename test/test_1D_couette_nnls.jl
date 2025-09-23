
@testset "couette test variable weight NNLS" begin

    # NNLS merging seems to be quite sensitive to actual architecture/OS
    # so we do only a few timesteps otherwise reference solution may diverge
    # from the one obtained on a different machine
    T_wall = 300.0
    v_wall = 500.0
    L = 5e-4
    ndens = 5e22
    nx = 50
    ppc = 1000
    Δt = 2.59e-9
    output_freq = 1
    n_timesteps = 20
    merge_threshold = 70
    merge_target = 60
    n_cons_moments = 5

    seed = 1234
    rng = StableRNG(seed)

    # load particle and interaction data
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")
    interaction_data_path = joinpath(@__DIR__, "..", "data", "vhs.toml")
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

    sample_particles_equal_weight!(rng, grid, particles[1], pia, 1,
                                   species_data, ndens, T_wall, Fnum)

    # create collision structs
    collision_data = CollisionData()
    
    # create struct for computation of physical properties
    phys_props = PhysProps(pia)


    # create struct for netCDF output
    sol_path = joinpath(@__DIR__, "data", "tmp_couette_nnls.nc")
    ds = NCDataHolder(sol_path, species_data, phys_props)

    # create struct for computation of surface properties
    surf_props = SurfProps(pia, grid)

    sol_path_surf = joinpath(@__DIR__, "data", "tmp_couette_nnls_surf.nc")
    ds_surf = NCDataHolderSurf(sol_path_surf, species_data, surf_props)

    # init collision factors
    collision_factors = create_collision_factors_array(pia, interaction_data, species_data, T_wall, Fnum)

    oc = OctreeN2Merge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel, max_Nbins=6000)


    mim = []
    n_moms = n_cons_moments
    for i in 1:n_moms
        append!(mim, compute_multi_index_moments(i))
    end

    # conserve center of mass and variance in x direction
    pos_moments = [[1,0,0],[2,0,0]]

    init_np = merge_threshold
    matrix_ncol_nprealloc = 15

    # the main NNLS
    mnnls = NNLSMerge(mim, init_np; multi_index_moments_pos=pos_moments, matrix_ncol_nprealloc=matrix_ncol_nprealloc)

    fails = 0

    for cell in 1:grid.n_cells
        if pia.indexer[cell,1].n_local > merge_threshold

            nnls_success_flag = merge_nnls_based!(rng, mnnls, particles[1], pia, cell, 1; centered_at_mean=false, v_multipliers=[], w_threshold=1e-12)

            if nnls_success_flag == -1
                fails += 1
                merge_octree_N2_based!(rng, oc, particles[1], pia, cell, 1, merge_target, grid)
            end
            # merge_octree_N2_based!(rng, oc, particles[1], pia, cell, 1, merge_target, grid)
        end
    end
    squash_pia!(particles, pia)

    # compute and write data at t=0
    compute_props!(particles, pia, species_data, phys_props)
    write_netcdf(ds, phys_props, 0)

    # n_avg = n_timesteps - avg_start + 1

    for t in 1:n_timesteps
        # collide particles
        for cell in 1:grid.n_cells
            ntc!(rng, collision_factors[1, 1, cell],
                 collision_data, interaction_data, particles[1], pia, cell, 1, Δt, grid.cells[cell].V)
        
            if pia.indexer[cell,1].n_local > merge_threshold
                nnls_success_flag = merge_nnls_based!(rng, mnnls, particles[1], pia, cell, 1; centered_at_mean=false, v_multipliers=[], w_threshold=1e-12)
    
                if nnls_success_flag == -1
                    merge_octree_N2_based!(rng, oc, particles[1], pia, cell, 1, merge_target, grid)
                end
                squash_pia!(particles, pia)
            end
        end

        # convect particles
        convect_particles!(rng, grid, boundaries, particles[1], pia, 1, species_data, surf_props, Δt)

        # sort particles
        sort_particles!(gridsorter, grid, particles[1], pia, 1)

        # if (t < avg_start)
        if (t % output_freq == 0)
            compute_props_sorted!(particles, pia, species_data, phys_props)
        end

        if (t % output_freq == 0)
            write_netcdf(ds, phys_props, t)
            write_netcdf(ds_surf, surf_props, t)
        end
    end

    close_netcdf(ds)
    close_netcdf(ds_surf)

    ref_sol_path = joinpath(@__DIR__, "data", "couette_0.0005_50_500.0_300.0_1000_nnls5moments.nc")
    ref_sol = NCDataset(ref_sol_path, "r")
    sol = NCDataset(sol_path, "r")
 
    ndens_conservation = true
    ref_ndens = 2.5e19

    tmax_test = 20

    for t in 1:tmax_test
        if abs(sum(sol["ndens"][:, 1, t]) - ref_ndens) / ref_ndens > 1e-15
            ndens_conservation = false
        end
    end

    @test fails == 0
    @test ndens_conservation == true
    @test maximum(abs.(ref_sol["ndens"][:, 1, 1:3] .- sol["ndens"][:, 1, 1:3])) < 2 * eps()
    @test maximum(abs.(ref_sol["T"][:, 1, 1:3] .- sol["T"][:, 1, 1:3])) < 2.4e-13

    close(sol)
    rm(sol_path)


    ref_sol_path = joinpath(@__DIR__, "data", "couette_0.0005_50_500.0_300.0_1000_nnls5moments_surf.nc")
    ref_sol = NCDataset(ref_sol_path, "r")
    sol = NCDataset(sol_path_surf, "r")
 
    @test maximum(abs.(ref_sol["np"][:, 1, 1:tmax_test] .- sol["np"][:, 1, 1:tmax_test])) < 2 * eps()
    @test maximum(abs.(ref_sol["kinetic_energy_flux"][:, 1, 1:tmax_test]
                       .- sol["kinetic_energy_flux"][:, 1, 1:tmax_test])) < 6.0e-11

    @test maximum(abs.(sol["flux_incident"][:, 1, 1:tmax_test]
                       - ref_sol["flux_incident"][:, 1, 1:tmax_test])) < 2 * eps()

    @test maximum(abs.(sol["flux_incident"][:, 1, 1:tmax_test]
                       + sol["flux_reflected"][:, 1, 1:tmax_test])) < 2 * eps()

    close(sol)
    rm(sol_path_surf)
end
