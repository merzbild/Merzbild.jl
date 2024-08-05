@testset "test I/O on grids" begin
    seed = 1234
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_list(particles_data_path, "Ar")

    grid = create_grid1D_uniform(4, 8)
    ppc = 1000
    T = 400.0
    phys_props::PhysProps = create_props(grid.n_cells, 1, [], Tref=1)

    particles = [create_particle_vector(ppc * grid.n_cells)]
    pia = create_particle_indexer_array(grid.n_cells, 1)

    ndens = 1e23
    Fnum = 1e20
    sample_particles_equal_weight!(rng, grid, particles[1], pia, 1,
                                   species_data, ndens, T, Fnum)

    compute_props!(particles, pia, species_data, phys_props)

    sol_path = joinpath(@__DIR__, "data", "tmp_1dgrid.nc")
    ds = create_netcdf_phys_props(sol_path, species_data, phys_props)
    write_netcdf_phys_props(ds, phys_props, 0)
    close_netcdf(ds)

    sol = NCDataset(sol_path, "r")

    @test length(sol["timestep"]) == 1

    @test size(sol["np"]) == (grid.n_cells, 1, 1)
    @test size(sol["T"]) == (grid.n_cells, 1, 1)
    @test size(sol["v"]) == (3, grid.n_cells, 1, 1)
    @test size(sol["ndens"]) == (grid.n_cells, 1, 1)

    @test maximum(abs.(sol["np"][:, 1, 1] - phys_props.np[:, 1])) == 0
    @test maximum(abs.(sol["ndens"][:, 1, 1] - phys_props.n[:, 1])) < 2 * eps()
    @test maximum(abs.(sol["T"][:, 1, 1] - phys_props.T[:, 1])) < 2 * eps()
    @test maximum(abs.(sol["v"][1, :, 1, 1] - phys_props.v[1, :, 1])) < 2 * eps()
    @test maximum(abs.(sol["v"][2, :, 1, 1] - phys_props.v[2, :, 1])) < 2 * eps()
    @test maximum(abs.(sol["v"][3, :, 1, 1] - phys_props.v[3, :, 1])) < 2 * eps()

    close(sol)
    rm(sol_path)
end