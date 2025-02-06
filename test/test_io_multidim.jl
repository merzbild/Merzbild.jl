@testset "test I/O on grids" begin
    seed = 1234
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")

    grid = Grid1DUniform(4, 8)
    ppc = 1000
    T = 400.0
    phys_props::PhysProps = PhysProps(grid.n_cells, 1, [], Tref=1)

    particles = [ParticleVector(ppc * grid.n_cells)]
    pia = ParticleIndexerArray(grid.n_cells, 1)

    ndens = 1e23
    Fnum = 1e20
    sample_particles_equal_weight!(rng, grid, particles[1], pia, 1,
                                   species_data, ndens, T, Fnum)

    compute_props!(particles, pia, species_data, phys_props)

    sol_path = joinpath(@__DIR__, "data", "tmp_1dgrid.nc")
    ds = NCDataHolder(sol_path, species_data, phys_props)
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

    grid_path = joinpath(@__DIR__, "data", "tmp_grid_info.nc")
    write_grid(grid_path, grid)

    grid_info =  NCDataset(grid_path, "r")

    @test grid_info.attrib["L"] == grid.L
    @test grid_info.attrib["dim"] == 1
    @test grid_info.attrib["dx"] == grid.Δx
    @test grid_info.dim["cells"] == grid.n_cells
    @test maximum(abs.(grid_info["xlo"][:] - [grid.cells[i].xlo for i in 1:grid.n_cells])) < 2 * eps()
    @test maximum(abs.(grid_info["xhi"][:] - [grid.cells[i].xhi for i in 1:grid.n_cells])) < 2 * eps()
    @test maximum(abs.(grid_info["cell_volume"][:] - [grid.cells[i].V for i in 1:grid.n_cells])) < 2 * eps()

    close(grid_info)
    rm(grid_path)
end