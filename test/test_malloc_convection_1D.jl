
@testset "malloc tests: 3D particles, 1D convection and surface computes" begin

    T_wall = 300.0
    v_wall = 500.0
    L = 5e-4
    ndens = 5e22
    nx = 50
    ppc = 1000
    Δt = 2.59e-9
    output_freq = 1000
    n_timesteps = 6000

    seed = 1234
    rng = StableRNG(seed)

    # load particle and interaction data
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")
    interaction_data_path = joinpath(@__DIR__, "..", "data", "vhs.toml")
    interaction_data::Array{Interaction, 2} = load_interaction_data(interaction_data_path, species_data)

    # create our grid and BCs
    grid = Grid1DUniform(L, nx)
    bc_list = (MaxwellWallBC1D(species_data, 1, T_wall, [0.0, -v_wall, 0.0], 1.0),
               MaxwellWallBC1D(species_data, 1, T_wall, [0.0, v_wall, 0.0], 1.0))

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

    surf_props = SurfProps(pia, grid)
    phys_props = PhysProps(pia)

    for t in 1:2
        # convect particles
        convect_particles!(rng, grid, bc_list, particles[1], pia, 1, species_data, Δt)

        # sort particles
        sort_particles!(gridsorter, grid, particles[1], pia, 1)

        compute_props_sorted!(particles, pia, species_data, phys_props)
    end

    for t in 1:2
        # convect particles
        bytes = @allocated convect_particles!(rng, grid, bc_list, particles[1], pia, 1, species_data, Δt)
        @test bytes == 0

        # sort particles
        bytes = @allocated sort_particles!(gridsorter, grid, particles[1], pia, 1)
        @test bytes == 0

        bytes = @allocated compute_props_sorted!(particles, pia, species_data, phys_props)
        @test bytes == 0
    end

    for t in 1:2
        # convect particles
        convect_particles!(rng, grid, bc_list, particles[1], pia, 1, species_data, surf_props, Δt)

        # sort particles
        sort_particles!(gridsorter, grid, particles[1], pia, 1)
    end

    for t in 1:2
        # convect particles
        bytes = @allocated convect_particles!(rng, grid, bc_list, particles[1], pia, 1, species_data, surf_props, Δt)
        @test bytes == 0

        # sort particles
        bytes = @allocated sort_particles!(gridsorter, grid, particles[1], pia, 1)
        @test bytes == 0
    end

    # now with precomputation of particle indices
    for t in 1:2
        # convect particles
        convect_particles_and_compute_cell!(rng, grid, bc_list, particles[1], pia, 1, species_data, Δt)

        # sort particles
        sort_particles!(gridsorter, particles[1], pia, 1)
    end

    for t in 1:2
        # convect particles
        bytes = @allocated convect_particles_and_compute_cell!(rng, grid, bc_list, particles[1], pia, 1, species_data, Δt)
        @test bytes == 0

        # sort particles
        bytes = @allocated sort_particles!(gridsorter, particles[1], pia, 1)
        @test bytes == 0
    end

    # now with precomputation of particle indices and computation of surface properties
    for t in 1:2
        # convect particles
        convect_particles_and_compute_cell!(rng, grid, bc_list, particles[1], pia, 1, species_data, surf_props, Δt)

        # sort particles
        sort_particles!(gridsorter, particles[1], pia, 1)
    end

    for t in 1:2
        # convect particles
        bytes = @allocated convect_particles_and_compute_cell!(rng, grid, bc_list, particles[1], pia, 1, species_data, surf_props, Δt)
        @test bytes == 0

        # sort particles
        bytes = @allocated sort_particles!(gridsorter, particles[1], pia, 1)
        @test bytes == 0
    end

    flux_props = FluxProps(pia)
    compute_flux_props!(particles, pia, species_data, phys_props, flux_props, grid)

    bytes = @allocated compute_flux_props!(particles, pia, species_data, phys_props, flux_props, grid)
    @test bytes == 0
end
