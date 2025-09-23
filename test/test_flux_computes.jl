@testset "flux_props compute" begin

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, "Ar")

    grid = Grid1DUniform(4.0, 8)

    particles = [ParticleVector(6)]

    # create 4 particles in cell 1
    velocities_x = [1.0, -1.0, -2.0, 6.0]  # mean_vx = 1.0
    velocities_y = [-3.0, -4.0, 3.0, 4.0]
    velocities_z = [0.0, 0.0, 0.0, 0.0]

    # cx^2 = 38.0 
    for i in 1:4
        Merzbild.add_particle!(particles[1], i,
                               1.0, [velocities_x[i], velocities_y[i], velocities_z[i]],
                               [0.3, 0.0, 0.0])
    end

    # create 2 particles in cell 6
    velocities_x = [8.0, -2.0]
    velocities_y = [0.0, 0.0]
    velocities_z = [2.0, -2.0]
    weights = [2.0, 8.0]
    for i in 5:6
        Merzbild.add_particle!(particles[1], i,
                               weights[i-4], [velocities_x[i-4], velocities_y[i-4], velocities_z[i-4]],
                               [2.8, 0.0, 0.0])
    end

    pia = ParticleIndexerArray(8,1)

    # 2 particles in group1 in cell1, 2 particles in group2 in cell1
    pia.n_total[1] = 8
    pia.indexer[1,1].n_local = 4
    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = 2
    pia.indexer[1,1].n_group1 = 2
    pia.indexer[1,1].start2 = 3
    pia.indexer[1,1].end2 = 4
    pia.indexer[1,1].n_group2 = 2

    pia.indexer[6,1].n_local = 2
    pia.indexer[6,1].start1 = 5
    pia.indexer[6,1].end1 = 6
    pia.indexer[6,1].n_group1 = 2

    phys_props = PhysProps(pia)
    flux_props = FluxProps(pia)
    @test flux_props.n_cells == 8
    @test flux_props.n_species == 1
    @test size(flux_props.kinetic_energy_flux) == (3,8,1)
    @test size(flux_props.diagonal_momentum_flux) == (3,8,1)
    @test size(flux_props.off_diagonal_momentum_flux) == (3,8,1)

    compute_props!(particles, pia, species_data, phys_props)
    compute_flux_props!(particles, pia, species_data, phys_props, flux_props, grid)

    # we don't write stuff to cells with no particles
    for cell in [2,3,4,5,7,8]
        @test flux_props.kinetic_energy_flux[:,cell,1] == [0.0, 0.0, 0.0]
        @test flux_props.diagonal_momentum_flux[:,cell,1] == [0.0, 0.0, 0.0]
        @test flux_props.off_diagonal_momentum_flux[:,cell,1] == [0.0, 0.0, 0.0]
    end

    kef = 0.5 * [38.0, 0.0, 0.0] * species_data[1].mass / 0.5
    @test maximum(abs.(flux_props.kinetic_energy_flux[:,1,1] - kef)) <= eps()
    dmf = [9.5, 12.5, 0.0] * species_data[1].mass / 0.5
    @test maximum(abs.(flux_props.diagonal_momentum_flux[:,1,1] - dmf)) <= eps()
    odmf = [4.75, 0.0, 0.0] * species_data[1].mass / 0.5
    @test maximum(abs.(flux_props.off_diagonal_momentum_flux[:,1,1] - odmf)) <= eps()

    kef2 = 0.5 * [960.0, 0.0, 7.68] * species_data[1].mass / 0.5
    @test maximum(abs.(flux_props.kinetic_energy_flux[:,6,1] - kef2)) <= eps()
    dmf2 = [160.0, 0.0, 6.4] * species_data[1].mass / 0.5
    @test maximum(abs.(flux_props.diagonal_momentum_flux[:,1,1] - dmf2)) <= eps()
    odmf2 = [0.0, 32.0, 0.0] * species_data[1].mass / 0.5
    @test maximum(abs.(flux_props.off_diagonal_momentum_flux[:,1,1] - odmf2)) <= eps()

    clear_props!(flux_props)
    for cell in 1:8
        @test flux_props.kinetic_energy_flux[:,cell,1] == [0.0, 0.0, 0.0]
        @test flux_props.diagonal_momentum_flux[:,cell,1] == [0.0, 0.0, 0.0]
        @test flux_props.off_diagonal_momentum_flux[:,cell,1] == [0.0, 0.0, 0.0]
    end

    flux_props_avg = FluxProps(pia)
    flux_props.kinetic_energy_flux[:,:,:] .= 1.0
    flux_props.diagonal_momentum_flux[:,:,:] .= 2.0
    flux_props.off_diagonal_momentum_flux[:,:,:] .= 3.0

    avg_props!(flux_props_avg, flux_props, 2)
    flux_props.kinetic_energy_flux[:,:,:] .= 3.0
    flux_props.diagonal_momentum_flux[:,:,:] .= 5.0
    flux_props.off_diagonal_momentum_flux[:,:,:] .= -3.0

    avg_props!(flux_props_avg, flux_props, 2)
    for cell in 1:8
        @test flux_props_avg.kinetic_energy_flux[:,cell,1] == [2.0, 2.0, 2.0]
        @test flux_props_avg.diagonal_momentum_flux[:,cell,1] == [3.5, 3.5, 3.5]
        @test flux_props_avg.off_diagonal_momentum_flux[:,cell,1] == [0.0, 0.0, 0.0]
    end
    clear_props!(flux_props)

    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = 4
    pia.indexer[1,1].n_group1 = 4
    pia.indexer[1,1].start2 = 0
    pia.indexer[1,1].end2 = -1
    pia.indexer[1,1].n_group2 = 0

    # sorted particles, whole grid
    compute_flux_props_sorted!(particles, pia, species_data, phys_props, flux_props, grid)

    # we don't write stuff to cells with no particles
    for cell in [2,3,4,5,7,8]
        @test flux_props.kinetic_energy_flux[:,cell,1] == [0.0, 0.0, 0.0]
        @test flux_props.diagonal_momentum_flux[:,cell,1] == [0.0, 0.0, 0.0]
        @test flux_props.off_diagonal_momentum_flux[:,cell,1] == [0.0, 0.0, 0.0]
    end

    @test maximum(abs.(flux_props.kinetic_energy_flux[:,1,1] - kef)) <= eps()
    @test maximum(abs.(flux_props.diagonal_momentum_flux[:,1,1] - dmf)) <= eps()
    @test maximum(abs.(flux_props.off_diagonal_momentum_flux[:,1,1] - odmf)) <= eps()

    @test maximum(abs.(flux_props.kinetic_energy_flux[:,6,1] - kef2)) <= eps()
    @test maximum(abs.(flux_props.diagonal_momentum_flux[:,1,1] - dmf2)) <= eps()
    @test maximum(abs.(flux_props.off_diagonal_momentum_flux[:,1,1] - odmf2)) <= eps()

    # finally, compute in chunk 3:7, so cell is skipped
    clear_props!(flux_props)

    compute_flux_props_sorted!(particles, pia, species_data, phys_props, flux_props, grid, 3:7)
    for cell in [1,2,3,4,5,7,8]
        @test flux_props.kinetic_energy_flux[:,cell,1] == [0.0, 0.0, 0.0]
        @test flux_props.diagonal_momentum_flux[:,cell,1] == [0.0, 0.0, 0.0]
        @test flux_props.off_diagonal_momentum_flux[:,cell,1] == [0.0, 0.0, 0.0]
    end

    @test maximum(abs.(flux_props.kinetic_energy_flux[:,6,1] - kef2)) <= eps()
    @test maximum(abs.(flux_props.diagonal_momentum_flux[:,1,1] - dmf2)) <= eps()
    @test maximum(abs.(flux_props.off_diagonal_momentum_flux[:,1,1] - odmf2)) <= eps()
end