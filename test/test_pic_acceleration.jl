@testset "PIC acceleration" begin
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, ["e-", "He+", "Ar"])

    L = 1.0
    nx = 4
    grid = Grid1DUniform(L, nx)
    Δx = grid.Δx
    Δt = 1e-9

    n_particles = 6
    particles = [ParticleVector(n_particles)]
    pia = ParticleIndexerArray(grid.n_cells, 1)
    gridsorter = GridSortInPlace(grid, n_particles)

    x_p = [0.13 * Δx, 1.62 * Δx, 2.0 * Δx, 2.85 * Δx, 3.5 * Δx, 3.99 * Δx]

    for i in 1:n_particles
        particles[1][i] = Particle(1.0, [0.0, 3.0, -4.0], [x_p[i], 0.0, 0.0])
    end
    pia.n_total[1] = n_particles
    pia.index_last[1] = n_particles
    pia.indexer[1,1].n_local = n_particles
    pia.indexer[1,1].n_group1 = n_particles
    pia.indexer[1,1].start1 = 1
    pia.indexer[1,1].end1 = n_particles
    pia.contiguous[1] = true

    sort_particles!(gridsorter, grid, particles[1], pia, 1)

    # linear interpolation is exact for a linear field
    field_props = ElectrostaticFieldProps(grid)
    E0 = 1500.0
    E_slope = -2.5e4
    for j in 1:field_props.n_nodes
        field_props.electric_field[j] = E0 + E_slope * (j - 1) * Δx
    end

    accelerate_electric_field_x!(grid, particles[1], pia, 1, species_data, field_props, Δt)

    charge_div_mass = species_data[1].charge_div_mass
    for i in 1:n_particles
        xp = particles[1][i].x[1]
        vx_exact = charge_div_mass * (E0 + E_slope * xp) * Δt
        @test abs(particles[1][i].v[1] - vx_exact) < 1e-13 * abs(charge_div_mass * E0 * Δt)
        @test particles[1][i].v[2] == 3.0
        @test particles[1][i].v[3] == -4.0
    end

    # a uniform field reproduces the constant-field acceleration exactly
    particles_const = [ParticleVector(n_particles)]
    pia_const = ParticleIndexerArray(grid.n_cells, 1)

    grid_single_cell = Grid1DUniform(L, 1)

    for i in 1:n_particles
        particles[1][i] = Particle(1.0, [0.0, 3.0, -4.0], [x_p[i], 0.0, 0.0])
        particles_const[1][i] = Particle(1.0, [0.0, 3.0, -4.0], [x_p[i], 0.0, 0.0])
    end

    for p in (pia, pia_const)
        for cell in 1:grid.n_cells
            p.indexer[cell,1] = ParticleIndexer()
        end
        p.n_total[1] = n_particles
        p.index_last[1] = n_particles
        p.indexer[1,1].n_local = n_particles
        p.indexer[1,1].n_group1 = n_particles
        p.indexer[1,1].start1 = 1
        p.indexer[1,1].end1 = n_particles
        p.contiguous[1] = true
    end
    sort_particles!(gridsorter, grid, particles[1], pia, 1)
    sort_particles!(gridsorter, grid, particles_const[1], pia_const, 1)

    E_const = 997.0
    fill!(field_props.electric_field, E_const)

    field_props_single_cell = ElectrostaticFieldProps(grid_single_cell)
    fill!(field_props_single_cell.electric_field, E_const)

    accelerate_electric_field_x!(grid, particles[1], pia, 1, species_data, field_props, Δt)

    for cell in 1:grid.n_cells
        accelerate_constant_field_x!(particles_const[1], pia_const, cell, 1, species_data, E_const, Δt)
    end

    for i in 1:n_particles
        @test particles[1][i].v[1] == particles_const[1][i].v[1]
    end

    # self-force: a single particle in an otherwise empty periodic box does not accelerate itself
    L_sf = 0.01
    nx_sf = 32
    grid_sf = Grid1DUniform(L_sf, nx_sf)

    particles_sf = [ParticleVector(1)]
    pia_sf = ParticleIndexerArray(grid_sf.n_cells, 1)
    gridsorter_sf = GridSortInPlace(grid_sf, 1)

    field_props_sf = ElectrostaticFieldProps(grid_sf)
    poisson_solver_sf = PoissonSolver1DUniform(grid_sf, PeriodicFieldBC1D(), PeriodicFieldBC1D())

    w_sf = 1e8
    Δt_sf = 1e-12

    for x_offset in [0.0, 0.27, 0.5, 0.83]
        particles_sf[1][1] = Particle(w_sf, [0.0, 0.0, 0.0],
                                      [(7 + x_offset) * grid_sf.Δx, 0.0, 0.0])
        pia_sf.n_total[1] = 1
        pia_sf.index_last[1] = 1
        for cell in 1:grid_sf.n_cells
            pia_sf.indexer[cell,1] = ParticleIndexer()
        end
        pia_sf.indexer[1,1].n_local = 1
        pia_sf.indexer[1,1].n_group1 = 1
        pia_sf.indexer[1,1].start1 = 1
        pia_sf.indexer[1,1].end1 = 1
        pia_sf.contiguous[1] = true

        sort_particles!(gridsorter_sf, grid_sf, particles_sf[1], pia_sf, 1)

        deposit_charge!(poisson_solver_sf, grid_sf, particles_sf, pia_sf, species_data, field_props_sf)
        solve_poisson!(poisson_solver_sf, field_props_sf)

        accelerate_electric_field_x!(grid_sf, particles_sf[1], pia_sf, 1, species_data, field_props_sf, Δt_sf)

        # scale of the velocity change a particle would get from the field it deposits itself
        dv_scale = abs(species_data[1].charge_div_mass * maximum(abs.(field_props_sf.electric_field)) * Δt_sf)
        @test dv_scale > 0.0
        @test abs(particles_sf[1][1].v[1]) < 1e-12 * dv_scale
    end
end
