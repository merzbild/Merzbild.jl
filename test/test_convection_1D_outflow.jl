@testset "convection 1D with outflow boundaries" begin
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")
    mass = species_data[1].mass

    seed = 1234
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    grid = Grid1DUniform(50.0, 100)
    gridsorter = GridSortInPlace(grid, 15000)

    # the 4 particles of test_convection_1D.jl, but the right wall is now an outflow boundary:
    # particle 1 just moves to x = 20.5
    # particle 2 reaches the right wall and is deleted
    # particle 3 reflects from the left wall to x = 23.0
    # particle 4 reflects from the left wall, reaches the right wall and is deleted
    function setup_4_particles()
        pv = ParticleVector(4)
        pv.nbuffer = 0  # all 4 slots are in use

        pia_local = ParticleIndexerArray(grid.n_cells, 1)
        pia_local.n_total[1] = 4
        pia_local.index_last[1] = 4
        pia_local.indexer[1,1].n_local = 4
        pia_local.indexer[1,1].n_group1 = 4
        pia_local.indexer[1,1].start1 = 1
        pia_local.indexer[1,1].end1 = 4

        pv[1] = Particle(1.0, [-1.25, -1.5, 4.0], [23.0, -8.0, 7.5])
        pv[2] = Particle(2.0, [11.0, -3.0, 1.0], [49.0, 6.0, -3.0])
        pv[3] = Particle(3.0, [-20.0, 0.0, 2.0], [17.0, 1.0, 3.0])
        pv[4] = Particle(4.0, [-49.0, -20.0, 13.0], [1.55, -1.0, 9.0])

        return [pv], pia_local
    end

    bc_list = (FullySpecularBC1D(), OutflowBC())

    particles, pia = setup_4_particles()

    n_deleted = convect_particles!(rng, grid, bc_list, particles[1], pia, 1, species_data, 2.0)

    @test n_deleted == 2
    @test pia.n_total[1] == 2
    @test pia.index_last[1] == 2
    @test particles[1].nbuffer == 2
    # deletion from a cell that is not the last one breaks contiguity
    @test pia.contiguous[1] == false

    sort_particles!(gridsorter, grid, particles[1], pia, 1)

    @test check_pia_is_correct(pia, 1) == (true, 0)
    @test check_unique_index(particles[1], pia, 1) == (true, 0)
    @test check_unique_buffer(particles[1]) == (true, 0)

    # the remaining particles are the ones that never reached the outflow boundary
    @test maximum(abs.(particles[1][1].x - [20.5, -8.0, 7.5])) < 2 * eps()
    @test particles[1][1].v == [-1.25, -1.5, 4.0]
    @test particles[1][1].w == 1.0

    @test maximum(abs.(particles[1][2].x - [23.0, 1.0, 3.0])) < 2 * eps()
    @test particles[1][2].v == [20.0, 0.0, 2.0]
    @test particles[1][2].w == 3.0

    phys_props = PhysProps(grid.n_cells, 1)
    compute_props!(particles, pia, species_data, phys_props)

    # variable weights: the deleted particles took their weights (2.0 and 4.0) with them
    @test abs(sum(phys_props.n[:, 1]) - 4.0) < eps()
    for i in 1:grid.n_cells
        if i == 42
            w = 1.0
        elseif i == 47
            w = 3.0
        else
            w = 0.0
        end
        @test abs(phys_props.n[i, 1] - w) < eps()
    end

    # a reflecting-only bc_list must be unaffected by the new code path
    particles, pia = setup_4_particles()
    n_deleted = convect_particles!(rng, grid, (FullySpecularBC1D(), FullySpecularBC1D()),
                                   particles[1], pia, 1, species_data, 2.0)
    @test n_deleted == 0
    @test pia.n_total[1] == 4
    @test pia.contiguous[1] == true
    sort_particles!(gridsorter, grid, particles[1], pia, 1)
    @test particles[1].index == [4, 1, 3, 2]

    # all particles leave the domain
    particles, pia = setup_4_particles()
    n_deleted = convect_particles!(rng, grid, (OutflowBC(), OutflowBC()),
                                   particles[1], pia, 1, species_data, 2.0)

    @test n_deleted == 3  # particle 1 does not reach a wall
    @test pia.n_total[1] == 1
    @test particles[1].nbuffer == 3

    n_deleted = convect_particles!(rng, grid, (OutflowBC(), OutflowBC()),
                                   particles[1], pia, 1, species_data, 100.0)
    @test n_deleted == 1
    @test pia.n_total[1] == 0
    @test pia.index_last[1] == 0
    @test particles[1].nbuffer == 4

    # sorting an emptied domain must work
    sort_particles!(gridsorter, grid, particles[1], pia, 1)
    @test check_pia_is_correct(pia, 1) == (true, 0)
    @test check_unique_buffer(particles[1]) == (true, 0)
    compute_props!(particles, pia, species_data, phys_props)
    @test abs(sum(phys_props.n[:, 1])) < eps()

    # surface properties: an outflowing particle is counted as incident, but not as reflected
    particles, pia = setup_4_particles()
    surf_props = SurfProps(pia, grid)
    Δt = 2.0

    n_deleted = convect_particles!(rng, grid, bc_list, particles[1], pia, 1, species_data, surf_props, Δt)
    @test n_deleted == 2

    factor = mass / Δt

    # right wall: particles 2 and 4 leave through it
    @test surf_props.np[2, 1] == 2
    @test surf_props.flux_reflected[2, 1] == 0.0
    @test abs(surf_props.flux_incident[2, 1] - (2.0 + 4.0) * factor) < eps(6.0 * factor)
    # momentum of particle 2 (w=2, vx=11) and particle 4 (w=4, vx=49 after the specular reflection)
    @test abs(surf_props.force[1, 2, 1] - (2.0 * 11.0 + 4.0 * 49.0) * factor) < 1e-8 * factor
    @test abs(surf_props.kinetic_energy_flux[2, 1] -
              (0.5 * 2.0 * (121.0 + 9.0 + 1.0) + 0.5 * 4.0 * (2401.0 + 400.0 + 169.0)) * factor) < 1e-8 * factor

    # left wall: particles 3 and 4 reflect from it, so it is unaffected by the outflow boundary
    @test surf_props.np[1, 1] == 2
    @test abs(surf_props.flux_incident[1, 1] - (3.0 + 4.0) * factor) < eps(7.0 * factor)
    @test abs(surf_props.flux_reflected[1, 1] + (3.0 + 4.0) * factor) < eps(7.0 * factor)

    # convection with in-routine cell computation: the 3-argument (cell-based) sort must
    # agree with the 4-argument (position-based) one, which requires particles.cell to be
    # kept in sync with particles.index when a deletion swaps particles around
    particles, pia = setup_4_particles()
    n_deleted = convect_particles_and_compute_cell!(rng, grid, bc_list, particles[1], pia, 1, species_data, 2.0)
    @test n_deleted == 2
    @test particles[1].cell[1] == 42
    @test particles[1].cell[2] == 47

    sort_particles!(gridsorter, particles[1], pia, 1)

    @test check_pia_is_correct(pia, 1) == (true, 0)
    @test check_unique_index(particles[1], pia, 1) == (true, 0)
    @test check_unique_buffer(particles[1]) == (true, 0)

    phys_props_cell = PhysProps(grid.n_cells, 1)
    compute_props!(particles, pia, species_data, phys_props_cell)

    particles_ref, pia_ref = setup_4_particles()
    convect_particles!(rng, grid, bc_list, particles_ref[1], pia_ref, 1, species_data, 2.0)
    sort_particles!(gridsorter, grid, particles_ref[1], pia_ref, 1)
    phys_props_ref = PhysProps(grid.n_cells, 1)
    compute_props!(particles_ref, pia_ref, species_data, phys_props_ref)

    @test phys_props_cell.n[:, 1] == phys_props_ref.n[:, 1]
    @test particles[1].index[1:2] == particles_ref[1].index[1:2]

    # non-contiguous pia: particles are deleted from both index groups, and the particles
    # the pia does not point to must be left alone
    function setup_noncontiguous()
        pv = ParticleVector(30)
        pia_local = ParticleIndexerArray(grid.n_cells, 1)

        for i in 1:30
            if i in 11:25
                # particles the pia does not point to
                Merzbild.add_particle!(pv, i, 10000.0, SVector{3,Float64}(-1000.0, 0.0, 0.0),
                                       SVector{3,Float64}(0.75, 0.0, 0.0))
            else
                # even-indexed tracked particles are fast enough to reach the outflow boundary
                vx = i % 2 == 0 ? 100.0 : 1.0
                Merzbild.add_particle!(pv, i, 1.0, SVector{3,Float64}(vx, 0.0, 0.0),
                                       SVector{3,Float64}(0.75, 0.0, 0.0))
            end
        end

        pia_local.n_total[1] = 15
        pia_local.index_last[1] = 30
        pia_local.indexer[2,1].n_local = 15
        pia_local.indexer[2,1].n_group1 = 10
        pia_local.indexer[2,1].start1 = 1
        pia_local.indexer[2,1].end1 = 10
        pia_local.indexer[2,1].n_group2 = 5
        pia_local.indexer[2,1].start2 = 26
        pia_local.indexer[2,1].end2 = 30
        pia_local.contiguous[1] = false

        return [pv], pia_local
    end

    # deleted: 5 from group1 (indices 2, 4, 6, 8, 10) and 3 from group2 (indices 26, 28, 30)
    function test_noncontiguous_deletion(particles, pia, n_deleted)
        @test n_deleted == 8
        @test pia.n_total[1] == 7
        @test pia.indexer[2,1].n_group1 == 5
        @test pia.indexer[2,1].n_group2 == 2
        @test particles[1].nbuffer == 8
        @test check_pia_is_correct(pia, 1) == (true, 0)
        @test check_unique_index(particles[1], pia, 1) == (true, 0)
        @test check_unique_buffer(particles[1]) == (true, 0)

        # the untracked particles have not been convected
        for i in 11:25
            @test particles[1].particles[i].v[1] == -1000.0
            @test particles[1].particles[i].x[1] == 0.75
        end
    end

    # the 7 remaining tracked particles have all moved from x=0.75 to x=1.75, i.e. cell 4
    function test_noncontiguous_sorted(particles, pia)
        @test check_pia_is_correct(pia, 1) == (true, 0)
        @test check_unique_index(particles[1], pia, 1) == (true, 0)

        compute_props!(particles, pia, species_data, phys_props)
        @test abs(phys_props.n[4, 1] - 7.0) < eps()
        @test abs(sum(phys_props.n[:, 1]) - 7.0) < eps()
    end

    # the cell indices of the remaining particles, including those swapped into the slot of a
    # deleted particle, have to be the post-convection ones
    function test_noncontiguous_cells(particles, pia)
        indexer = pia.indexer[2,1]
        for i in indexer.start1:indexer.end1
            @test particles[1].cell[i] == 4
        end
        for i in indexer.start2:indexer.end2
            @test particles[1].cell[i] == 4
        end
    end

    particles, pia = setup_noncontiguous()
    @test check_pia_is_correct(pia, 1) == (true, 0)

    n_deleted = convect_particles!(rng, grid, bc_list, particles[1], pia, 1, species_data, 1.0)
    test_noncontiguous_deletion(particles, pia, n_deleted)

    sort_particles!(gridsorter, grid, particles[1], pia, 1)
    test_noncontiguous_sorted(particles, pia)

    particles, pia = setup_noncontiguous()
    n_deleted = convect_particles_and_compute_cell!(rng, grid, bc_list, particles[1], pia, 1, species_data, 1.0)
    test_noncontiguous_deletion(particles, pia, n_deleted)
    test_noncontiguous_cells(particles, pia)

    sort_particles!(gridsorter, particles[1], pia, 1)
    test_noncontiguous_sorted(particles, pia)

    particles, pia = setup_noncontiguous()
    surf_props = SurfProps(pia, grid)
    n_deleted = convect_particles_and_compute_cell!(rng, grid, bc_list, particles[1], pia, 1, species_data,
                                                    surf_props, 1.0)
    test_noncontiguous_deletion(particles, pia, n_deleted)
    test_noncontiguous_cells(particles, pia)

    # all 8 deleted particles (w=1, vx=100) leave through the right wall, none reach the left one
    @test surf_props.np[1, 1] == 0
    @test surf_props.np[2, 1] == 8
    @test surf_props.flux_reflected[2, 1] == 0.0
    @test abs(surf_props.flux_incident[2, 1] - 8.0 * mass) < eps(8.0 * mass)
    @test abs(surf_props.force[1, 2, 1] - 800.0 * mass) < 1e-8 * mass

    sort_particles!(gridsorter, particles[1], pia, 1)
    test_noncontiguous_sorted(particles, pia)

    # free-molecular effusion: the mass leaving a domain with one outflow boundary over a
    # single timestep must match the analytic number flux n * v_mean / 4
    T0 = 300.0
    L = 5e-4
    ndens = 5e22
    nx = 50
    ppc = 2000

    grid_eff = Grid1DUniform(L, nx)
    rng_eff = StableRNG(seed)

    n_particles = ppc * nx
    particles_eff = [ParticleVector(n_particles)]
    pia_eff = ParticleIndexerArray(grid_eff.n_cells, 1)
    gridsorter_eff = GridSortInPlace(grid_eff, n_particles)

    Fnum = grid_eff.cells[1].V * ndens / ppc
    sample_particles_equal_weight!(rng_eff, grid_eff, particles_eff[1], pia_eff, 1,
                                   species_data, ndens, T0, Fnum)

    v_mean = sqrt(8 * Merzbild.k_B * T0 / (pi * mass))
    Δt_eff = 0.04 * L / v_mean  # a particle travels far less than a domain length per step

    n_deleted = convect_particles!(rng_eff, grid_eff, (FullySpecularBC1D(), OutflowBC()),
                                   particles_eff[1], pia_eff, 1, species_data, Δt_eff)

    n_expected = n_particles * v_mean * Δt_eff / (4 * L)
    @test abs(n_deleted - n_expected) / n_expected < 0.1
    @test pia_eff.n_total[1] == n_particles - n_deleted
    @test particles_eff[1].nbuffer == n_deleted

    sort_particles!(gridsorter_eff, grid_eff, particles_eff[1], pia_eff, 1)
    @test check_pia_is_correct(pia_eff, 1) == (true, 0)
    @test check_unique_index(particles_eff[1], pia_eff, 1) == (true, 0)
    @test check_unique_buffer(particles_eff[1]) == (true, 0)
end
