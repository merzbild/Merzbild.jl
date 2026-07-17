@testset "malloc: grid N:2 merging on 1D grid with particles with dim(x)=1,2,3" begin

    # Important!
    # The time scaling in the analytical solution is different
    # Tref = 273.0
    # mref = 66.3e-27 
    # mcd = mref / 2.0
    # dref = 4.11e-10
    # nref = 1e23
    # Lref = 1.0 / (nref * constants.pi * dref**2)
    # vref = ((2 * constants.k * Tref) / mref)**0.5
    # time_ref = Lref / vref

    # kappa_mult = constants.pi * dref**2 * (mcd / (2 * constants.k * tref))**(-0.5) / gamma(5/2 - 1.0)
    # ttt_bkw = 1 / (4 * constants.pi * n * kappa_mult)
    # magic_factor = time_ref / ttt_bkw / (4 * constants.pi)
    # print(magic_factor)  # approximately 1.59577 for Argon, 1.5963 for N
        
    # def analytic(time, N):    
    #     C = 1. - 0.4 * np.exp(-time * magic_factor / 6)
    #     kk = N // 2
    #     return C**(kk - 1) * (kk - (kk - 1) * C)
        
    seed = 1234
    rng = StableRNG(seed)

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, "Ar")

    interaction_data_path = joinpath(@__DIR__, "..", "data", "pseudo_maxwell.toml")
    interaction_data::Array{Interaction, 2} = load_interaction_data(interaction_data_path, species_data)

    dt_scaled = 0.025

    np_base = 2000

    nx = 10
    grid = Grid1DUniform(0.5, 10)
    phys_props = PhysProps(grid.n_cells, 1)

    # init particle vector, particle indexer, grid particle sorter
    n_particles = np_base * nx

    # sample particles
    ndens = 5000.0
    T = 300.0

    Fnum = grid.cells[1].V * ndens / np_base

    @testset "1D particles, 1D merging" begin
        # 1d particle x vector
        mg1 = GridN2Merge{1}(6, 6, 6, 3.0)
        target = 1000

        particles1 = [ParticleVector{1}(n_particles)]
        pia = ParticleIndexerArray(grid.n_cells, 1)

        sample_particles_equal_weight!(rng, grid, particles1[1], pia, 1,
                                       species_data, ndens, T, Fnum)

        merges = false
        compute_props!(particles1, pia, species_data, phys_props)
        for cell in 1:nx
            if pia.indexer[cell,1].n_local > target
                merge_grid_based!(rng, mg1, particles1[1], pia, 1, 1, species_data, phys_props, grid)
                merges = true
            end
        end
        @test merges == true

        # reset particles
        particles1 = [ParticleVector{1}(n_particles)]
        pia = ParticleIndexerArray(grid.n_cells, 1)

        sample_particles_equal_weight!(rng, grid, particles1[1], pia, 1,
                                       species_data, ndens, T, Fnum)

        for cell in 1:nx
            if pia.indexer[cell,1].n_local > target
                bytes_merge = @allocated merge_grid_based!(rng, mg1, particles1[1], pia, 1, 1, species_data, phys_props, grid)
                merges = true
                @test bytes_merge <= 80
            end
        end

        @test merges == true

        # now try merging with explicit velocity bounds
        # reset particles
        particles1 = [ParticleVector{1}(n_particles)]
        pia = ParticleIndexerArray(grid.n_cells, 1)

        sample_particles_equal_weight!(rng, grid, particles1[1], pia, 1,
                                       species_data, ndens, T, Fnum)

        vlim = SVector{2,Float64}(-1000.0, 1000.0)
        for cell in 1:nx
            if pia.indexer[cell,1].n_local > target
                merge_grid_based!(rng, mg1, particles1[1], pia, 1, 1, species_data, vlim,
                                  vlim, vlim, grid)
                merges = true
            end
        end

        @test merges == true
        # reset particles
        particles1 = [ParticleVector{1}(n_particles)]
        pia = ParticleIndexerArray(grid.n_cells, 1)

        sample_particles_equal_weight!(rng, grid, particles1[1], pia, 1,
                                       species_data, ndens, T, Fnum)

        for cell in 1:nx
            if pia.indexer[cell,1].n_local > target
                bytes_merge = @allocated merge_grid_based!(rng, mg1, particles1[1], pia, 1, 1, species_data, vlim,
                                                           vlim, vlim, grid)
                @test bytes_merge <= 180
                merges = true
            end
        end

        @test merges == true
    end

    @testset "2D particles, 2D merging" begin
        # 2d particle x vector
        mg2 = GridN2Merge{2}(6, 6, 6, 3.0)
        target = 1000

        particles2 = [ParticleVector{2}(n_particles)]
        pia = ParticleIndexerArray(grid.n_cells, 1)

        sample_particles_equal_weight!(rng, grid, particles2[1], pia, 1,
                                       species_data, ndens, T, Fnum)

        merges = false
        compute_props!(particles2, pia, species_data, phys_props)
        for cell in 1:nx
            if pia.indexer[cell,1].n_local > target
                merge_grid_based!(rng, mg2, particles2[1], pia, 1, 1, species_data, phys_props, grid)
                merges = true
            end
        end
        @test merges == true

        # reset particles
        particles2 = [ParticleVector{2}(n_particles)]
        pia = ParticleIndexerArray(grid.n_cells, 1)

        sample_particles_equal_weight!(rng, grid, particles2[1], pia, 1,
                                       species_data, ndens, T, Fnum)

        for cell in 1:nx
            if pia.indexer[cell,1].n_local > target
                bytes_merge = @allocated merge_grid_based!(rng, mg2, particles2[1], pia, 1, 1, species_data, phys_props, grid)
                merges = true
                @test bytes_merge <= 80
            end
        end

        @test merges == true

        # now try merging with explicit velocity bounds
        # reset particles
        particles2 = [ParticleVector{2}(n_particles)]
        pia = ParticleIndexerArray(grid.n_cells, 1)

        sample_particles_equal_weight!(rng, grid, particles2[1], pia, 1,
                                       species_data, ndens, T, Fnum)

        vlim = SVector{2,Float64}(-1000.0, 1000.0)
        for cell in 1:nx
            if pia.indexer[cell,1].n_local > target
                merge_grid_based!(rng, mg2, particles2[1], pia, 1, 1, species_data, vlim,
                                  vlim, vlim, grid)
                merges = true
            end
        end

        @test merges == true
        # reset particles
        particles2 = [ParticleVector{2}(n_particles)]
        pia = ParticleIndexerArray(grid.n_cells, 1)

        sample_particles_equal_weight!(rng, grid, particles2[1], pia, 1,
                                       species_data, ndens, T, Fnum)

        for cell in 1:nx
            if pia.indexer[cell,1].n_local > target
                bytes_merge = @allocated merge_grid_based!(rng, mg2, particles2[1], pia, 1, 1, species_data, vlim,
                                                           vlim, vlim, grid)
                @test bytes_merge <= 180
                merges = true
            end
        end

        @test merges == true
    end

    @testset "3D particles, 3D merging" begin
        # 3d particle x vector
        mg3 = GridN2Merge{3}(6, 6, 6, 3.0)
        target = 1000

        particles3 = [ParticleVector{3}(n_particles)]
        pia = ParticleIndexerArray(grid.n_cells, 1)

        sample_particles_equal_weight!(rng, grid, particles3[1], pia, 1,
                                       species_data, ndens, T, Fnum)

        merges = false
        compute_props!(particles3, pia, species_data, phys_props)
        for cell in 1:nx
            if pia.indexer[cell,1].n_local > target
                merge_grid_based!(rng, mg3, particles3[1], pia, 1, 1, species_data, phys_props, grid)
                merges = true
            end
        end
        @test merges == true

        # reset particles
        particles3 = [ParticleVector{3}(n_particles)]
        pia = ParticleIndexerArray(grid.n_cells, 1)

        sample_particles_equal_weight!(rng, grid, particles3[1], pia, 1,
                                       species_data, ndens, T, Fnum)

        for cell in 1:nx
            if pia.indexer[cell,1].n_local > target
                bytes_merge = @allocated merge_grid_based!(rng, mg3, particles3[1], pia, 1, 1, species_data, phys_props, grid)
                merges = true
                @test bytes_merge <= 80
            end
        end

        @test merges == true

        # now try merging with explicit velocity bounds
        # reset particles
        particles3 = [ParticleVector{3}(n_particles)]
        pia = ParticleIndexerArray(grid.n_cells, 1)

        sample_particles_equal_weight!(rng, grid, particles3[1], pia, 1,
                                       species_data, ndens, T, Fnum)

        vlim = SVector{2,Float64}(-1000.0, 1000.0)
        for cell in 1:nx
            if pia.indexer[cell,1].n_local > target
                merge_grid_based!(rng, mg3, particles3[1], pia, 1, 1, species_data, vlim,
                                  vlim, vlim, grid)
                merges = true
            end
        end

        @test merges == true
        # reset particles
        particles3 = [ParticleVector{3}(n_particles)]
        pia = ParticleIndexerArray(grid.n_cells, 1)

        sample_particles_equal_weight!(rng, grid, particles3[1], pia, 1,
                                       species_data, ndens, T, Fnum)

        for cell in 1:nx
            if pia.indexer[cell,1].n_local > target
                bytes_merge = @allocated merge_grid_based!(rng, mg3, particles3[1], pia, 1, 1, species_data, vlim,
                                                           vlim, vlim, grid)
                @test bytes_merge <= 180
                merges = true
            end
        end

        @test merges == true
    end
end