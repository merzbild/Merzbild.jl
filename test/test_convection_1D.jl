@testset "convection 1D" begin
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")

    seed = 1234
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    grid = Grid1DUniform(50.0, 100)

    gridsorter = GridSortInPlace(grid, 15000)
    particles = [ParticleVector(4)]

    pia = ParticleIndexerArray(grid.n_cells, 1)
    pia.n_total[1] = 4

    # will just move: new x_coord = 20.5
    particles[1][1] = Particle(1.0, [-1.25, -1.5, 4.0], [23.0, -8.0, 7.5])

    # will reflect from right wall: new x_cord = 50.0 - (10.0 + 11.0) = 29.0
    particles[1][2] = Particle(2.0, [11.0, -3.0, 1.0], [49.0, 6.0, -3.0])

    # will reflect from left wall: new x_cord = 3.0 + 20.0 = 23.0
    particles[1][3] = Particle(3.0, [-20.0, 0.0, 2.0], [17.0, 1.0, 3.0])

    # will reflect from left wall and from right wall: new x_coord = 3.5
    particles[1][4] = Particle(4.0, [-49.0, -20.0, 13.0], [1.5, -1.0, 9.0])

    # 1D specularly reflecting boundaries
    boundaries = MaxwellWalls1D(species_data, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0)

    convect_particles!(rng, grid, boundaries, particles[1], pia, 1, species_data, 2.0)
    sort_particles!(gridsorter, grid, particles[1], pia, 1)

    @test particles[1].index == [4, 1, 3, 2]

    @test maximum(abs.(particles[1][1].x - [3.5, -1.0, 9.0])) < 2 * eps()
    @test particles[1][1].v == [-49.0, -20.0, 13.0]
    @test particles[1][1].w == 4.0

    @test maximum(abs.(particles[1][2].x - [20.5, -8.0, 7.5])) < 2 * eps()
    @test particles[1][2].v == [-1.25, -1.5, 4.0]
    @test particles[1][2].w == 1.0

    @test maximum(abs.(particles[1][3].x - [23.0, 1.0, 3.0])) < 2 * eps()
    @test particles[1][3].v == [20.0, 0.0, 2.0]
    @test particles[1][3].w == 3.0

    @test maximum(abs.(particles[1][4].x - [29.0, 6.0, -3.0])) < 2 * eps()
    @test particles[1][4].v == [-11.0, -3.0, 1.0]
    @test particles[1][4].w == 2.0
    phys_props = PhysProps(grid.n_cells, 1, [], Tref=1)
    compute_props!(particles, pia, species_data, phys_props)

    for i in 1:grid.n_cells
        if i == 8
            w = 4.0
        elseif i == 42
            w = 1.0
        elseif i == 47
            w = 3.0
        elseif i == 59
            w = 2.0
        else
            w = 0.0
        end
        @test abs(phys_props.n[i, 1] - w) < eps()
    end

    # sample a lot of particles, move them so that they hit the wall
    # and sample from a half-Maxwellian and test the resulting distribution
    n_particles = 10000
    particles = [ParticleVector(n_particles)]
    pia = ParticleIndexerArray(grid.n_cells, 1)
    pia.n_total[1] = n_particles
    n = 1e10

    for i in 1:n_particles
        particles[1][i] = Particle(n / n_particles, [-1000.0, 0.0, 0.0], [0.999e-4, 0.0, 0.0])
    end

    Δt = 1e-7
    vy_left_wall = 1100.0
    vy_right_wall = -820.0
    boundaries_acc = MaxwellWalls1D(species_data, 2000.0, 500.0, vy_left_wall, vy_right_wall, 1.0, 1.0)

    convect_particles!(rng, grid, boundaries_acc, particles[1], pia, 1, species_data, Δt)
    sort_particles!(gridsorter, grid, particles[1], pia, 1)

    sign_pos = true
    for i in 1:n_particles
        if particles[1][i].v[1] < 0
            sign_pos = false
        end
    end

    # no reflected particles have a negative velocity
    @test sign_pos == true

    compute_props!(particles, pia, species_data, phys_props)
    @test abs(phys_props.n[1, 1] - n) < eps()

    # we sampled a lot of particles, reflected y velocity should be very close to prescribed
    @test abs((phys_props.v[2, 1, 1] - vy_left_wall) / vy_left_wall) < 2.25e-3

    # z-velocity should be almost zero
    @test abs(phys_props.v[3, 1, 1] - 0.0) < 3.0
    
    # now we test reflection from right wall
    particles[1].index = Vector(1:n_particles)

    for i in 1:n_particles
        particles[1][i] = Particle(n / n_particles, [1000.0, 0.0, 0.0], [50. - 0.999e-4, 0.0, 0.0])
    end

    convect_particles!(rng, grid, boundaries_acc, particles[1], pia, 1, species_data, Δt)
    sort_particles!(gridsorter, grid, particles[1], pia, 1)

    sign_neg = true
    for i in 1:n_particles
        if particles[1][i].v[1] > 0
            sign_neg = false
        end
    end

    # no reflected particles have a positive velocity
    @test sign_neg == true

    compute_props!(particles, pia, species_data, phys_props)
    @test abs(phys_props.n[100, 1] - n) < eps()

    # we sampled a lot of particles, reflected y velocity should be very close to prescribed
    @test abs((phys_props.v[2, 100, 1] - vy_right_wall) / vy_right_wall) < 2.25e-3

    # z-velocity should be almost zero
    @test abs(phys_props.v[3, 100, 1] - 0.0) < 3.0


    # now we try out a very cold accomodating wall
    # cold so that the chances of a particle reflected from the wall with a velocity of 1000.0
    # are virtually zero
    for i in 1:n_particles
        particles[1][i] = Particle(n / n_particles, [-1000.0, 0.0, 0.0], [0.999e-4, 0.0, 0.0])
    end

    Δt = 1e-7
    vy_left_wall = 1100.0
    vy_right_wall = -820.0
    # accomodation coefficient of 0.2 - 80% chance of specular reflection
    boundaries_acc = MaxwellWalls1D(species_data, 10.0, 10.0, vy_left_wall, vy_right_wall, 0.2, 1.0)

    convect_particles!(rng, grid, boundaries_acc, particles[1], pia, 1, species_data, Δt)
    sort_particles!(gridsorter, grid, particles[1], pia, 1)

    sign_pos = true
    nspecular = 0
    for i in 1:n_particles
        if particles[1][i].v[1] < 0
            sign_pos = false
        end
        if abs(particles[1][i].v[1] - 1000.0) < 2 * eps()
            nspecular += 1
        end
    end

    # no reflected particles have a negative velocity
    @test sign_pos == true

    # approximately 80% of particles are reflected with specular reflection
    @test abs(nspecular - 8000) < 100
end