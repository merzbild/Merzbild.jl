@testset "merging_grid_indexing" begin

    species_list::Vector{Species} = load_species_list("data/particles.toml", "Ar")

    seed = 1234
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    Nx = 4
    Ny = 3
    Nz = 2
    mg = create_merging_grid(Nx, Ny, Nz, 1.0)

    phys_props::PhysProps = create_props(1, 1, [], Tref=1)

    phys_props.T[1,1] = species_list[1].mass / (2 * k_B)

    @test mg.Ntotal == 32  # Nx * Ny * Nz + 8
    @test mg.NyNz == 6  # Ny * Nz

    Merzbild.compute_velocity_extent!(1, 1, mg, phys_props, species_list)

    for i in 1:3
        @test abs(mg.extent_v_lower[i] - (-1)) <= eps(Float64)
        @test abs(mg.extent_v_upper[i] - (1)) <= eps(Float64)
    end
    
    counter = 0
    for i in 1:Nx
        for j in 1:Ny
            for k in 1:Nz
                counter += 1
                v = [mg.extent_v_lower[1] + 1e-5 + (i-1) * mg.Δv[1],
                     mg.extent_v_lower[2] + 1e-5 + (j-1) * mg.Δv[2],
                     mg.extent_v_lower[3] + 1e-5 + (k-1) * mg.Δv[3]]
                index = Merzbild.compute_grid_index(mg, v)
                @test index == counter
            end
        end
    end
 
    v0 = [-2.0, -0.5, -0.5]
    index = Merzbild.compute_grid_index(mg, v0)
    @test index == mg.Ntotal + 1
    v0 = [2.0, -0.5, -0.5]
    index = Merzbild.compute_grid_index(mg, v0)
    @test index == mg.Ntotal + 2  # y, z < middle

    v0 = [-2.0, 0.5, -0.5]
    index = Merzbild.compute_grid_index(mg, v0)
    @test index == mg.Ntotal + 3
    v0 = [2.0, 1.5, -0.5]
    index = Merzbild.compute_grid_index(mg, v0)
    @test index == mg.Ntotal + 4  # y > middle, z < middle

    v0 = [-2.0, -0.5, 3.5]
    index = Merzbild.compute_grid_index(mg, v0)
    @test index == mg.Ntotal + 5
    v0 = [2.0, -0.5, 10.5]
    index = Merzbild.compute_grid_index(mg, v0)
    @test index == mg.Ntotal + 6  # y < middle, z > middle

    v0 = [-1.1, 0.5, 0.7]
    index = Merzbild.compute_grid_index(mg, v0)
    @test index == mg.Ntotal + 7
    v0 = [30.0, 0.7, 0.9]
    index = Merzbild.compute_grid_index(mg, v0)
    @test index == mg.Ntotal + 8  # y, z > middle

    Δabs = 2.5
    Δrel_xsmall = 5e-13

    n_particles = 40000
    particles::Vector{Vector{Particle}} = [Vector{Particle}(undef, n_particles)]
    T0 = 300.0
    Fnum = 1e20
    vx0 = 2000.0
    vy0 = 500.0
    vz0 = -400.0
    sample_particles_equal_weight!(rng, particles[1], n_particles, T0, species_list[1].mass, Fnum, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
    vx0=vx0, vy0=vy0, vz0=vz0)


    particle_indexer::Array{ParticleIndexer, 2} = Array{ParticleIndexer, 2}(undef, 1, 1)
    particle_indexer[1,1] = create_particle_indexer(n_particles)

    compute_props!(phys_props, particle_indexer, particles, species_list)
    mg2 = create_merging_grid(Nx, Ny, Nz, 3.5)

    Merzbild.compute_velocity_extent!(1, 1, mg2, phys_props, species_list)

    @test abs(vx0 - mg2.extent_v_mid[1]) <= Δabs
    @test abs(vy0 - mg2.extent_v_mid[2]) <= Δabs
    @test abs(vz0 - mg2.extent_v_mid[3]) <= Δabs
    Merzbild.compute_grid!(1, 1, mg2, particles, particle_indexer)

    ntot = 0.0
    nptot = 0
    for i in 1:mg2.Ntotal
        ntot += mg2.cells[i].w
        nptot += mg2.cells[i].np
    end

    @test nptot == n_particles
    @test abs((ntot - n_particles * Fnum)) / (n_particles * Fnum) < Δrel_xsmall

    mg3 = create_merging_grid(1, 1, 1, 500.0)
    Merzbild.compute_velocity_extent!(1, 1, mg3, phys_props, species_list)
    Merzbild.compute_grid!(1, 1, mg3, particles, particle_indexer)

    @test mg3.cells[1].np == n_particles
    @test abs((mg3.cells[1].w - n_particles * Fnum)) / (n_particles * Fnum) < Δrel_xsmall

    @test abs(phys_props.v[1,1,1] - mg3.cells[1].v_mean[1]) <= eps(Float64)
    @test abs(phys_props.v[2,1,1] - mg3.cells[1].v_mean[2]) <= eps(Float64)
    @test abs(phys_props.v[3,1,1] - mg3.cells[1].v_mean[3]) <= eps(Float64)

    @test mg3.cells[1].x_mean[1] <= 1.0
    @test mg3.cells[1].x_mean[1] >= 0.0

    @test mg3.cells[1].x_mean[2] <= 1.0
    @test mg3.cells[1].x_mean[2] >= 0.0

    @test mg3.cells[1].x_mean[3] <= 1.0
    @test mg3.cells[1].x_mean[3] >= 0.0

    etot = (mg3.cells[1].v_std_sq[1] + mg3.cells[1].v_std_sq[2] + mg3.cells[1].v_std_sq[3]) * species_list[1].mass / (3 * k_B)

    @test abs(etot - phys_props.T[1,1]) / phys_props.T[1,1] <= Δrel_xsmall
end