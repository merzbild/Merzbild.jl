@testset "phys_props compute" begin

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, "Ar")

    n_particles = 2000

    n_dens = 1e25
    Δsmall = 1e-8
    Δverysmall = 1e-11

    Fnum::Float64 = n_dens / n_particles

    v0 = [-1.0, 2.0, -4.0]
    x0 = [10.0, 20.0, 30.0]

    particles::Vector{ParticleVector} = [ParticleVector(n_particles)]

    for i in 1:n_particles
        particles[1][i] = Particle(Fnum, v0, x0)
    end

    pia = ParticleIndexerArray(n_particles)

    phys_props::PhysProps = PhysProps(1, 1)
    phys_props_no_moments::PhysProps = PhysProps(1, 1)
    compute_props!(particles, pia, species_data, phys_props)
    @test abs((phys_props.np[1,1] - n_particles) / n_particles) <= eps()
    @test abs((phys_props.n[1,1] - n_dens) / n_dens) < Δsmall
    @test abs((phys_props.v[1,1,1] - v0[1])) < Δsmall
    @test abs((phys_props.v[2,1,1] - v0[2])) < Δsmall
    @test abs((phys_props.v[3,1,1] - v0[3])) < Δsmall
    @test abs((phys_props.T[1,1])) < eps()

    compute_props_sorted!(particles, pia, species_data, phys_props_no_moments)
    @test abs(phys_props_no_moments.np[1,1] - phys_props.np[1,1]) < eps()
    @test abs(phys_props_no_moments.n[1,1] - phys_props.n[1,1]) < eps()
    @test abs((phys_props_no_moments.v[1,1,1] - phys_props.v[1,1,1])) < Δsmall
    @test abs((phys_props_no_moments.v[2,1,1] - phys_props.v[2,1,1])) < Δsmall
    @test abs((phys_props_no_moments.v[3,1,1] - phys_props.v[3,1,1])) < Δsmall
    @test abs(phys_props_no_moments.T[1,1]) < eps()

    mom000 = Merzbild.compute_mixed_moment(particles, pia, 1, 1, [0, 0, 0])
    n0 = phys_props.n[1,1]
    @test abs(mom000 - phys_props.n[1,1]) < eps()
    mom100 = Merzbild.compute_mixed_moment(particles, pia, 1, 1, [1, 0, 0])
    @test abs((mom100 - phys_props.v[1,1,1] * n0)) < Δsmall
    mom010 = Merzbild.compute_mixed_moment(particles, pia, 1, 1, [0, 1, 0])
    @test abs((mom010 - phys_props.v[2,1,1] * n0)) < Δsmall
    mom001 = Merzbild.compute_mixed_moment(particles, pia, 1, 1, [0, 0, 1])
    @test abs((mom001 - phys_props.v[3,1,1] * n0)) < Δsmall

    mom110 = Merzbild.compute_mixed_moment(particles, pia, 1, 1, [1, 1, 0])
    mom110 /= n0
    @test abs(mom110 - (-2)) < Δverysmall
    
    mom101 = Merzbild.compute_mixed_moment(particles, pia, 1, 1, [1, 0, 1])
    mom101 /= n0
    @test abs(mom101 - (4)) < Δverysmall

    mom011 = Merzbild.compute_mixed_moment(particles, pia, 1, 1, [0, 1, 1])
    mom011 /= n0
    @test abs(mom011 - (-8)) < Δverysmall

    sum_scaler = 1.0 / n0
    mom132 = Merzbild.compute_mixed_moment(particles, pia, 1, 1, [1, 3, 2]; sum_scaler=sum_scaler)
    @test abs(mom132 - (-128)) < Δverysmall

    # test different constructors of PhysProps
    phys_props2 = PhysProps(pia)
    @test phys_props2.n_cells == 1
    @test phys_props2.n_species == 1

    pia2 = ParticleIndexerArray(40, 3)
    phys_props3 = PhysProps(pia2)
    @test phys_props3.n_cells == 40
    @test phys_props3.n_species == 3

    # test averaging of PhysProps
    phys_props4 = PhysProps(pia2)
    phys_props_avg = PhysProps(pia2)

    phys_props3.T .= 100.0
    phys_props3.n .= 2e11
    phys_props3.np .= 100
    phys_props3.lpa .= 1e7
    phys_props3.v .= 1e3

    phys_props4.T .= 200.0
    phys_props4.n .= 4e11
    phys_props4.np .= 200
    phys_props4.lpa .= 2e7
    phys_props4.v .= 2e3

    avg_props!(phys_props_avg, phys_props3, 2)
    avg_props!(phys_props_avg, phys_props4, 2)

    @test maximum(abs.(phys_props_avg.T .- 150.0)) < 2 * eps()
    @test maximum(abs.(phys_props_avg.n .- 3e11)) < 2 * eps()
    @test maximum(abs.(phys_props_avg.np .- 150)) < 2 * eps()
    @test maximum(abs.(phys_props_avg.lpa .- 1.5e7)) < 2 * eps()
    @test maximum(abs.(phys_props_avg.v .- 1.5e3)) < 2 * eps()

    clear_props!(phys_props_avg)
    @test maximum(abs.(phys_props_avg.T)) == 0.0
    @test maximum(abs.(phys_props_avg.n)) == 0.0
    @test maximum(abs.(phys_props_avg.np)) == 0.0
    @test maximum(abs.(phys_props_avg.lpa)) == 0.0
    @test maximum(abs.(phys_props_avg.v)) == 0.0

    phys_props_ndens = PhysProps(pia; ndens_not_Np=true)
    @test_throws ErrorException avg_props!(phys_props_avg, phys_props_ndens, 2)
end

@testset "phys_props compute bias correction" begin
    rng = StableRNG(1234)

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data::Vector{Species} = load_species_data(particles_data_path, "Ar")

    T0 = 500.0
    Fnum = 1e10
    vscale = sqrt(k_B * T0 / species_data[1].mass)

    ppc = 3
    n_cells = 20000
    n_particles = ppc * n_cells
    particles = [ParticleVector(n_particles)]
    pia = ParticleIndexerArray(n_cells, 1)
    for cell in 1:n_cells
        pia.indexer[cell, 1] = ParticleIndexer(ppc, ppc * (cell - 1) + 1, ppc * cell, ppc, 0, -1, 0)
        for i in ppc * (cell - 1) + 1:ppc * cell
            particles[1][i] = Particle(Fnum * (1.0 + 4.0 * rand(rng)),
                                       [100.0 + vscale * randn(rng), -50.0 + vscale * randn(rng), vscale * randn(rng)],
                                       [0.0, 0.0, 0.0])
        end
    end
    pia.n_total[1] = n_particles

    pp_biased = PhysProps(pia)
    pp_corr = PhysProps(pia)
    compute_props!(particles, pia, species_data, pp_biased)
    compute_props!(particles, pia, species_data, pp_corr; bias_correction=true)

    max_ratio_err = 0.0
    for cell in 1:n_cells
        n = 0.0
        w2 = 0.0
        for i in ppc * (cell - 1) + 1:ppc * cell
            n += particles[1][i].w
            w2 += particles[1][i].w^2
        end
        ratio = pp_corr.T[cell, 1] / pp_biased.T[cell, 1]
        max_ratio_err = max(max_ratio_err, abs(ratio - n / (n - w2 / n)))
    end
    @test max_ratio_err < 1e-12
    @test maximum(abs.(pp_corr.n .- pp_biased.n)) == 0.0
    @test maximum(abs.(pp_corr.v .- pp_biased.v)) == 0.0
    @test maximum(abs.(pp_corr.np .- pp_biased.np)) == 0.0

    @test abs(sum(pp_corr.T) / n_cells - T0) / T0 < 0.02
    @test sum(pp_biased.T) / n_cells < 0.8 * T0

    pp_sorted = PhysProps(pia)
    compute_props_sorted!(particles, pia, species_data, pp_sorted; bias_correction=true)
    @test maximum(abs.(pp_sorted.T .- pp_corr.T)) < 1e-10

    pp_chunk = PhysProps(pia)
    compute_props_sorted!(particles, pia, species_data, pp_chunk, 1:2; bias_correction=true)
    @test maximum(abs.(pp_chunk.T[1:2, 1] .- pp_corr.T[1:2, 1])) < 1e-10
    @test maximum(abs.(pp_chunk.T[3:end, 1])) == 0.0

    grid = Grid1DUniform(2.0, n_cells)
    pp_grid = PhysProps(pia)
    compute_props_sorted!(particles, pia, species_data, pp_grid, grid; bias_correction=true)
    @test maximum(abs.(pp_grid.T .- pp_corr.T)) < 1e-10

    pp_grid_chunk = PhysProps(pia)
    compute_props_sorted!(particles, pia, species_data, pp_grid_chunk, grid, 1:2; bias_correction=true)
    @test maximum(abs.(pp_grid_chunk.T[1:2, 1] .- pp_corr.T[1:2, 1])) < 1e-10
    @test maximum(abs.(pp_grid_chunk.T[3:end, 1])) == 0.0

    pp_grid_ndens = PhysProps(pia; ndens_not_Np=true)
    compute_props_sorted!(particles, pia, species_data, pp_grid_ndens, grid; bias_correction=true)
    @test maximum(abs.(pp_grid_ndens.T .- pp_corr.T)) < 1e-10

    # equal weights: correction factor is Np / (Np - 1)
    n_eq = 10
    particles_eq = [ParticleVector(n_eq)]
    for i in 1:n_eq
        particles_eq[1][i] = Particle(Fnum, [vscale * randn(rng), vscale * randn(rng), vscale * randn(rng)], [0.0, 0.0, 0.0])
    end
    pia_eq = ParticleIndexerArray(n_eq)
    pp_eq_biased = PhysProps(pia_eq)
    pp_eq_corr = PhysProps(pia_eq)
    compute_props!(particles_eq, pia_eq, species_data, pp_eq_biased)
    compute_props!(particles_eq, pia_eq, species_data, pp_eq_corr; bias_correction=true)
    @test abs(pp_eq_corr.T[1, 1] / pp_eq_biased.T[1, 1] - n_eq / (n_eq - 1)) < 1e-12

    # single particle: temperature is not defined, set to 0
    pia_single = ParticleIndexerArray(1)
    pp_single = PhysProps(pia_single)
    compute_props!(particles_eq, pia_single, species_data, pp_single; bias_correction=true)
    @test pp_single.T[1, 1] == 0.0
    @test pp_single.n[1, 1] == Fnum
    compute_props_sorted!(particles_eq, pia_single, species_data, pp_single; bias_correction=true)
    @test pp_single.T[1, 1] == 0.0
end