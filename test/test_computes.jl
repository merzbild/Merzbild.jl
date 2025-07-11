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

    particles::Vector{Vector{Particle}} = [Vector{Particle}(undef, n_particles)]

    for i in 1:n_particles
        particles[1][i] = Particle(Fnum, v0, x0)
    end

    pia = ParticleIndexerArray(n_particles)

    phys_props::PhysProps = PhysProps(1, 1, [], Tref=1)
    phys_props_no_moments::PhysProps = PhysProps(1, 1, [], Tref=1)
    compute_props!(particles, pia, species_data, phys_props)
    @test abs((phys_props.n[1,1] - n_dens) / n_dens) < Δsmall
    @test abs((phys_props.v[1,1,1] - v0[1])) < Δsmall
    @test abs((phys_props.v[2,1,1] - v0[2])) < Δsmall
    @test abs((phys_props.v[3,1,1] - v0[3])) < Δsmall
    @test abs((phys_props.T[1,1])) < eps()

    compute_props_sorted!(particles, pia, species_data, phys_props_no_moments)
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
    phys_props2 = PhysProps(pia, [2, 3, 4], Tref=600.0)
    @test phys_props2.n_cells == 1
    @test phys_props2.n_species == 1
    @test phys_props2.Tref == 600.0
    @test phys_props2.moment_powers == [2,3,4]

    pia2 = ParticleIndexerArray(40, 3)
    phys_props3 = PhysProps(pia2)
    @test phys_props3.n_cells == 40
    @test phys_props3.n_species == 3
    @test phys_props3.Tref == 300.0
    @test phys_props3.moment_powers == []

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

    phys_props_ndens = PhysProps(pia, [2, 3, 4], Tref=600.0; ndens_not_Np=true)
    @test_throws ErrorException avg_props!(phys_props_avg, phys_props_ndens, 2)
end