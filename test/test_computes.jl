@testset "phys_props compute" begin

    species_list::Vector{Species} = load_species_list("data/particles.toml", "Ar")

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

    particle_indexer = create_particle_indexer_array(n_particles)

    phys_props::PhysProps = create_props(1, 1, [], Tref=1)
    phys_props_no_moments::PhysProps = create_props(1, 1, [], Tref=1)
    compute_props!(phys_props, particle_indexer, particles, species_list)
    @test abs((phys_props.n[1,1] - n_dens) / n_dens) < Δsmall
    @test abs((phys_props.v[1,1,1] - v0[1])) < Δsmall
    @test abs((phys_props.v[2,1,1] - v0[2])) < Δsmall
    @test abs((phys_props.v[3,1,1] - v0[3])) < Δsmall
    @test abs((phys_props.T[1,1])) < eps()

    compute_props_sorted_without_moments!(phys_props_no_moments, particle_indexer, particles, species_list)
    @test abs(phys_props_no_moments.n[1,1] - phys_props.n[1,1]) < eps()
    @test abs((phys_props_no_moments.v[1,1,1] - phys_props.v[1,1,1])) < Δsmall
    @test abs((phys_props_no_moments.v[2,1,1] - phys_props.v[2,1,1])) < Δsmall
    @test abs((phys_props_no_moments.v[3,1,1] - phys_props.v[3,1,1])) < Δsmall
    @test abs(phys_props_no_moments.T[1,1]) < eps()

    mom000 = Merzbild.compute_mixed_moment(particle_indexer, particles, 1, 1, [0, 0, 0])
    n0 = phys_props.n[1,1]
    @test abs(mom000 - phys_props.n[1,1]) < eps()
    mom100 = Merzbild.compute_mixed_moment(particle_indexer, particles, 1, 1, [1, 0, 0])
    @test abs((mom100 - phys_props.v[1,1,1] * n0)) < Δsmall
    mom010 = Merzbild.compute_mixed_moment(particle_indexer, particles, 1, 1, [0, 1, 0])
    @test abs((mom010 - phys_props.v[2,1,1] * n0)) < Δsmall
    mom001 = Merzbild.compute_mixed_moment(particle_indexer, particles, 1, 1, [0, 0, 1])
    @test abs((mom001 - phys_props.v[3,1,1] * n0)) < Δsmall

    mom110 = Merzbild.compute_mixed_moment(particle_indexer, particles, 1, 1, [1, 1, 0])
    mom110 /= n0
    @test abs(mom110 - (-2)) < Δverysmall
    
    mom101 = Merzbild.compute_mixed_moment(particle_indexer, particles, 1, 1, [1, 0, 1])
    mom101 /= n0
    @test abs(mom101 - (4)) < Δverysmall

    mom011 = Merzbild.compute_mixed_moment(particle_indexer, particles, 1, 1, [0, 1, 1])
    mom011 /= n0
    @test abs(mom011 - (-8)) < Δverysmall

    sum_scaler = 1.0 / n0
    mom132 = Merzbild.compute_mixed_moment(particle_indexer, particles, 1, 1, [1, 3, 2]; sum_scaler=sum_scaler)
    @test abs(mom132 - (-128)) < Δverysmall
end