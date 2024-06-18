# Test that sampling from a Maxwellian, doing multiple timesteps of collisions and computing physical properties gives reasonable results

include("../../src/merzbild.jl")

using ..Merzbild
using Random
using TimerOutputs

function run(n_t, n_particles, seed)
    Random.seed!(seed)
    rng = Xoshiro(seed)

    species_list = load_species_list("data/particles.toml", "Ar")
    interaction_data = load_interaction_data("data/vhs.toml", species_list)

    reset_timer!()
    # println([species.name for species in species_list])
    # println(interaction_data)

    # n_particles = 5000
    T0 = 500.0
    ndens = 5e20
    Fnum = ndens / n_particles
    println(Fnum)

    particles = [Vector{Particle}(undef, n_particles)]
    sample_particles_equal_weight!(rng, particles[1], n_particles, T0, species_list[1].mass, Fnum, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)

    pia = create_particle_indexer_array(n_particles)

    phys_props = create_props(1, 1, [], Tref=T0)
    compute_props!(phys_props, pia, particles, species_list)
    # println(phys_props.n)
    # println(phys_props.v)
    # println(phys_props.T)

    ds = create_netcdf_phys_props("test_maxwellian_$n_particles.nc", phys_props, species_list)
    write_netcdf_phys_props(ds, phys_props, 0)

    collision_factors = create_collision_factors()
    collision_data = create_collision_data()

    collision_factors.sigma_g_w_max = estimate_sigma_g_w_max(interaction_data[1,1], species_list[1], T0, Fnum)

    Δt = 1e-6
    V = 1.0

    @timeit "outer loop" for ts in 1:n_t
        @timeit "collisions" ntc!(1, 1, rng, collision_factors, pia, collision_data, interaction_data[1,1], particles[1],
            Δt, V)
            # ntc!(1, 1, rng, collision_factors, pia, collision_data, interaction_data[1,1], particles[1],
            # Δt, V)
        # println(collision_factors)
        
        @timeit "props" compute_props!(phys_props, pia, particles, species_list)
        @timeit "I/O" write_netcdf_phys_props(ds, phys_props, ts)
    end
    # println(phys_props.n)
    # println(phys_props.v)
    # println(phys_props.T)
    close(ds)
    print_timer()
end

run(1, 100, 1234)
# run(10000, 100, 1234)
# run(10000, 1000, 1234)
# run(10000, 50000, 1234)
run(10000, 200000, 1234)