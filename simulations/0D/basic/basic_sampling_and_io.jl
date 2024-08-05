# Test that sampling from a Maxwellian and computing physical properties gives reasonable results

include("../../../src/Merzbild.jl")

using ..Merzbild
using Random

function run(seed)
    species_list = load_species_list("data/particles.toml", "Ar")
    println([species.name for species in species_list])
    rng = Xoshiro(seed)

    n_particles = 5000

    particles = [Vector{Particle}(undef, n_particles)]
    sample_particles_equal_weight!(rng, particles[1], n_particles, species_list[1].mass, 500.0, 1e10, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
    
    pia = create_particle_indexer_array(n_particles)

    phys_props = create_props(1, 1, [])
    compute_props!(particles, pia, species_list, phys_props)
    println(phys_props.n)
    println(phys_props.v)
    println(phys_props.T)

    ds = create_netcdf_phys_props("scratch/data/test_sampling_io.nc", species_list, phys_props)
    write_netcdf_phys_props(ds, phys_props, 1)
    close_netcdf(ds)
end

run(1234)