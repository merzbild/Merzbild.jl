# Test that sampling from a Maxwellian and computing physical properties gives reasonable results

include("../../../src/Merzbild.jl")

using ..Merzbild
using Random

function run(seed)
    species_data = load_species_data("data/particles.toml", "Ar")
    println([species.name for species in species_data])
    rng = Xoshiro(seed)

    n_particles = 5000

    particles = [ParticleVector(n_particles)]
    pia = ParticleIndexerArray(0)

    sample_particles_equal_weight!(rng, particles[1], pia, 1, 1,
                                   n_particles, species_data[1].mass, 500.0, 1e10, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
    
    phys_props = PhysProps(1, 1, [])
    compute_props!(particles, pia, species_data, phys_props)
    println(phys_props.n)
    println(phys_props.v)
    println(phys_props.T)

    ds = NCDataHolder("scratch/data/test_sampling_io.nc", species_data, phys_props)
    write_netcdf(ds, phys_props, 1)
    close_netcdf(ds)
end

run(1234)