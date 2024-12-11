# Test that sampling two different species from a Maxwellian, doing multiple timesteps of collisions and computing physical properties gives reasonable results

# SPARTA: 12k particles 13.29s user 0.08s system 98% cpu 13.539 total
# this code: 6.90s user 1.78s system 127% cpu 6.829 total
# julia --project=. simulations/basic/basic_collisions_multispecies.jl  7.47s user 1.53s system 126% cpu 7.142 total
# julia --project=. simulations/basic/basic_collisions_multispecies.jl  7.28s user 1.43s system 127% cpu 6.839 total
# julia --project=. simulations/basic/basic_collisions_multispecies.jl  4.92s user 1.57s system 140% cpu 4.634 total
# without moments
# julia --project=. simulations/basic/basic_collisions_multispecies.jl  4.80s user 1.74s system 139% cpu 4.677 total

include("../../../src/Merzbild.jl")

using ..Merzbild
using Random

function run(seed)
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    species_data::Vector{Species} = load_species_data("data/particles.toml", ["Ar", "He"])
    interaction_data::Array{Interaction, 2} = load_interaction_data("data/vhs.toml", species_data)
    n_species = length(species_data)

    println([species.name for species in species_data])
    println(interaction_data)

    n_t = 5000

    n_particles_Ar = 2000
    n_particles_He = 10000
    T0_Ar = 300.0
    T0_He = 2000.0
    T0_list = [T0_Ar, T0_He]
    Fnum = 5e11

    particles::Vector{ParticleVector} = [ParticleVector(n_particles_Ar), ParticleVector(n_particles_He)]
    pia = ParticleIndexerArray([0, 0])

    sample_particles_equal_weight!(rng, particles[1], pia, 1, 1, n_particles_Ar, species_data[1].mass, T0_Ar, Fnum, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
    sample_particles_equal_weight!(rng, particles[2], pia, 1, 2, n_particles_He, species_data[2].mass, T0_He, Fnum, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)

    phys_props::PhysProps = PhysProps(1, 2, [], Tref=T0_Ar)
    compute_props!(particles, pia, species_data, phys_props)
    println(phys_props.n)
    println(phys_props.v)
    println(phys_props.T)

    ds = NCDataHolder("scratch/data/2species.nc", species_data, phys_props)
    write_netcdf_phys_props(ds, phys_props, 0)

    collision_factors::Array{CollisionFactors, 3} = create_collision_factors_array(n_species)
    collision_data::CollisionData = CollisionData()

    estimate_sigma_g_w_max!(collision_factors, interaction_data, species_data, T0_list, Fnum)

    Δt = 2.5e-3
    V = 1.0

    @time for ts in 1:n_t
        for s2 in 1:n_species
            for s1 in s2:n_species
                if (s1 == s2)
                    ntc!(rng, collision_factors[s1,s1,1], collision_data, interaction_data, particles[s1], pia, 1, s1, Δt, V)
                else
                    ntc!(rng, collision_factors[s1,s2,1], collision_data, interaction_data, particles[s1], particles[s2],
                         pia, 1, s1, s2, Δt, V)
                end
            end
        end

        compute_props_sorted!(particles, pia, species_data, phys_props)
        write_netcdf_phys_props(ds, phys_props, ts)
    end
    println(phys_props.n)
    println(phys_props.v)
    println(phys_props.T)

    close_netcdf(ds)
end

run(1234)