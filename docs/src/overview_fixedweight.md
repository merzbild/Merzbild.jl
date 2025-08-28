# Fixed-weight DSMC simulations

Now that we can load species data, sample and index particles and compute macroscopic properties,
we now look at performing collisions. Currently, the No Time Counter (NTC) algorithm of Bird
is implemented for fixed-weight DSMC simulations.

## Loading interaction data: Interaction
First, we need to load interaction data for the species the collisions of which we want to model.
This data includes things like the collision-reduced mass and VHS parameters. The VHS
parameters are collision species pair-specific. These data are stored in an `Interaction` instance.
These `Interaction` instances are stored in a `n_species x n_species` matrix, where
the `(i,j)`-th element contains the interaction data for collisions of particles of species `i` with 
particles of species `j`.

We can use the [`load_interaction_data`](@ref) function to create such a matrix of `Interaction` instances.
It which reads a TOML file with the relevant interaction information, and uses the (already loaded)
species data to compute quantities such as the collision-reduced mass.
One can also use [`load_interaction_data_with_dummy`](@ref) function, which will not throw an error
in case data for a specific interaction
is missing in the TOML file, but will just create an interaction using the passed dummy
parameters. This is relevant for electron-neutral interactions, for example, since 
VHS collision parameters don't really make sense for such interactions, but are required to fill in the fields.

An additional utility function `load_species_and_interaction_data` is also available, which loads both
the species' and interaction data for those species at the same time.

## Storing intermediate collision data: CollisionData
During collisions, multiple vector quantities (such as the center of mass velocity,
pre- and post-collision relative velocity) are computed for multiple collision pairs. In order to avoid excessive allocations
and keep track of these quantities, the `CollisionData` structure exists.
It does not need to be initialized in any special way; an instance of `CollisionData` must be created and simply
passed to the collision routine, which will use it as required.

## NTC-specific data: CollisionFactors
The No-Time Counter (NTC) collision algorithm requires an estimate of ``(\sigma g w)_{max}`` for each species pair
for each cell in the flow. The `CollisionFactors` data structure stores the values of this factor, as well as 
other NTC-related parameters (number of collisions performed, number of collision partners).
One needs to initialize a 3-dimensional array of `CollisionFactors` of shape 
`n_species x n_species x n_cells` in order to store the required factors for each species' pair for all cells in the flow;
this is done by calling `create_collision_factors_array(n_species, n_cells)`.
Some initial estimate for the values of ``(\sigma g w)_{max}``; the simplest way to do so for a fixed-weight DSMC simulation
is to call the `estimate_sigma_g_w_max!(collision_factors, interactions, species_data, T_list, Fnum; mult_factor=1.0)`
function.
Version of the `create_collision_factors_array` function that automatically estimate ``(\sigma g w)_{max}`` are also
available.

This will pre-compute  ``(\sigma g w)_{max}`` based on a list of temperatures for each species by computing the
average thermal velocities and calculating the value of the VHS collision cross-section:
```julia
g_thermal1 = sqrt(2 * T1 * k_B / species1.mass)
g_thermal2 = sqrt(2 * T2 * k_B / species2.mass)
g_thermal = 0.5 * (g_thermal1 + g_thermal2)
return mult_factor * sigma_vhs(interaction, g_thermal) * g_thermal * Fnum
```
Here `mult_factor` (default value of 1.0) is used to optionally increase or decrease the computed estimate; and the value ``w``
is assumed to be constant and is designated by `Fnum`. This approach computes the same value of ``(\sigma g w)_{max}`` for all
the cells in the domain; however, unless an initial condition with large gradients is used, this should not be an issue.

## NTC collisions
Now that the interaction data has been loaded, the `CollisionData` instance to store intermediate collision-related quantities
has been instantiated, the 3-dimensional array of `CollisionFactors` has been created, and the ``(\sigma g w)_{max}``
values precomputed, one can perform collisions.

Single-species elastic VHS collisions can be performed by calling
`ntc!(rng, collision_factors, collision_data, interaction, particles, pia, cell, species, Δt, V)`.
Here `collision_factors` is the specific instance of `CollisionFactors`, i.e. a specific element of the
3-dimensional array of `CollisionFactors`. `Δt` is the timestep, and `V` is the volume of the physical cell
(for spatially homogeneous simulations this can be set to 1.0).

Multi-species elastic VHS collisions are performed in a similar fashion, by calling
`ntc!(rng, collision_factors, collision_data, interaction, particles_1, particles_2, pia, cell, species1, species2, Δt, V)`.

## Example: bringing it all together
An example of computation of collisions for a two-species mixture is presented here.
```julia
using Merzbild
using Random

seed = 1
Random.seed!(seed)
rng = Xoshiro(seed)

species_data = load_species_data("data/particles.toml", ["Ar", "He"])
n_species = length(species_data)

# number of timesteps to run simulation for
n_t = 5000

# we will have the number density of argon = 1e15, the number density of helium = 5e15
n_particles_Ar = 2000
n_particles_He = 10000
Fnum = 5e11

# initial temperatures
T0_Ar = 300.0
T0_He = 2000.0
T0_list = [T0_Ar, T0_He]

# create the 2-element Vector of ParticleVectors for the 2 species
particles = [ParticleVector(n_particles_Ar), ParticleVector(n_particles_He)]

# create the 2-species 1-cell particle indexer array filled with zeros
# as we haven't sampled any particles yet
pia = ParticleIndexerArray([0, 0])

# sample particles from a Maxwellian distribution
sample_particles_equal_weight!(rng, particles[1], pia, 1, 2, n_particles_Ar,
                               species_data[1].mass, T0_Ar, Fnum,
                               0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
sample_particles_equal_weight!(rng, particles[2], pia, 1, 2, n_particles_He,
                               species_data[2].mass, T0_He, Fnum,
                               0.0, 1.0, 0.0, 1.0, 0.0, 1.0)


# create the PhysProps instance to store computed properties
phys_props = PhysProps(1, 2, [], Tref=T0_Ar)

# create struct for I/O
ds = NCDataHolder("2species.nc", species_data, phys_props)

# compute and write physical properties at t=0
compute_props!(particles, pia, species_data, phys_props)
write_netcdf(ds, phys_props, 0)

# load interaction data
interaction_data = load_interaction_data("data/vhs.toml", species_data)

# create the 3-D array of collision factors
collision_factors::Array{CollisionFactors, 3} = create_collision_factors_array(n_species)

# create the structure to store temporary collision data
collision_data = CollisionData()

# estimate  (σ(g) * g * w)_max
estimate_sigma_g_w_max!(collision_factors, interaction_data, species_data, T0_list, Fnum)

# set our timestep
Δt = 2.5e-3

# set cell volume to 1.0 as we're doing a 0-D simulation
V = 1.0

for ts in 1:n_t  # loop over time
    for s2 in 1:n_species  # loop over first species
        for s1 in s2:n_species  # loop over second species
            if (s1 == s2)
                # collisions between particles of same species
                ntc!(rng, collision_factors[s1,s1,1], collision_data,
                     interaction_data, particles[s1], pia, 1, s1, Δt, V)
            else
                # collisions between particles of different species
                ntc!(rng, collision_factors[s1,s2,1], collision_data,
                     interaction_data, particles[s1], particles[s2],
                     pia, 1, s1, s2, Δt, V)
            end
        end
    end

    # we don't do convection, so we can take advantage of the fact that the particles stay sorted
    # and compute the physical properties slightly faster without having to sort the particles
    compute_props_sorted!(particles, pia, species_data, phys_props)

    # write to output
    write_netcdf(ds, phys_props, ts)
end

# close output file
close_netcdf(ds)
```

## Summary
Now we have an overview of how to

1. Load interaction-related data
2. Estimate the factor required for the NTC collision algorithm
3. Perform elastic collisions between particles of alike and unalike species

In the next section, an overview of simulating collisions with variable-weight particles will be given.