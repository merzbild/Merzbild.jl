# Variable-weight DSMC simulations

For variable-weight DSMC simulation not much is different compared to fixed-weight DSMC simulations.
Existing routines for computation of macroscopic properties, I/O, and even collisions can be re-used, since
they all do not explicitly assume fixed particle weights.

Two additional things are required, however: a way to sample particles with variable weights, and
an algorithm to merge particles to avoid a blow-up of the number of particles due to collisions.

## Sampling variable-weight particles
In order to obtain variable-weight particles, one can either use the [`sample_particles_equal_weight!`](@ref)
multiple times with different `Fnum` values, or one can use the 
[`sample_on_grid!`](@ref) function and its more specific version [`sample_maxwellian_on_grid!`](@ref), as mentioned
in the [Overview of basic building blocks](@ref) section.

These approaches compute the values of a VDF function on a discrete velocity grid and then re-interpret the grid points
and associated VDF values as particles. Thus one can get a good resolution of high-velocity tails (with the particles there
having small weights). The function constructs the velocity grid via a Cartesian product of discrete grids for each velocity direction.
In each direction, the interval ``[-v_m\sqrt{2kT/m},v_m\sqrt{2kT/m}]`` is discretized into `nv` subintervals. Here
the multiplier ``v_m`` is specified by the `v_mult` parameter; so the discrete velocity grid spans several mean thermal velocities.
In addition, only particles with a speed less than ``c_m \sqrt{2kT/m}`` are actually created, where ``c_m`` is specified by
the `cutoff_mult` parameter. Finally, numerical noise may be added to the velocity via the `noise` parameter,
and a velocity offset via the `v_offset` parameter; the spatial distribution of the particles is randomized in the cell.
The function returns the total number of particles sampled (which may be less than `nv^3` depending on the choice of the value
of `cutoff_mult`).

An example of such sampling is shown below:
```julia
using Merzbild
using Random

seed = 1
Random.seed!(seed)
rng = Xoshiro(seed)

# our number density
ndens = 1e23 

# our temperature
T0 = 300.0  

# load species and interaction data
species_data = load_species_data("data/particles.toml", "Ar")
interaction_data = load_interaction_data("data/pseudo_maxwell.toml", species_data)

# number of velocity grid points in each velocity direction
nv = 20  

# some initial guess on # of particles in simulation
np_base = 40^3  

particles = [ParticleVector(np_base)]

# sample from a BKW distribution at t=0
vdf0 = (vx, vy, vz) -> bkw(vx, vy, vz, species_data[1].mass, T0, 0.0)

# returns number of particles sampled
n_sampled = sample_on_grid!(rng, vdf0, particles[1], nv, species_data[1].mass, T0, ndens,
                            0.0, 1.0, 0.0, 1.0, 0.0, 1.0;
                            v_mult=3.5, cutoff_mult=3.5, noise=0.0, v_offset=[0.0, 0.0, 0.0])

# create the ParticleIndexerArray
pia = ParticleIndexerArray(n_sampled)
```

## Colliding variable-weight particles
As mentioned above, the existing DSMC routines take care of the required particle splitting, so no specific adaptation is required
for the variable-weight case. However, one needs to estimate ``(\sigma g w)_{max}``. The simplest approach is to
use `estimate_sigma_g_w_max` function, but it requires a (fixed) value of `Fnum` - so one can simply compute it as
`ndens/n_sampled` (i.e. the number density divided by the total number of sampled particles).

## Merging variable-weight particles
In order to merge variable-weight particles, one needs to set up a merging algorithm and the associated
data structures and parameters. A more detailed overview of merging algorithms will appear later in a Tutorials
section; here, the Octree merging approach of [Martin and Cambier (2016)](https://doi.org/10.1016/j.jcp.2016.01.020)
is used. The algorithm groups particles into bins in velocity space, recursively refining the bins until
the number of post-merge particles reaches the prescribed value. It then performs an ``N:2`` merge in each bin,
replacing all particles in a bin with 2 particles.

To set up this merging algorithm, one needs to create a `OctreeN2Merge` instance, specifying
- The extent of the first (root) bin - whether it accounts for the extent of the particles or whether it is just taken to be very large
- How bins are split (along the middle velocity, the mean velocity, or the median velocity)
- Whether the velocity bounds of each sub-bin are recomputed based on the particles in the bin or are based purely on the bounds of the parent bin and the splitting velocity
- Maximum number of bins
- Maximum refinement depth

For example, we can create an octree merging instance that splits velocity bins across the middle,
sets the root bin bounds to the bounding box of the particle velocities, and inherits the bin bounds of the parent
bin when splitting a bin. We also immediately merge our particles, setting a target particle number of 100.

```julia
# set up the merging algorithm
oc = OctreeN2Merge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel,
                   bin_bounds_compute=OctreeBinBoundsInherit, max_Nbins=6000)

# set Ntarget
Ntarget = 100

# perform merging
merge_octree_N2_based!(rng, oc, particles[1], pia, 1, 1, Ntarget)
```

We can check the number of particles after the merging procedure by looking at `pia.indexer[1,1].n_local` (the number of particles
of species 1 in cell 1), which should be equal to 100.

## Example: bringing it all together
Now we can put the variable-weight collisions and merging into a timestep loop, performing merging by when `pia.indexer[1,1].n_local > Nthreshold`,
where `Nthreshold` is the threshold number of particles. Alternatively, one can check `phys_props.np[1,1] > Nthreshold`, in case the physical
properties are computed frequently enough. The following simulation computes the BKW relaxation problem.
```julia
using Merzbild
using Random

seed = 1
Random.seed!(seed)
rng = Xoshiro(seed)

# our number density
ndens = 1e23 

# our temperature
T0 = 300.0  

# load species and interaction data
species_data = load_species_data("data/particles.toml", "Ar")
interaction_data = load_interaction_data("data/pseudo_maxwell.toml", species_data)

# number of velocity grid points in each velocity direction
nv = 20  

# some initial guess on # of particles in simulation
np_base = 40^3  

particles = [ParticleVector(np_base)]

# sample from a BKW distribution at t=0
vdf0 = (vx, vy, vz) -> bkw(vx, vy, vz, species_data[1].mass, T0, 0.0)

# returns number of particles sampled
n_sampled = sample_on_grid!(rng, vdf0, particles[1], nv, species_data[1].mass, T0, ndens,
                            0.0, 1.0, 0.0, 1.0, 0.0, 1.0;
                            v_mult=3.5, cutoff_mult=3.5, noise=0.0, v_offset=[0.0, 0.0, 0.0])

# create the ParticleIndexerArray
pia = ParticleIndexerArray(n_sampled)

# set up the merging algorithm
oc = OctreeN2Merge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel,
                   bin_bounds_compute=OctreeBinBoundsInherit, max_Nbins=6000)

# set Ntarget for merging
Ntarget = 100

# set Nthreshold for merging
Nthreshold = 120

# perform initial merge
merge_octree_N2_based!(rng, oc, particles[1], pia, 1, 1, Ntarget)

# set some reference values
sigma_ref = π * (interaction_data[1,1].vhs_d^2)
vref = sqrt(2 * k_B * T0 / species_data[1].mass)
Lref = 1.0 / (ndens * sigma_ref)
tref = Lref / vref

# set up computation of physical properties
phys_props = PhysProps(1, 1, [], Tref=T0)
compute_props!(particles, pia, species_data, phys_props)

ds = NCDataHolder("output_bkw_octree.nc", species_data, phys_props)
write_netcdf(ds, phys_props, 0)

collision_factors = create_collision_factors_array(1)  # 1-species, 0-D
collision_data = CollisionData()

# estimate average Fnum
Fnum = ndens/n_sampled

# estimate (sigma_g_w)_max
estimate_sigma_g_w_max!(collision_factors, interaction_data, species_data, [T0], Fnum)


# set number of timesteps and scaled timestep
n_t = 100
dt_scaled = 0.025

# set timestep and cell volume (0D so V=1)
Δt = dt_scaled * tref
V = 1.0

# start time loop
for ts in 1:n_t

    # collide
    ntc!(rng, collision_factors[1,1,1], collision_data, interaction_data, particles[1], pia, 1, 1, Δt, V)

    # check if we need to merge
    if pia.indexer[1,1].n_local > Nthreshold
        merge_octree_N2_based!(oc, particles[1], pia, 1, 1, Ntarget)
    end
    if ts % 10 == 0
        println("t=$ts, np=$(pia.indexer[1,1].n_local)")
    end
    
    compute_props!(particles, pia, species_data, phys_props)
    write_netcdf(ds, phys_props, ts)
end
close_netcdf(ds)
```