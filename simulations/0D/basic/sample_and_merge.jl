include("../../../src/Merzbild.jl")

using ..Merzbild
using Random

"""
    tail_functional(particles, pia, v_cutoff)

Count total weight of particles that have a speed higher than a given value.

Positional arguments:
* `particles`: the `ParticleVector` instance
* `pia`: the `ParticleIndexerArray` structure
* `v_cutoff`: the threshold value of the speed

Returns:
* The computational weight of particles that have a speed higher than `v_cutoff`.
"""
function tail_functional(particles, pia, v_cutoff)
    w = 0.0
    v_c_sq = v_cutoff^2
    @inbounds for pin in 1:pia.indexer[1,1].n_local
        if particles[pin].v[1]^2 + particles[pin].v[2]^2 + particles[pin].v[3]^2 >= v_c_sq
            w += particles[pin].w
        end
    end
    return w
end

"""
    run(seed, merge_method, merge_parameter, Nsamples; sampling_method=:equal_weight)

Sample 500 particles with total weight 1, merge them using either octree merging or NNLS merging,
compute how the tail functionals and weight distribution are affected by the merging
by repeating the sample-and-merge processed over `Nsamples` independent runs.
The function prints the various statistics after the computation is finished.

Positional arguments:
* `seed`: the initial random seed
* `merge_method`: if `:octree`, then octree N:2 merging is used; if `:nnls`, NNLS merging is used
* `merge_parameter`: if octree merging is used, this is the number of particles to merge down to;
if NNLS merging is used, this is the maximum total order of the mixed velocity moments to conserve
* `Nsamples`: number of independent sample-and-merge simulations to run

Keyword arguments:
* `sampling_method`: if `:equal_weight`, equal_weight particles are sampled by sampling their velocities
from a Maxwell--Boltzmann distribution. If `:weighted_samples`, the particle velocities are sampled
from a uniform distribution in a cube, and particles' weights are computed as being proportional
to the value of the Maxwell--Boltzmann distribution at their velocities.
"""
function run(seed, merge_parameter, Nsamples; sampling_method=:equal_weight)

end