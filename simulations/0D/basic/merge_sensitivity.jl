# Measure how much a merging scheme moves its output particles in response to a tiny change in its
# input, i.e. how well-conditioned the merge map is.
#
# For each sample: draw a set of particles, draw an identical second copy, perturb the velocities of
# the copy with a small Gaussian noise, merge both sets, and compute the weighted Chamfer distance
# between the two post-merge sets. Both merges are given the same random stream, so the distance
# measures the response to the perturbation alone and not the stochasticity of the merging itself
# (with zero noise the distance is exactly zero, which is a useful check).
#
# The quantity to look at is the amplification: the post-merge Chamfer RMS divided by the RMS of the
# input displacement. A scheme that places its output particles smoothly in the input has an
# amplification of order one; a scheme whose output support jumps discontinuously as the input moves
# has a large one, even when every conserved moment is reproduced exactly.

using Merzbild
using Random

"""
    group_ranges(indexer)

Index ranges of a species' particles; they may be split across two groups, the second being empty
when unused.
"""
@inline group_ranges(indexer) = (indexer.start1:indexer.end1,
                                 indexer.n_group2 > 0 ? (indexer.start2:indexer.end2) : (1:0))

"""
    sample!(rng, pv, pia, np_base, species_data, T0, Fnum, ndens, sampling_method)

Draw a fresh set of particles. Called twice with identically seeded generators to produce two
identical sets, which avoids having to deep-copy a `ParticleVector`.
"""
function sample!(rng, pv, pia, np_base, species_data, T0, Fnum, ndens, sampling_method)
    if sampling_method == :equal_weight
        sample_particles_equal_weight!(rng, pv, pia, 1, 1, np_base, species_data[1].mass, T0, Fnum,
                                       0.0, 1.0, 0.0, 1.0, 0.0, 1.0; distribution=:Maxwellian)
    else
        sample_particles_phase_box_weighted!(rng, pv, pia, 1, 1, np_base, species_data[1].mass, T0,
                                             ndens, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0; v_mult=4)
    end
end

"""
    perturb_velocities!(rng, pv, pia, sigma)

Add an isotropic Gaussian displacement of standard deviation `sigma` to every particle velocity.
The weights are left alone, so the perturbed set carries exactly the same total weight.
"""
function perturb_velocities!(rng, pv, pia, sigma)
    @inbounds for r in group_ranges(pia.indexer[1, 1]), i in r
        p = pv[i]
        p.v = typeof(p.v)(p.v[1] + sigma * randn(rng),
                          p.v[2] + sigma * randn(rng),
                          p.v[3] + sigma * randn(rng))
    end

    return nothing
end

"""
    extract(pv, pia, vscale)

Copy the particles' velocities (scaled by `1/vscale`) and weights out into plain arrays, so the
distance computation does not have to know about the particle indexing.
"""
function extract(pv, pia, vscale)
    n = pia.indexer[1, 1].n_local
    v = zeros(3, n)
    w = zeros(n)

    k = 0
    @inbounds for r in group_ranges(pia.indexer[1, 1]), i in r
        k += 1
        p = pv[i]
        v[1, k] = p.v[1] / vscale
        v[2, k] = p.v[2] / vscale
        v[3, k] = p.v[3] / vscale
        w[k] = p.w
    end

    return v, w
end

"""
    weighted_chamfer(vA, wA, vB, wB)

Weighted Chamfer distance between two weighted point sets in velocity space:

    (1/WA) sum_a w_a min_b |v_a - v_b|^2  +  (1/WB) sum_b w_b min_a |v_b - v_a|^2

Each half is a weighted mean of the squared distance to the nearest point of the other set, so the
result has units of velocity squared and is symmetric. The sets here hold at most a few hundred
particles, so the O(NM) search is done directly.
"""
function weighted_chamfer(vA, wA, vB, wB)
    nA = length(wA)
    nB = length(wB)

    total = 0.0

    for (v1, w1, v2, n1, n2) in ((vA, wA, vB, nA, nB), (vB, wB, vA, nB, nA))
        acc = 0.0
        wsum = 0.0

        @inbounds for i in 1:n1
            dmin = Inf

            for j in 1:n2
                d = (v1[1, i] - v2[1, j])^2 + (v1[2, i] - v2[2, j])^2 + (v1[3, i] - v2[3, j])^2
                dmin = min(dmin, d)
            end

            acc += w1[i] * dmin
            wsum += w1[i]
        end

        total += acc / wsum
    end

    return total
end

"""
    run(seed, merge_method, merge_parameter, Nsamples, noise_level, io_handle;
        sampling_method=:equal_weight)

Sample `np_base` particles with total weight 1, merge them, then merge a copy whose velocities have
been perturbed by a Gaussian of standard deviation `noise_level * vref`, and report the weighted
Chamfer distance between the two merged sets over `Nsamples` independent repetitions.

Positional arguments:
* `seed`: the initial random seed
* `merge_method`: if `:octree`, octree N:2 merging is used; if `:nnls`, NNLS merging is used
* `merge_parameter`: for octree merging, the number of particles to merge down to; for NNLS merging,
the maximum total order of the mixed velocity moments to conserve
* `Nsamples`: number of independent sample-perturb-merge repetitions
* `noise_level`: standard deviation of the velocity perturbation, as a fraction of `vref`
* `io_handle`: handle of file to write output to

Keyword arguments:
* `sampling_method`: if `:equal_weight`, equal-weight particles are sampled with velocities drawn
from a Maxwell--Boltzmann distribution. If `:weighted_samples`, velocities are drawn uniformly in a
cube and the weights are set proportional to the Maxwell--Boltzmann distribution at those velocities
"""
function run(seed, merge_method, merge_parameter, Nsamples, noise_level, io_handle;
             sampling_method=:equal_weight)
    species_data::Vector{Species} = load_species_data(joinpath(MERZBILD_DATA_PATH, "particles.toml"), "Ar")

    oc = OctreeMerge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel, max_Nbins=6000)

    mim = []
    for i in 1:merge_parameter
        append!(mim, compute_multi_index_moments(i))
    end
    mnnls = NNLSMerge(mim, 40)

    ndens = 1.0
    np_base = 500
    T0::Float64 = 300.0

    vref = sqrt(2 * k_B * T0 / species_data[1].mass)
    Fnum = ndens / np_base
    sigma = noise_level * vref

    chamfer_sum = 0.0
    n_post_sum = 0.0
    Nsamples_end = 0

    for i in 1:Nsamples
        particles_a = [ParticleVector(np_base)]
        particles_b = [ParticleVector(np_base)]
        pia_a = ParticleIndexerArray(0)
        pia_b = ParticleIndexerArray(0)

        # the two sets are drawn from identically seeded generators, so they start out identical
        sample!(Xoshiro(seed + i), particles_a[1], pia_a, np_base, species_data, T0, Fnum, ndens,
                sampling_method)
        sample!(Xoshiro(seed + i), particles_b[1], pia_b, np_base, species_data, T0, Fnum, ndens,
                sampling_method)

        perturb_velocities!(Xoshiro(seed + i + Nsamples), particles_b[1], pia_b, sigma)

        # both merges get the same random stream, so only the perturbation differs between them
        flag_a = 1
        flag_b = 1

        if merge_method == :octree
            merge_octree!(Xoshiro(seed + i), oc, particles_a[1], pia_a, 1, 1, merge_parameter)
            merge_octree!(Xoshiro(seed + i), oc, particles_b[1], pia_b, 1, 1, merge_parameter)
        else
            flag_a = merge_nnls_based!(Xoshiro(seed + i), mnnls, particles_a[1], pia_a, 1, 1;
                                       vref=vref, scaling=:variance, iteration_mult=8)
            flag_b = merge_nnls_based!(Xoshiro(seed + i), mnnls, particles_b[1], pia_b, 1, 1;
                                       vref=vref, scaling=:variance, iteration_mult=8)
        end

        # a sample only counts if both merges succeeded, otherwise the two sets are not comparable
        if flag_a != -1 && flag_b != -1
            vA, wA = extract(particles_a[1], pia_a, vref)
            vB, wB = extract(particles_b[1], pia_b, vref)

            chamfer_sum += weighted_chamfer(vA, wA, vB, wB)
            n_post_sum += pia_a.indexer[1, 1].n_local
            Nsamples_end += 1
        end

        if i % 1000 == 0
            println("$i/$Nsamples")
        end
    end

    if Nsamples_end == 0
        println("all merges failed for $(merge_method), parameter $(merge_parameter)")
        return nothing
    end

    # the Chamfer distance is a mean squared displacement, so take its root to get a length; the
    # input displacement has RMS sqrt(3)*sigma since the noise is added to all three components
    chamfer_rms = sqrt(chamfer_sum / Nsamples_end)
    input_rms = sqrt(3.0) * noise_level
    n_post = n_post_sum / Nsamples_end

    for out in (stdout, io_handle)
        write(out, "method: $(merge_method), merge parameter: $(merge_parameter), " *
                   "sampling: $(sampling_method)\n")
        write(out, "  samples used: $(Nsamples_end)/$(Nsamples)\n")
        write(out, "  Npost: $(n_post)\n")
        write(out, "  noise sigma / vref: $(noise_level)\n")
        write(out, "  input RMS displacement / vref: $(input_rms)\n")
        write(out, "  post-merge Chamfer RMS / vref: $(chamfer_rms)\n")
        write(out, "  amplification: $(chamfer_rms / input_rms)\n\n")
    end

    return nothing
end

# number of independent sample-perturb-merge repetitions per configuration
n_t = 2000

# first element is the maximum total order of moments conserved when using NNLS merging,
# second element is the target number of particles for when octree merging is used. the two are
# paired so that the post-merge counts are comparable: NNLS returns exactly C(L+3,3) particles,
# i.e. 35, 56, 84, 120, 165, 220, while octree lands somewhat below its target. these are the
# parameter sets used for "Moment-preserving particle merging via non-negative least squares"
params = [[4, 36], [5, 55], [6, 85], [7, 120], [8, 164], [9, 220]]

# perturbation sizes, as a fraction of vref.
# noise_level == 0 is a check with no perturbation the two merges are handed identical inputs and
# identical random streams, so the distance has to come out exactly zero
noise_levels = [0.001, 0.01]

pref = "scratch/data"

for (sm, fname) in zip([:equal_weight, :weighted_samples], ["equalweight", "weighted"])
    io_tmp = open("$(pref)/merge_sensitivity_$(fname).log", "w")

    for paramset in params
        for noise_level in noise_levels
            run(1, :nnls, paramset[1], n_t, noise_level, io_tmp; sampling_method=sm)
            run(1, :octree, paramset[2], n_t, noise_level, io_tmp; sampling_method=sm)
        end
    end

    close(io_tmp)
end
