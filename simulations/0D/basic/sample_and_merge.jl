# include("../../../src/Merzbild.jl")

using Merzbild
using Random

"""
    tail_function(particles, pia, v_cutoff)

Count total weight of particles that have a speed higher than a given value.

Positional arguments:
* `particles`: the `ParticleVector` instance
* `pia`: the `ParticleIndexerArray` structure
* `v_cutoff`: the threshold value of the speed

Returns:
* The computational weight of particles that have a speed higher than `v_cutoff`.
"""
function tail_function(particles, pia, v_cutoff)
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
    ratio(particles, pia)

Compute statistics of particle weights: ratio of largest-to-smallest weight, standard deviation,
log-standard deviation.


Positional arguments:
* `particles`: the `ParticleVector` of particles
* `pia`: the `ParticleIndexerArray` instance
"""
function ratio(particles, pia)
    w_min = 1.0
    w_max = -1.0

    w_log_mean = 0.0
    w_mean = 0.0
    w_std = 0.0
    w_log_std = 0.0
    @inbounds for pin in 1:pia.indexer[1,1].n_local
        w_min = min(particles[pin].w, w_min)
        w_max = max(particles[pin].w, w_max)
        w_mean += particles[pin].w
        w_log_mean += log(particles[pin].w)
    end

    w_mean /= pia.indexer[1,1].n_local
    w_log_mean /= pia.indexer[1,1].n_local
    @inbounds for pin in 1:pia.indexer[1,1].n_local
        w_std += (particles[pin].w - w_mean)^2
        w_log_std += (log(particles[pin].w) - w_log_mean)^2
    end

    return w_max / w_min, sqrt(w_std / pia.indexer[1,1].n_local), sqrt(w_log_std / pia.indexer[1,1].n_local)
end


"""
    run(seed, merge_method, merge_parameter, Nsamples, io_handle; sampling_method=:equal_weight)

Sample 500 particles with total weight 1, merge them using either octree merging or NNLS merging,
compute how the tail functionals and weight distribution are affected by the merging
by repeating the sample-and-merge processed over `Nsamples` independent runs.
The function prints the various statistics after the computation is finished and writes them to a file as well.

Positional arguments:
* `seed`: the initial random seed
* `merge_method`: if `:octree`, then octree N:2 merging is used; if `:nnls`, NNLS merging is used
* `merge_parameter`: if octree merging is used, this is the number of particles to merge down to;
if NNLS merging is used, this is the maximum total order of the mixed velocity moments to conserve
* `Nsamples`: number of independent sample-and-merge simulations to run
* `io_handle`: handle of file to write output to

Keyword arguments:
* `sampling_method`: if `:equal_weight`, equal_weight particles are sampled by sampling their velocities
from a Maxwell--Boltzmann distribution. If `:weighted_samples`, the particle velocities are sampled
from a uniform distribution in a cube, and particles' weights are computed as being proportional
to the value of the Maxwell--Boltzmann distribution at their velocities
"""
function run(seed, merge_method, merge_parameter, Nsamples, io_handle; sampling_method=:equal_weight)
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    species_data::Vector{Species} = load_species_data("data/particles.toml", "Ar")
    oc = OctreeN2Merge(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinMinMaxVel, max_Nbins=6000)

    mim = []
    n_moms = merge_parameter
    for i in 1:n_moms
        append!(mim, compute_multi_index_moments(i))
    end
    mnnls = NNLSMerge(mim, 40)
    
    ndens = 1.0
    np_base = 500

    T0::Float64 = 300.0

    vref = sqrt(2 * k_B * T0 / species_data[1].mass)
    Fnum = ndens / np_base

    r_w_pre = 0.0
    std_w_pre = 0.0
    std_logw_pre = 0.0

    r_w = 0.0
    std_w = 0.0
    std_logw = 0.0

    tail_vel_1 = 500.0
    tail_vel_2 = 750.0

    w_tail_pre_1 = 0.0
    w_tail_pre_2 = 0.0

    w_tail_post_1 = 0.0
    w_tail_post_2 = 0.0

    n_post = 0.0

    Nsamples_end = 0
    for i in 1:Nsamples
        particles::Vector{ParticleVector} = [ParticleVector(np_base)]

        pia = ParticleIndexerArray(0)

        if sampling_method == :equal_weight
            sample_particles_equal_weight!(rng, particles[1], pia, 1, 1, np_base, species_data[1].mass, T0, Fnum,
                                           0.0, 1.0, 0.0, 1.0, 0.0, 1.0; distribution=:Maxwellian)
        else
            sample_particles_phase_box_weighted!(rng, particles[1], pia, 1, 1, np_base, species_data[1].mass, T0, ndens,
                                                 0.0, 1.0, 0.0, 1.0, 0.0, 1.0; v_mult=4)
        end

        w_tail_pre_1 += tail_function(particles[1], pia, tail_vel_1)
        w_tail_pre_2 += tail_function(particles[1], pia, tail_vel_2)

        rwo = ratio(particles[1], pia)
        r_w_pre += rwo[1]
        std_w_pre += rwo[2]
        std_logw_pre += rwo[3]

        flag = 1
        if merge_method == :octree
            merge_octree_N2_based!(rng, oc, particles[1], pia, 1, 1, merge_parameter)
            Nsamples_end += 1
        else
            flag = merge_nnls_based!(rng, mnnls, particles[1], pia, 1, 1;
                                     vref=vref, scaling=:variance, centered_at_mean=false, v_multipliers=[], iteration_mult=8)

            if flag != -1
                Nsamples_end += 1
            end
        end

        if flag != -1
            w_tail_post_1 += tail_function(particles[1], pia, tail_vel_1)
            w_tail_post_2 += tail_function(particles[1], pia, tail_vel_2)

            rwo = ratio(particles[1], pia)
            r_w += rwo[1]
            std_w += rwo[2]
            std_logw += rwo[3]
            n_post += pia.indexer[1,1].n_local
        end

        if i%10000 == 0
            println("$i/$Nsamples")
        end
    end

    r_w_pre /= Nsamples_end
    std_w_pre /= Nsamples_end
    std_logw_pre /= Nsamples_end

    r_w /= Nsamples_end
    std_w /= Nsamples_end
    std_logw /= Nsamples_end

    w_tail_pre_1 /= Nsamples_end
    w_tail_pre_2 /= Nsamples_end
    w_tail_post_1 /= Nsamples_end
    w_tail_post_2 /= Nsamples_end

    n_post /= Nsamples_end

    println("Npost = $n_post")
    println("pre-merge weight ratio: ", r_w_pre)
    println("pre-merge weight std: ", std_w_pre)
    println("pre-merge weight log std: ", std_logw_pre)
    println("weight ratio: ", r_w)
    println("weight std: ", std_w)
    println("weight log std: ", std_logw)
    println("f_tail($(tail_vel_1)): $(w_tail_pre_1) -> $(w_tail_post_1)")
    println("f_tail($(tail_vel_2)): $(w_tail_pre_2) -> $(w_tail_post_2)")
    println()

    write(io_handle, "Npost = $n_post\n")
    write(io_handle, "pre-merge weight ratio: $(r_w_pre)\n")
    write(io_handle, "pre-merge weight std: $(std_w_pre)\n")
    write(io_handle, "pre-merge weight log std: $(std_logw_pre)\n")
    write(io_handle, "weight ratio: $(r_w)\n")
    write(io_handle, "weight std: $(std_w)\n")
    write(io_handle, "weight log std: $(std_logw)\n")
    write(io_handle, "f_tail($(tail_vel_1)): $(w_tail_pre_1) -> $(w_tail_post_1)\n")
    write(io_handle, "f_tail($(tail_vel_2)): $(w_tail_pre_2) -> $(w_tail_post_2)\n\n")
end

# n_t is the number of samples for each run
n_t = 10000

# first element in a list in params is the maximum total order of moments conserved when using NNLS merging
# second element is the target number of particles for when octree merging is used
params = [[4, 36]]

#  # uncomment lines below to run over the parameter sets used for "Moment-preserving particle merging via non-negative least squares"
# params = [[4, 36], [5, 55], [6, 85], [7, 120], [8, 164], [9, 220]]
# n_t = 100000  # number of samples for each run

# pref = "scratch/data"

# for (sm, fname) in zip([:equal_weight, :weighted_samples], ["equalweight", "weighted"])
#     io_nnls = open("$(pref)/nnls_$(fname).log", "w")
#     io_octree = open("$(pref)/octree_$(fname).log", "w")

#     for paramset in params
#         println("NNLS merging: $paramset; sampling method: $(sm)")
#         run(1, :nnls, paramset[1], n_t, io_nnls; sampling_method=sm)

#         println("Octree merging: $paramset; sampling method: $(sm)")
#         run(1, :octree, paramset[2], n_t, io_octree; sampling_method=sm)
#     end

#     close(io_nnls)
#     close(io_octree)
# end
