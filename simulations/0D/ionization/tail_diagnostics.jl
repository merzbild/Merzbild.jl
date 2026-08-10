# Diagnostics for the effect of merging on the electron energy tail.
#
# The tail function is the total computational weight of particles whose kinetic energy exceeds a
# cutoff, as in simulations/0D/basic/sample_and_merge.jl. Cutoffs are placed relative to the
# ionization threshold so that one sits in the bulk (the part that sets the temperature), one just
# below the threshold, and two above it. A merge that conserves the low-order velocity moments can
# still move weight across the threshold, and it is the weight above the threshold that the
# ionization rate responds to, so a scheme can be moment-accurate and still bias the rate.
#
# Values are recorded immediately before and immediately after the electron merge of each timestep.
# Nothing else runs in between, so the difference is caused by the merge and by nothing else -- in
# particular this needs none of the controls required when merges are located after the fact from a
# drop in the particle count (merging triggers on the count crossing a threshold, which only happens
# on a step that just ionized, so merge steps found that way are selected for having ionized).

using NCDatasets

# merge_kind codes recorded per timestep
const TAIL_NO_MERGE = Int8(0)
const TAIL_MERGE_PRIMARY = Int8(1)  # octree, or the requested NNLS variant
const TAIL_MERGE_BACKUP = Int8(2)   # NNLS with fewer moments conserved
const TAIL_MERGE_OCTREE_FALLBACK = Int8(3)

"""
    TailRecorder(energies_eV, mass, n_t)

Per-timestep record of the tail functions of a species of mass `mass`, evaluated at the kinetic
energies `energies_eV` (in eV), before and after merging, over `n_t` timesteps.

The total weight is recorded alongside the tail weights: merging conserves it exactly, so the
tail functions are best read as fractions of it (the electron density grows by three orders of
magnitude over one of these runs, which makes the raw weights meaningless to average over time),
and `w_total_post / w_total_pre` doubles as a conservation check on the merge.

The spread of the particle weights is recorded at the same two points, together with the
electron-neutral NTC majorant. At matched particle count the NNLS runs create close to twice as
many computational particles per unit time as the octree ones, which is what drives their higher
merge frequency; the number of NTC candidate pairs goes as `N_e * N_n * (sigma g w)_max`, so if the
weight spread left behind by a merge is what inflates the majorant, these series show it directly.
"""
struct TailRecorder
    energies_eV::Vector{Float64}
    v_cutoffs_sq::Vector{Float64}
    tail_pre::Matrix{Float64}   # [cutoff, timestep]
    tail_post::Matrix{Float64}
    w_total_pre::Vector{Float64}
    w_total_post::Vector{Float64}
    wstats_pre::Matrix{Float64}  # [statistic, timestep]: w_max/w_min, sigma_w, sigma_ln_w
    wstats_post::Matrix{Float64}
    sigma_g_w_max::Vector{Float64}
    n_coll::Vector{Int64}
    n_eq_w_coll::Vector{Int64}
    merge_kind::Vector{Int8}

    function TailRecorder(energies_eV, mass, n_t)
        # E = m v^2 / 2, with E in eV
        v_cutoffs_sq = [2.0 * E * Merzbild.eV_J / mass for E in energies_eV]

        return new(collect(energies_eV), v_cutoffs_sq,
                   zeros(length(energies_eV), n_t), zeros(length(energies_eV), n_t),
                   zeros(n_t), zeros(n_t),
                   zeros(3, n_t), zeros(3, n_t), zeros(n_t),
                   zeros(Int64, n_t), zeros(Int64, n_t), zeros(Int8, n_t))
    end
end

# particle indices of a species may be split across two groups; the second is empty when unused
@inline group_ranges(indexer) = (indexer.start1:indexer.end1,
                                 indexer.n_group2 > 0 ? (indexer.start2:indexer.end2) : (1:0))

@inline function accumulate_particle!(out, p, v_cutoffs_sq)
    v_sq = p.v[1]^2 + p.v[2]^2 + p.v[3]^2

    @inbounds for k in eachindex(v_cutoffs_sq)
        if v_sq >= v_cutoffs_sq[k]
            out[k] += p.w
        end
    end

    return p.w
end

"""
    tail_functions!(out, pv, pia, cell, species, v_cutoffs_sq)

Accumulate the weight above each squared-speed cutoff into `out`, returning the total weight.
"""
function tail_functions!(out, pv, pia, cell, species, v_cutoffs_sq)
    fill!(out, 0.0)
    w_total = 0.0

    @inbounds for r in group_ranges(pia.indexer[cell, species]), i in r
        w_total += accumulate_particle!(out, pv[i], v_cutoffs_sq)
    end

    return w_total
end

"""
    weight_stats!(out, pv, pia, cell, species)

Spread of the computational weights of a species: `out` receives the max-to-min weight ratio, the
standard deviation of the weights, and the standard deviation of their logarithms, as in the
`ratio` function of simulations/0D/basic/sample_and_merge.jl. The log-based measure is the one to
watch, since the weights themselves span orders of magnitude.
"""
function weight_stats!(out, pv, pia, cell, species)
    ranges = group_ranges(pia.indexer[cell, species])

    w_min = Inf
    w_max = -Inf
    w_mean = 0.0
    logw_mean = 0.0
    n = 0
    n_pos = 0

    @inbounds for r in ranges, i in r
        w = pv[i].w
        w_min = min(w_min, w)
        w_max = max(w_max, w)
        w_mean += w
        n += 1

        # NNLS returns non-negative weights, so guard the logarithm rather than assuming w > 0
        if w > 0.0
            logw_mean += log(w)
            n_pos += 1
        end
    end

    if n == 0
        fill!(out, NaN)
        return nothing
    end

    w_mean /= n
    logw_mean /= max(n_pos, 1)

    w_var = 0.0
    logw_var = 0.0

    @inbounds for r in ranges, i in r
        w = pv[i].w
        w_var += (w - w_mean)^2

        if w > 0.0
            logw_var += (log(w) - logw_mean)^2
        end
    end

    out[1] = w_max / w_min
    out[2] = sqrt(w_var / n)
    out[3] = n_pos > 0 ? sqrt(logw_var / n_pos) : NaN

    return nothing
end

"""
    record_tail_pre!(rec, ts, pv, pia, cell, species, coll_factors)

Record the tail functions and weight spread as they stand before merging on timestep `ts`.

`coll_factors` is the `CollisionFactors` instance of the electron-neutral interaction, read after
the collision step of this timestep, so its counters describe the collisions just sampled: the NTC
majorant, the number of candidate pairs tested, and how many of those were between particles of
(near-)equal weight. The last is the one to divide by the second: collisions between equal-weight
particles do not split, so `n_eq_w_coll / n_coll` is the fraction of candidate pairs that created
no new particles.
"""
function record_tail_pre!(rec::TailRecorder, ts, pv, pia, cell, species, coll_factors)
    @views rec.w_total_pre[ts] = tail_functions!(rec.tail_pre[:, ts], pv, pia, cell, species,
                                                 rec.v_cutoffs_sq)
    @views weight_stats!(rec.wstats_pre[:, ts], pv, pia, cell, species)

    rec.sigma_g_w_max[ts] = coll_factors.sigma_g_w_max
    rec.n_coll[ts] = coll_factors.n_coll
    rec.n_eq_w_coll[ts] = coll_factors.n_eq_w_coll_performed

    return nothing
end

"""
    record_tail_post!(rec, ts, pv, pia, cell, species, merge_kind)

Record the tail functions as they stand after merging on timestep `ts`. If no merge happened the
pre-merge values are copied over, so that the post-merge series is a valid time series in its own
right and `merge_kind` alone selects the steps on which a merge occurred.
"""
function record_tail_post!(rec::TailRecorder, ts, pv, pia, cell, species, merge_kind)
    rec.merge_kind[ts] = merge_kind

    if merge_kind == TAIL_NO_MERGE
        @views rec.tail_post[:, ts] .= rec.tail_pre[:, ts]
        @views rec.wstats_post[:, ts] .= rec.wstats_pre[:, ts]
        rec.w_total_post[ts] = rec.w_total_pre[ts]
    else
        @views rec.w_total_post[ts] = tail_functions!(rec.tail_post[:, ts], pv, pia, cell, species,
                                                      rec.v_cutoffs_sq)
        @views weight_stats!(rec.wstats_post[:, ts], pv, pia, cell, species)
    end

    return nothing
end

"""
    write_tail_netcdf(fname, rec)

Write the recorded tail functions to a NetCDF file. Written in one go at the end of the run rather
than incrementally: the record is ~10 Float64 per timestep, which is small next to the main output
file, but note that it is held in memory for the whole run.
"""
function write_tail_netcdf(fname, rec::TailRecorder)
    ds = NCDataset(fname, "c")

    defDim(ds, "time", length(rec.merge_kind))
    defDim(ds, "cutoff", length(rec.energies_eV))

    ds.attrib["merge_kind"] = "0 = no merge, 1 = primary merge, 2 = backup NNLS, 3 = octree fallback"
    ds.attrib["COMMENT"] = "tail_pre/tail_post hold the computational weight above each cutoff " *
                           "energy immediately before/after the electron merge of each timestep; " *
                           "normalise by w_total_pre/w_total_post to get weight fractions. " *
                           "w_ratio/sigma_w/sigma_logw describe the spread of the electron weights " *
                           "at the same two points; sigma_g_w_max, n_coll and n_eq_w_coll are the " *
                           "electron-neutral NTC majorant, the number of candidate pairs tested " *
                           "and the number of those that were between (near-)equal-weight " *
                           "particles, read once per timestep after the collision step. " *
                           "equal-weight collisions do not split particles, so n_eq_w_coll/n_coll " *
                           "is the fraction of candidate pairs that created no new particles"

    defVar(ds, "cutoff_energy_eV", Float64, ("cutoff",))[:] = rec.energies_eV
    defVar(ds, "tail_pre", Float64, ("cutoff", "time"))[:, :] = rec.tail_pre
    defVar(ds, "tail_post", Float64, ("cutoff", "time"))[:, :] = rec.tail_post
    defVar(ds, "w_total_pre", Float64, ("time",))[:] = rec.w_total_pre
    defVar(ds, "w_total_post", Float64, ("time",))[:] = rec.w_total_post

    for (row, name) in enumerate(["w_ratio", "sigma_w", "sigma_logw"])
        defVar(ds, name * "_pre", Float64, ("time",))[:] = rec.wstats_pre[row, :]
        defVar(ds, name * "_post", Float64, ("time",))[:] = rec.wstats_post[row, :]
    end

    defVar(ds, "sigma_g_w_max", Float64, ("time",))[:] = rec.sigma_g_w_max
    defVar(ds, "n_coll", Int64, ("time",))[:] = rec.n_coll
    defVar(ds, "n_eq_w_coll", Int64, ("time",))[:] = rec.n_eq_w_coll
    defVar(ds, "merge_kind", Int8, ("time",))[:] = rec.merge_kind

    close(ds)

    println("Wrote tail diagnostics to $(fname)")
    return nothing
end

"""
    tail_cutoff_energies(E_ionization)

Cutoff energies, in eV, placed relative to the ionization threshold: in the bulk, just below the
threshold, just above it, and far above it. For Ar (`E_ionization` = 15.76 eV) these come out at
roughly 5.2, 11.8, 19.7 and 31.5 eV, against an electron temperature of ~6.8 eV at 400 Td.
"""
tail_cutoff_energies(E_ionization) = [0.33, 0.75, 1.25, 2.0] .* E_ionization
