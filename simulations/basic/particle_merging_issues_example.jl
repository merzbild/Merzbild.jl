include("../../src/merzbild.jl")

using ..Merzbild
using Random

function create_particles()
    vp = Vector{Particle}(undef, 64)
    v_arr = [-3.0, 0.9, 1.0, 1.1]

    ind = 1
    for v_z in v_arr
        for v_y in v_arr
            for v_x in v_arr
                vp[ind] = Particle(1.0, [v_x, v_y, v_z], [0.0, 0.0, 0.0])
                ind += 1
            end
        end
    end
    return vp
end


function run(seed)
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)
    species_list::Vector{Species} = load_species_list("data/particles.toml", "Ar")

    octree = create_merging_octree(OctreeBinMidSplit; init_bin_bounds=OctreeInitBinC)
    particles::Vector{Vector{Particle}} = [create_particles()]
    pia = create_particle_indexer_array(64)

    phys_props::PhysProps = create_props(1, 1, [], Tref=1)
    compute_props!(phys_props, pia, particles, species_list)

    println(phys_props.v[1,1])

    g_max_sq = 0.0
    m4 = 0.0
    wtot = 0.0
    for pi in 1:64
        p = particles[1][pi]
        tmp_sq = sum(p.v .* p.v)
        if tmp_sq > g_max_sq
            g_max_sq = tmp_sq
        end
        m4 += p.w * (p.v[1]^4 + p.v[2]^4 + p.v[3]^4)
        wtot += p.w
    end
    println("g_max: ", g_max_sq)
    println("T: ", phys_props.T[1,1])
    println("wtot: ", wtot)
    println("M4: ", m4/64.0)

    # merge_octree_N2_based!(1, 1, octree, particles, pia, 16)
    # println("N_post merge: ", pia.n_total[1])
    # compute_props!(phys_props, pia, particles, species_list)

    # println(phys_props.v[1,1])

    # g_max_sq = 0.0
    # for pi in 1:pia.n_total[1]
    #     p = particles[1][pi]
    #     tmp_sq = sum(p.v .* p.v)
    #     if tmp_sq > g_max_sq
    #         g_max_sq = tmp_sq
    #     end
    # end
    # println("g_max: ", g_max_sq)
    # println("T: ", phys_props.T[1,1])

    merge_octree_N2_based!(1, 1, octree, particles, pia, 2)
    println("N_post merge: ", pia.n_total[1])
    compute_props!(phys_props, pia, particles, species_list)

    println(phys_props.v[1,1])

    g_max_sq = 0.0
    m4 = 0.0
    wtot = 0.0
    for pi in 1:pia.n_total[1]
        p = particles[1][pi]
        tmp_sq = sum(p.v .* p.v)
        if tmp_sq > g_max_sq
            g_max_sq = tmp_sq
            println(p.v)
        end
        m4 += p.w * (p.v[1]^4 + p.v[2]^4 + p.v[3]^4)
        wtot += p.w
    end
    println("g_max: ", g_max_sq)
    println("T: ", phys_props.T[1,1])
    println("wtot: ", wtot)
    println("M4: ", m4/wtot)
end

run(1234)