using StaticArrays
using LinearAlgebra

# iterate over species, then over cells - should be cells x species then
mutable struct PhysProps
    n_cells::Int64
    n_species::Int64
    n_moments::Int64
    lpa::Vector{Int64}  # length of particle array: species
    np::Array{Int64,2}  # number of particles: cells x species
    n::Array{Float64,2}  # number density: cells x species
    v::Array{Float64,3}  # velocity: velocity component x cells x species
    T::Array{Float64,2}  # temperature: cells x species
    moment_powers::Vector{Int8}  # which moments we compute
    moments::Array{Float64,3}  # moment_id x cells x species
    Tref::Float64  # used to scale moments
end

function create_props(n_cells, n_species, moments_list; Tref=300.0)
    # println("Computing $(length(moments_list)) moments")
    return PhysProps(n_cells, n_species, length(moments_list), zeros(n_species), zeros(n_cells, n_species),
    zeros(n_cells, n_species), zeros(3, n_cells, n_species), zeros(n_cells, n_species),
    moments_list, zeros(length(moments_list), n_cells, n_species), Tref)
end

function compute_props!(phys_props, pia, particles, species_data)
    for species in 1:phys_props.n_species
        if phys_props.n_moments > 0
            moment_factor = 4 * π * (species_data[species].mass / (twopi * k_B * phys_props.Tref))^(1.5) * 0.5
            moment_vref = (species_data[species].mass / (2 * k_B * phys_props.Tref))^0.5
        end

        for cell in 1:phys_props.n_cells
            np = 0
            n = 0.0
            E = 0.0
            T = 0.0
            v = SVector{3,Float64}(0.0, 0.0, 0.0)
            phys_props.moments[:, cell, species] .= 0.0

            for i in pia.indexer[cell,species].start1:pia.indexer[cell,species].end1
                n += particles[species][i].w
                v = v + particles[species][i].v * particles[species][i].w
                np += 1
            end

            if pia.indexer[cell,species].start2 > 0
                for i in pia.indexer[cell,species].start2:pia.indexer[cell,species].end2
                    n += particles[species][i].w
                    v = v + particles[species][i].v * particles[species][i].w
                    np += 1
                end
            end

            if (n > 0.0)
                v /= n
                for i in pia.indexer[cell,species].start1:pia.indexer[cell,species].end1
                    normv = norm(particles[species][i].v - v)
                    # TODO: make normv optional, only if we include moments!
                    E += particles[species][i].w * normv^2

                    for (n_mom, m) in enumerate(phys_props.moment_powers)
                        phys_props.moments[n_mom,cell,species] += particles[species][i].w * normv^m
                    end
                end
            
                if pia.indexer[cell,species].start2 > 0
                    for i in pia.indexer[cell,species].start2:pia.indexer[cell,species].end2
                        # TODO: make normv optional, only if we include moments!
                        normv = norm(particles[species][i].v - v)
                        E += particles[species][i].w * normv^2

                        for (n_mom, m) in enumerate(phys_props.moment_powers)
                            phys_props.moments[n_mom,cell,species] += particles[species][i].w * normv^m
                        end
                    end
                end
                
                E *= 0.5 * species_data[species].mass / (n * k_B)
                T = (2.0/3.0) * E
            end

            if phys_props.n_moments > 0
                for (n_mom, m) in enumerate(phys_props.moment_powers)
                    # TODO: precompute!
                    moment_scaling = moment_factor * moment_vref^(-(3+m)) * gamma((3+m)/2)
                    phys_props.moments[n_mom,cell,species] /= (moment_scaling * n)
                end
            end
            
            phys_props.lpa[species] = length(particles[species])
            phys_props.np[cell,species] = np
            phys_props.n[cell,species] = n
            phys_props.v[:,cell,species] = v
            phys_props.T[cell,species] = T
        end
    end
end



function compute_props_sorted_without_moments!(phys_props, pia, particles, species_data)
    for species in 1:phys_props.n_species
        for cell in 1:phys_props.n_cells
            n = 0.0
            E = 0.0
            T = 0.0
            v = SVector{3,Float64}(0.0, 0.0, 0.0)

            for i in pia.indexer[cell,species].start1:pia.indexer[cell,species].end1
                n += particles[species][i].w
                v = v + particles[species][i].v * particles[species][i].w
            end

            if (n > 0.0)
                v /= n
                for i in pia.indexer[cell,species].start1:pia.indexer[cell,species].end1
                    E += particles[species][i].w * ((particles[species][i].v[1] - v[1])^2
                                                  + (particles[species][i].v[2] - v[2])^2
                                                  + (particles[species][i].v[3] - v[3])^2)
                end
                E *= 0.5 * species_data[species].mass / (n * k_B)
                T = (2.0/3.0) * E
            end
    
            phys_props.n[cell,species] = n
            phys_props.v[:,cell,species] = v
            phys_props.T[cell,species] = T
        end
    end
end