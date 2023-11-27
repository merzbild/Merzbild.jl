using StaticArrays
using LinearAlgebra

# iterate over species, then over cells - should be cells x species then
mutable struct PhysProps
    n::Array{Float64,2}  # cells x species
    v::Array{Float64,3}  # velocity component x cells x species
    T::Array{Float64,2}  # cells x species
    n_moments::Vector{Int8}  # which moments we compute
    moments::Array{Float64,3}  # moment id x cells x species
end

function create_props(n_cells, n_species, moments_list)
    println("Computing $(length(moments_list)) moments")
    return PhysProps(zeros(n_cells, n_species), zeros(3, n_cells, n_species), zeros(n_cells, n_species),
    moments_list, zeros(length(moments_list), n_cells, n_species))
end

function compute_props!(phys_props, particle_indexer_array, particles, n_cells, n_species, species_data)
    for species in 1:n_species
        for cell in 1:n_cells
            n = 0.0
            v = MVector{3,Float64}(0.0, 0.0, 0.0)
            E = 0.0
            T = 0.0
            phys_props.moments[:, cell, species] .= 0.0

            for i in particle_indexer_array[species,cell].start1:particle_indexer_array[species,cell].end1
                n += particles[i].w
                v += particles[i].v * particles[i].w
            end

            if particle_indexer_array[species,cell].start2 > 0
                for i in particle_indexer_array[species,cell].start2:particle_indexer_array[species,cell].end2
                    n += particles[i].w
                    v += particles[i].v * particles[i].w
                end
            end

            if (n > 0.0)
                v /= n
                for i in particle_indexer_array[species,cell].start1:particle_indexer_array[species,cell].end1
                    normv = norm(particles[i].v - v)
                    # TODO: make normv optional, only if we include moments!
                    E += particles[i].w * normv^2

                    for (n_mom, m) in enumerate(phys_props.n_moments)
                        phys_props.moments[n_mom] += particles[i].w * normv^m
                    end
                end
            
                if particle_indexer_array[species,cell].start2 > 0
                    for i in particle_indexer_array[species,cell].start2:particle_indexer_array[species,cell].end2
                        # TODO: make normv optional, only if we include moments!
                        normv = norm(particles[i].v - v)
                        E += particles[i].w * normv^2

                        for (n_mom, m) in enumerate(phys_props.n_moments)
                            phys_props.moments[n_mom,cell,species] += particles[i].w * normv^m
                        end
                    end
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