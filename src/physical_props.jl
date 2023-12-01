using StaticArrays
using LinearAlgebra

# iterate over species, then over cells - should be cells x species then
mutable struct PhysProps
    n_cells::Int64
    n_species::Int64
    n_moments::Int64
    n::Array{Float64,2}  # cells x species
    v::Array{Float64,3}  # cells x species x velocity component
    T::Array{Float64,2}  # cells x species
    moment_powers::Vector{Int8}  # which moments we compute
    moments::Array{Float64,3}  # cells x species x moment_id
    Tref::Float64  # used to scale moments
end

function create_props(n_cells, n_species, moments_list; Tref=300.0)
    println("Computing $(length(moments_list)) moments")
    return PhysProps(n_cells, n_species, length(moments_list), zeros(n_cells, n_species), zeros(n_cells, n_species, 3), zeros(n_cells, n_species),
    moments_list, zeros(n_cells, n_species, length(moments_list)), Tref)
end

function compute_props!(phys_props, particle_indexer_array, particles, species_data)
    for species in 1:phys_props.n_species
        if phys_props.n_moments > 0
            moment_factor = 4 * π * (species_data[species].mass / (2 * π * k_B * phys_props.Tref))^(1.5) * 0.5
            moment_vref = (species_data[species].mass / (2 * k_B * phys_props.Tref))^0.5
        end

        for cell in 1:phys_props.n_cells
            n = 0.0
            v = MVector{3,Float64}(0.0, 0.0, 0.0)
            E = 0.0
            T = 0.0
            phys_props.moments[cell, species, :] .= 0.0

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

                    for (n_mom, m) in enumerate(phys_props.moment_powers)
                        phys_props.moments[cell,species,n_mom] += particles[i].w * normv^m
                    end
                end
            
                if particle_indexer_array[species,cell].start2 > 0
                    for i in particle_indexer_array[species,cell].start2:particle_indexer_array[species,cell].end2
                        # TODO: make normv optional, only if we include moments!
                        normv = norm(particles[i].v - v)
                        E += particles[i].w * normv^2

                        for (n_mom, m) in enumerate(phys_props.moment_powers)
                            phys_props.moments[cell,species,n_mom] += particles[i].w * normv^m
                        end
                    end
                end
                
                E *= 0.5 * species_data[species].mass / (n * k_B)
                T = (2.0/3.0) * E
            end

            if phys_props.n_moments > 0
                for (n_mom, m) in enumerate(phys_props.moment_powers)
                    moment_scaling = moment_factor * moment_vref^(-(3+m)) * gamma((3+m)/2)
                    phys_props.moments[cell,species,n_mom] /= (moment_scaling * n)
                end
            end
    

            phys_props.n[cell,species] = n
            phys_props.v[cell,species,:] = v
            phys_props.T[cell,species] = T
        end
    end
end