@muladd begin

"""
    convect_single_particle!(rng, grid::Grid1DUniform, boundaries::MaxwellWalls1D, particle, species, Δt)

Convect a singe particle on a 1-D uniform grid.

# Positional arguments
* `rng`: the random number generator
* `grid`: the grid on which the convection is performed
* `boundaries`: the `MaxwellWalls1D` struct describing the boundaries (it is assumed that the wall with index 1 is the left wall and
    the wall with index 2 is the right wall)
* `particles`: the particle to be convected
* `species`: the index of the species being convected
* `Δt`: the convection timestep
"""
@inline function convect_single_particle!(rng, grid::Grid1DUniform, boundaries::MaxwellWalls1D, particle, species, Δt)
    t_rest = Δt
    @inbounds x_old = particle.x[1]
    @inbounds x_new = particle.x[1] + particle.v[1] * Δt

    @inbounds while (x_new >= grid.L) || (x_new <= 0.0)
        
        if (x_new >= grid.L)
            t_rest -= abs((grid.L - x_old) / particle.v[1])
            bc_id = 2
            wall_normal = -1.0
            x_old = grid.L
        else
            t_rest -= abs(x_old / particle.v[1])
            bc_id = 1
            wall_normal = 1.0
            x_old = 0.0
        end

        reflect_particle_x!(rng, particle, boundaries.reflection_velocities_sq[bc_id, species],
                                                wall_normal,
                                                boundaries.boundaries[bc_id].v,
                                                boundaries.boundaries[bc_id].accommodation)

        x_new = x_old + particle.v[1] * t_rest
    end

    # if a particle is too near a wall, we offset it a bit to avoid particles
    # that are exactly at a wall screwing up counters, etc.
    x_new = clamp(x_new, grid.min_x, grid.max_x)
    # if x_new < grid.min_x
    #     x_new = grid.min_x
    # elseif x_new > grid.max_x
    #     x_new = grid.max_x
    # end

    @inbounds particle.x = SVector{3,Float64}(x_new, particle.x[2], particle.x[3])
end

"""
    convect_single_particle!(rng, grid::Grid1DUniform, boundaries::MaxwellWalls1D, particle, species, surf_props::SurfProps, mass, Δt)

Convect a singe particle on a 1-D uniform grid, updating surface properties if it collides with a wall.

# Positional arguments
* `rng`: the random number generator
* `grid`: the grid on which the convection is performed
* `boundaries`: the `MaxwellWalls1D` struct describing the boundaries (it is assumed that the wall with index 1 is the left wall and
    the wall with index 2 is the right wall)
* `particles`: the particle to be convected
* `species`: the index of the species being convected
* `surf_props`: the `SurfProps` struct where the computed surface properties will be stored
* `mass`: the molecular mass of the species
* `Δt`: the convection timestep
"""
@inline function convect_single_particle!(rng, grid::Grid1DUniform, boundaries::MaxwellWalls1D, particle, species, surf_props::SurfProps, mass, Δt)
    t_rest = Δt
    @inbounds x_old = particle.x[1]
    @inbounds x_new = particle.x[1] + particle.v[1] * Δt

    while (x_new >= grid.L) || (x_new <= 0.0)
        
        if (x_new >= grid.L)
            @inbounds t_rest -= abs((grid.L - x_old) / particle.v[1])
            bc_id = 2
            wall_normal = -1.0
            x_old = grid.L
        else
            @inbounds t_rest -= abs(x_old / particle.v[1])
            bc_id = 1
            wall_normal = 1.0
            x_old = 0.0
        end

        update_surface_incident!(particle, species, surf_props, bc_id)

        @inbounds reflect_particle_x!(rng, particle, boundaries.reflection_velocities_sq[bc_id, species],
                                    wall_normal,
                                    boundaries.boundaries[bc_id].v,
                                    boundaries.boundaries[bc_id].accommodation)

        update_surface_reflected!(particle, species, surf_props, bc_id)

        @inbounds x_new = x_old + particle.v[1] * t_rest
    end

    # if a particle is too near a wall, we offset it a bit to avoid particles
    # that are exactly at a wall screwing up counters, etc.
    if x_new < grid.min_x
        x_new = grid.min_x
    elseif x_new > grid.max_x
        x_new = grid.max_x
    end

    @inbounds particle.x = SVector{3,Float64}(x_new, particle.x[2], particle.x[3])
end

"""
    convect_particles!(rng, grid::Grid1DUniform, boundaries::MaxwellWalls1D, particles, pia, species, species_data, Δt)

Convect particles on a 1-D uniform grid.

# Positional arguments
* `rng`: the random number generator
* `grid`: the grid on which the convection is performed
* `boundaries`: the `MaxwellWalls1D` struct describing the boundaries (it is assumed that the wall with index 1 is the left wall and
    the wall with index 2 is the right wall)
* `particles`: the `ParticleVector` of particles to be convected
* `pia`: the `ParticleIndexerArray` instance
* `species`: the index of the species being convected
* `species_data`: the vector of `Species` data
* `Δt`: the convection timestep
"""
function convect_particles!(rng, grid::Grid1DUniform, boundaries::MaxwellWalls1D, particles, pia, species, species_data, Δt)
    # @inbounds @simd for i in 1:pia.n_total[species]
    
    @inbounds if pia.contiguous[species]
        @inbounds n_tot = pia.n_total[species]
        @inbounds @simd for i in 1:n_tot
            convect_single_particle!(rng, grid, boundaries, particles[i], species, Δt) 
        end
    else
        for cell in 1:grid.n_cells
            @inbounds s = pia.indexer[cell, species].start1
            @inbounds e = pia.indexer[cell, species].end1
            
            @inbounds @simd for i in s:e
                convect_single_particle!(rng, grid, boundaries, particles[i], species, Δt) 
            end

            @inbounds if pia.indexer[cell, species].n_group2 > 0
                @inbounds s = pia.indexer[cell, species].start2
                @inbounds e = pia.indexer[cell, species].end2
            
                @inbounds @simd for i in s:e
                    convect_single_particle!(rng, grid, boundaries, particles[i], species, Δt) 
                end
            end
        end
    end
end

"""
    convect_particles!(rng, grid::Grid1DUniform, boundaries::MaxwellWalls1D, surf_props::SurfProps, particles, pia, species, species_data, Δt)

Convect particles on a 1-D uniform grid, computing surface properties if particles hit a surface.

# Positional arguments
* `rng`: the random number generator
* `grid`: the grid on which the convection is performed
* `boundaries`: the `MaxwellWalls1D` struct describing the boundaries (it is assumed that the wall with index 1 is the left wall and
    the wall with index 2 is the right wall)
* `particles`: the `ParticleVector` of particles to be convected
* `pia`: the `ParticleIndexerArray` instance
* `species`: the index of the species being convected
* `species_data`: the vector of `Species` data
* `surf_props`: the `SurfProps` struct where the computed surface properties will be stored
* `Δt`: the convection timestep
"""
function convect_particles!(rng, grid::Grid1DUniform, boundaries::MaxwellWalls1D, particles, pia, species, species_data, surf_props::SurfProps, Δt)
    # @inbounds @simd for i in 1:pia.n_total[species]
    
    clear_props!(surf_props)
    @inbounds if pia.contiguous[species]
        @inbounds n_tot = pia.n_total[species]
        @inbounds @simd for i in 1:n_tot
            convect_single_particle!(rng, grid, boundaries, particles[i], species, surf_props, species_data[species].mass, Δt) 
        end
    else
        for cell in 1:grid.n_cells
            @inbounds s = pia.indexer[cell, species].start1
            @inbounds e = pia.indexer[cell, species].end1
            
            @inbounds @simd for i in s:e
                convect_single_particle!(rng, grid, boundaries, particles[i], species, surf_props, species_data[species].mass, Δt) 
            end

            @inbounds if pia.indexer[cell, species].n_group2 > 0
                @inbounds s = pia.indexer[cell, species].start2
                @inbounds e = pia.indexer[cell, species].end2
            
                @inbounds @simd for i in s:e
                    convect_single_particle!(rng, grid, boundaries, particles[i], species, surf_props, species_data[species].mass, Δt) 
                end
            end
        end
    end

    surface_props_scale!(species, surf_props, species_data, Δt)
end

end