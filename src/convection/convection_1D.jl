@muladd begin

using StaticArrays

"""
    convect_single_particle!(rng, grid::Grid1DUniform, bc_list, particle::Particle{D}, species, Δt) where D

Convect a singe particle on a 1-D uniform grid.

# Positional arguments
* `rng`: the random number generator
* `grid`: the grid on which the convection is performed
* `bc_list`: the `Tuple` of boundary conditions (left and right wall)
* `particle`: the particle to be convected
* `species`: the index of the species being convected
* `Δt`: the convection timestep
"""
@inline function convect_single_particle!(rng, grid::Grid1DUniform, bc_list, particle::Particle{D}, species, Δt) where D
    t_rest = Δt
    @inbounds x_old = particle.x[1]
    @inbounds x_new = x_old + particle.v[1] * Δt

    @inbounds while (x_new >= grid.L) || (x_new <= 0.0)
        if (x_new >= grid.L)
            t_rest -= abs((grid.L - x_old) / particle.v[1])
            bc_id = 2
            x_old = grid.L
        else
            t_rest -= abs(x_old / particle.v[1])
            bc_id = 1
            x_old = 0.0
        end

        t_rest = apply_bc_dispatched!(rng, particle, bc_list, bc_id, grid.surface_normals[bc_id], t_rest)

        x_new = x_old + particle.v[1] * t_rest
    end

    # if a particle is too near a wall, we offset it a bit to avoid particles
    # that are exactly at a wall screwing up counters, etc.
    x_clamped = clamp(x_new, grid.min_x, grid.max_x)
    particle.x = set_x(particle.x, x_clamped)
end

"""
    convect_single_particle_periodic!(grid::Grid1DUniform, particle::Particle{D}, Δt) where D

Convect a singe particle on a 1-D uniform grid assuming a periodic grid.

# Positional arguments
* `grid`: the grid on which the convection is performed
* `particle`: the particle to be convected
* `Δt`: the convection timestep
"""
@inline function convect_single_particle_periodic!(grid::Grid1DUniform, particle::Particle{D}, Δt) where D
    @inbounds x_old = particle.x[1]
    @inbounds x_new = mod(x_old + particle.v[1] * Δt, grid.L)

    # if a particle is too near a wall, we offset it a bit to avoid particles
    # that are exactly at a wall screwing up counters, etc.
    x_clamped = clamp(x_new, grid.min_x, grid.max_x)
    particle.x = set_x(particle.x, x_clamped)
end

"""
    convect_single_particle!(rng, grid::Grid1DUniform, bc_list, particle::Particle{D}, species, surf_props::SurfProps, mass, Δt) where D

Convect a singe particle on a 1-D uniform grid, updating surface properties if it collides with a wall.

# Positional arguments
* `rng`: the random number generator
* `grid`: the grid on which the convection is performed
* `bc_list`: the `Tuple` of boundary conditions (left and right wall)
* `particle`: the particle to be convected
* `species`: the index of the species being convected
* `surf_props`: the `SurfProps` struct where the computed surface properties will be stored
* `mass`: the molecular mass of the species
* `Δt`: the convection timestep
"""
@inline function convect_single_particle!(rng, grid::Grid1DUniform, bc_list, particle::Particle{D}, species, surf_props::SurfProps, mass, Δt) where D
    t_rest = Δt
    @inbounds x_old = particle.x[1]
    @inbounds x_new = x_old + particle.v[1] * Δt

    @inbounds while (x_new >= grid.L) || (x_new <= 0.0) 
        if (x_new >= grid.L)
            t_rest -= abs((grid.L - x_old) / particle.v[1])
            bc_id = 2
            x_old = grid.L
        else
            t_rest -= abs(x_old / particle.v[1])
            bc_id = 1
            x_old = 0.0
        end

        update_surface_incident!(particle, species, surf_props, bc_id)
        t_rest = apply_bc_dispatched!(rng, particle, bc_list, bc_id, grid.surface_normals[bc_id], t_rest)
        update_surface_reflected!(particle, species, surf_props, bc_id)

        x_new = x_old + particle.v[1] * t_rest
    end

    # if a particle is too near a wall, we offset it a bit to avoid particles
    # that are exactly at a wall screwing up counters, etc.
    x_clamped = clamp(x_new, grid.min_x, grid.max_x)
    particle.x = set_x(particle.x, x_clamped)
end

"""
    convect_particles!(rng, grid::Grid1DUniform, bc_list, particles::ParticleVector{D}, pia, species, species_data, Δt) where D

Convect particles on a 1-D uniform grid.

# Positional arguments
* `rng`: the random number generator
* `grid`: the grid on which the convection is performed
* `bc_list`: the `Tuple` of boundary conditions (left and right wall)
* `particles`: the `ParticleVector` of particles to be convected
* `pia`: the `ParticleIndexerArray` instance
* `species`: the index of the species being convected
* `species_data`: the vector of `Species` data
* `Δt`: the convection timestep
"""
function convect_particles!(rng, grid::Grid1DUniform, bc_list, particles::ParticleVector{D}, pia, species, species_data, Δt) where D
    @inbounds if pia.contiguous[species]
        @inbounds n_tot = pia.n_total[species]
        @inbounds for i in 1:n_tot
            convect_single_particle!(rng, grid, bc_list, particles[i], species, Δt) 
        end
    else
        @inbounds for cell in 1:grid.n_cells
            s = pia.indexer[cell, species].start1
            e = pia.indexer[cell, species].end1
            
            for i in s:e
                convect_single_particle!(rng, grid, bc_list, particles[i], species, Δt) 
            end

            if pia.indexer[cell, species].n_group2 > 0
                s = pia.indexer[cell, species].start2
                e = pia.indexer[cell, species].end2
            
                for i in s:e
                    convect_single_particle!(rng, grid, bc_list, particles[i], species, Δt) 
                end
            end
        end
    end
end

"""
    convect_particles_periodic!(grid::Grid1DUniform, particles::ParticleVector{D}, pia, species, Δt) where D

Convect particles on a 1-D uniform grid assuming a periodic grid.

# Positional arguments
* `grid`: the grid on which the convection is performed
* `particles`: the `ParticleVector` of particles to be convected
* `pia`: the `ParticleIndexerArray` instance
* `species`: the index of the species being convected
* `Δt`: the convection timestep
"""
function convect_particles_periodic!(grid::Grid1DUniform, particles::ParticleVector{D}, pia, species, Δt) where D
    @inbounds if pia.contiguous[species]
        @inbounds n_tot = pia.n_total[species]
        @inbounds for i in 1:n_tot
            convect_single_particle_periodic!(grid, particles[i], Δt) 
        end
    else
        @inbounds for cell in 1:grid.n_cells
            s = pia.indexer[cell, species].start1
            e = pia.indexer[cell, species].end1
            
            for i in s:e
                convect_single_particle_periodic!(grid, particles[i], Δt) 
            end

            if pia.indexer[cell, species].n_group2 > 0
                s = pia.indexer[cell, species].start2
                e = pia.indexer[cell, species].end2
            
                for i in s:e
                    convect_single_particle_periodic!(grid, particles[i], Δt) 
                end
            end
        end
    end
end

"""
    convect_particles!(rng, grid::Grid1DUniform, bc_list, surf_props::SurfProps, particles::ParticleVector{D}, pia, species, species_data, Δt) where D

Convect particles on a 1-D uniform grid, computing surface properties if particles hit a surface.

# Positional arguments
* `rng`: the random number generator
* `grid`: the grid on which the convection is performed
* `bc_list`: the `Tuple` of boundary conditions (left and right wall)
* `particles`: the `ParticleVector` of particles to be convected
* `pia`: the `ParticleIndexerArray` instance
* `species`: the index of the species being convected
* `species_data`: the vector of `Species` data
* `surf_props`: the `SurfProps` struct where the computed surface properties will be stored
* `Δt`: the convection timestep
"""
function convect_particles!(rng, grid::Grid1DUniform, bc_list, particles::ParticleVector{D}, pia, species, species_data, surf_props::SurfProps, Δt) where D
    @inbounds mass = species_data[species].mass
    clear_props!(surf_props)
    @inbounds if pia.contiguous[species]
        @inbounds n_tot = pia.n_total[species]
        @inbounds for i in 1:n_tot
            convect_single_particle!(rng, grid, bc_list, particles[i], species, surf_props, mass, Δt) 
        end
    else
        @inbounds for cell in 1:grid.n_cells
            s = pia.indexer[cell, species].start1
            e = pia.indexer[cell, species].end1
            
            for i in s:e
                convect_single_particle!(rng, grid, bc_list, particles[i], species, surf_props, mass, Δt) 
            end

            if pia.indexer[cell, species].n_group2 > 0
                s = pia.indexer[cell, species].start2
                e = pia.indexer[cell, species].end2
            
                for i in s:e
                    convect_single_particle!(rng, grid, bc_list, particles[i], species, surf_props, mass, Δt) 
                end
            end
        end
    end

    surface_props_scale!(species, species_data, surf_props, Δt)
end

"""
    convect_particles_and_compute_cell!(rng, grid::Grid1DUniform, bc_list, particles::ParticleVector{D}, pia, species, species_data, Δt) where D

Convect particles on a 1-D uniform grid and write post-convection cell index to `particles.cell`.

# Positional arguments
* `rng`: the random number generator
* `grid`: the grid on which the convection is performed
* `bc_list`: the `Tuple` of boundary conditions (left and right wall)
* `particles`: the `ParticleVector` of particles to be convected
* `pia`: the `ParticleIndexerArray` instance
* `species`: the index of the species being convected
* `species_data`: the vector of `Species` data
* `Δt`: the convection timestep
"""
function convect_particles_and_compute_cell!(rng, grid::Grid1DUniform, bc_list, particles::ParticleVector{D}, pia, species, species_data, Δt) where D
    @inbounds if pia.contiguous[species]
        @inbounds n_tot = pia.n_total[species]
        @inbounds for i in 1:n_tot
            convect_single_particle!(rng, grid, bc_list, particles[i], species, Δt)
            particles.cell[i] = get_cell(grid, particles[i].x)
        end
    else
        n_cells = grid.n_cells
        @inbounds for cell in 1:n_cells
            s = pia.indexer[cell, species].start1
            e = pia.indexer[cell, species].end1
            
            for i in s:e
                convect_single_particle!(rng, grid, bc_list, particles[i], species, Δt)
                particles.cell[i] = get_cell(grid, particles[i].x)
            end

            if pia.indexer[cell, species].n_group2 > 0
                s = pia.indexer[cell, species].start2
                e = pia.indexer[cell, species].end2
            
                for i in s:e
                    convect_single_particle!(rng, grid, bc_list, particles[i], species, Δt)
                    particles.cell[i] = get_cell(grid, particles[i].x)
                end
            end
        end
    end
end

"""
    convect_particles_and_compute_cell_periodic!(grid::Grid1DUniform, particles::ParticleVector{D}, pia, species, Δt) where D

Convect particles on a 1-D uniform grid and write post-convection cell index to `particles.cell`
assuming a periodic grid.

# Positional arguments
* `grid`: the grid on which the convection is performed
* `particles`: the `ParticleVector` of particles to be convected
* `pia`: the `ParticleIndexerArray` instance
* `species`: the index of the species being convected
* `Δt`: the convection timestep
"""
function convect_particles_and_compute_cell!(grid::Grid1DUniform, particles::ParticleVector{D}, pia, species, Δt) where D
    @inbounds if pia.contiguous[species]
        @inbounds n_tot = pia.n_total[species]
        @inbounds for i in 1:n_tot
            convect_single_particle_periodic!(grid, particles[i], Δt)
            particles.cell[i] = get_cell(grid, particles[i].x)
        end
    else
        n_cells = grid.n_cells
        @inbounds for cell in 1:n_cells
            s = pia.indexer[cell, species].start1
            e = pia.indexer[cell, species].end1
            
            for i in s:e
                convect_single_particle_periodic!(grid, particles[i], Δt)
                particles.cell[i] = get_cell(grid, particles[i].x)
            end

            if pia.indexer[cell, species].n_group2 > 0
                s = pia.indexer[cell, species].start2
                e = pia.indexer[cell, species].end2
            
                for i in s:e
                    convect_single_particle_periodic!(grid, particles[i], Δt)
                    particles.cell[i] = get_cell(grid, particles[i].x)
                end
            end
        end
    end
end

"""
    convect_particles_and_compute_cell!(rng, grid::Grid1DUniform, bc_list, surf_props::SurfProps, particles::ParticleVector{D}, pia, species, species_data, Δt) where D

Convect particles on a 1-D uniform grid and write post-convection cell index to `particles.cell`, computing surface properties if particles hit a surface.

# Positional arguments
* `rng`: the random number generator
* `grid`: the grid on which the convection is performed
* `bc_list`: the `Tuple` of boundary conditions (left and right wall)
* `particles`: the `ParticleVector` of particles to be convected
* `pia`: the `ParticleIndexerArray` instance
* `species`: the index of the species being convected
* `species_data`: the vector of `Species` data
* `surf_props`: the `SurfProps` struct where the computed surface properties will be stored
* `Δt`: the convection timestep
"""
function convect_particles_and_compute_cell!(rng, grid::Grid1DUniform, bc_list, particles::ParticleVector{D}, pia, species, species_data, surf_props::SurfProps, Δt) where D
    @inbounds mass = species_data[species].mass
    clear_props!(surf_props)
    @inbounds if pia.contiguous[species]
        @inbounds n_tot = pia.n_total[species]
        @inbounds for i in 1:n_tot
            convect_single_particle!(rng, grid, bc_list, particles[i], species, surf_props, mass, Δt)
            particles.cell[i] = get_cell(grid, particles[i].x)
        end
    else
        n_cells = grid.n_cells
        for cell in 1:n_cells
            @inbounds s = pia.indexer[cell, species].start1
            @inbounds e = pia.indexer[cell, species].end1
            
            @inbounds for i in s:e
                convect_single_particle!(rng, grid, bc_list, particles[i], species, surf_props, mass, Δt)
                particles.cell[i] = get_cell(grid, particles[i].x)
            end

            @inbounds if pia.indexer[cell, species].n_group2 > 0
                @inbounds s = pia.indexer[cell, species].start2
                @inbounds e = pia.indexer[cell, species].end2
            
                @inbounds for i in s:e
                    convect_single_particle!(rng, grid, bc_list, particles[i], species, surf_props, mass, Δt) 
                    particles.cell[i] = get_cell(grid, particles[i].x)
                end
            end
        end
    end

    surface_props_scale!(species, species_data, surf_props, Δt)
end

end