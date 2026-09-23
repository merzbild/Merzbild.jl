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

# Returns
Whether the particle left the domain through an outflow boundary and is to be deleted.
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

        if bc_is_outflow(bc_list, bc_id)
            return true
        end

        t_rest = apply_bc_dispatched!(rng, particle, bc_list, bc_id, grid.surface_normals[bc_id], t_rest)

        x_new = x_old + particle.v[1] * t_rest
    end

    # if a particle is too near a wall, we offset it a bit to avoid particles
    # that are exactly at a wall screwing up counters, etc.
    x_clamped = clamp(x_new, grid.min_x, grid.max_x)
    particle.x = set_x(particle.x, x_clamped)
    return false
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

# Returns
Whether the particle left the domain through an outflow boundary and is to be deleted.
A particle leaving through an outflow boundary is counted as incident on the surface,
but produces no reflected contribution.
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

        if bc_is_outflow(bc_list, bc_id)
            return true
        end

        t_rest = apply_bc_dispatched!(rng, particle, bc_list, bc_id, grid.surface_normals[bc_id], t_rest)
        update_surface_reflected!(particle, species, surf_props, bc_id)

        x_new = x_old + particle.v[1] * t_rest
    end

    # if a particle is too near a wall, we offset it a bit to avoid particles
    # that are exactly at a wall screwing up counters, etc.
    x_clamped = clamp(x_new, grid.min_x, grid.max_x)
    particle.x = set_x(particle.x, x_clamped)
    return false
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

# Returns
The number of particles deleted at outflow boundaries (0 if `bc_list` contains no [`OutflowBC`](@ref)).
"""
function convect_particles!(rng, grid::Grid1DUniform, bc_list, particles::ParticleVector{D}, pia, species, species_data, Δt) where D
    if has_outflow_bc(bc_list)
        n_deleted = 0

        @inbounds for cell in 1:grid.n_cells
            indexer = pia.indexer[cell, species]

            # iterating in reverse: delete_particle! swaps the deleted particle with the one at the
            # end of the group, which going backwards has already been convected, so nothing is
            # skipped or convected twice
            # no need for delete_particle_and_cell! since cell index is not tracked in this routine
            if indexer.n_group2 > 0
                for i in indexer.end2:-1:indexer.start2
                    if convect_single_particle!(rng, grid, bc_list, particles[i], species, Δt)
                        delete_particle!(particles, pia, cell, species, i)
                        n_deleted += 1
                    end
                end
            end

            for i in indexer.end1:-1:indexer.start1
                if convect_single_particle!(rng, grid, bc_list, particles[i], species, Δt)
                    delete_particle!(particles, pia, cell, species, i)
                    n_deleted += 1
                end
            end
        end

        @inbounds if n_deleted > 0
            pia.contiguous[species] = false
        end

        return n_deleted
    end

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
    return 0
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

# Returns
0 (the number of particles deleted at outflow boundaries, to be consistent with other convection routines).
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
    return 0
end

"""
    convect_particles!(rng, grid::Grid1DUniform, bc_list, particles::ParticleVector{D}, pia, species, species_data, surf_props::SurfProps, Δt) where D

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

# Returns
The number of particles deleted at outflow boundaries (0 if `bc_list` contains no [`OutflowBC`](@ref)).
"""
function convect_particles!(rng, grid::Grid1DUniform, bc_list, particles::ParticleVector{D}, pia, species, species_data, surf_props::SurfProps, Δt) where D
    @inbounds mass = species_data[species].mass
    clear_props!(surf_props)

    if has_outflow_bc(bc_list)
        n_deleted = 0

        @inbounds for cell in 1:grid.n_cells
            indexer = pia.indexer[cell, species]

            # iterating in reverse: delete_particle! swaps the deleted particle with the one at the
            # end of the group, which going backwards has already been convected, so nothing is
            # skipped or convected twice
            # no need for delete_particle_and_cell! since cell index is not tracked in this routine
            if indexer.n_group2 > 0
                for i in indexer.end2:-1:indexer.start2
                    if convect_single_particle!(rng, grid, bc_list, particles[i], species, surf_props, mass, Δt)
                        delete_particle!(particles, pia, cell, species, i)
                        n_deleted += 1
                    end
                end
            end

            for i in indexer.end1:-1:indexer.start1
                if convect_single_particle!(rng, grid, bc_list, particles[i], species, surf_props, mass, Δt)
                    delete_particle!(particles, pia, cell, species, i)
                    n_deleted += 1
                end
            end
        end

        @inbounds if n_deleted > 0
            pia.contiguous[species] = false
        end

        surface_props_scale!(species, species_data, surf_props, Δt)
        return n_deleted
    end

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
    return 0
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

# Returns
The number of particles deleted at outflow boundaries (0 if `bc_list` contains no [`OutflowBC`](@ref)).
"""
function convect_particles_and_compute_cell!(rng, grid::Grid1DUniform, bc_list, particles::ParticleVector{D}, pia, species, species_data, Δt) where D
    if has_outflow_bc(bc_list)
        n_deleted = 0

        @inbounds for cell in 1:grid.n_cells
            indexer = pia.indexer[cell, species]

            if indexer.n_group2 > 0
                for i in indexer.end2:-1:indexer.start2
                    p = particles[i]
                    if convect_single_particle!(rng, grid, bc_list, p, species, Δt)
                        # the cell index is swapped along with the particle index, so the particle
                        # moved into slot i keeps the cell index computed when it was convected
                        delete_particle_and_cell!(particles, pia, cell, species, i)
                        n_deleted += 1
                    else
                        particles.cell[i] = get_cell(grid, p.x)
                    end
                end
            end

            for i in indexer.end1:-1:indexer.start1
                p = particles[i]
                if convect_single_particle!(rng, grid, bc_list, p, species, Δt)
                    delete_particle_and_cell!(particles, pia, cell, species, i)
                    n_deleted += 1
                else
                    particles.cell[i] = get_cell(grid, p.x)
                end
            end
        end

        @inbounds if n_deleted > 0
            pia.contiguous[species] = false
        end

        return n_deleted
    end

    @inbounds if pia.contiguous[species]
        @inbounds n_tot = pia.n_total[species]
        @inbounds for i in 1:n_tot
            p = particles[i]
            convect_single_particle!(rng, grid, bc_list, p, species, Δt)
            particles.cell[i] = get_cell(grid, p.x)
        end
    else
        n_cells = grid.n_cells
        @inbounds for cell in 1:n_cells
            s = pia.indexer[cell, species].start1
            e = pia.indexer[cell, species].end1
            
            for i in s:e
                p = particles[i]
                convect_single_particle!(rng, grid, bc_list, p, species, Δt)
                particles.cell[i] = get_cell(grid, p.x)
            end

            if pia.indexer[cell, species].n_group2 > 0
                s = pia.indexer[cell, species].start2
                e = pia.indexer[cell, species].end2
            
                for i in s:e
                    p = particles[i]
                    convect_single_particle!(rng, grid, bc_list, p, species, Δt)
                    particles.cell[i] = get_cell(grid, p.x)
                end
            end
        end
    end
    return 0
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

# Returns
0 (the number of particles deleted at outflow boundaries, to be consistent with other convection routines).
"""
function convect_particles_and_compute_cell_periodic!(grid::Grid1DUniform, particles::ParticleVector{D}, pia, species, Δt) where D
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
    return 0
end

"""
    convect_particles_and_compute_cell!(rng, grid::Grid1DUniform, bc_list, particles::ParticleVector{D}, pia, species, species_data, surf_props::SurfProps, Δt) where D

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

# Returns
The number of particles deleted at outflow boundaries (0 if `bc_list` contains no [`OutflowBC`](@ref)).
"""
function convect_particles_and_compute_cell!(rng, grid::Grid1DUniform, bc_list, particles::ParticleVector{D}, pia, species, species_data, surf_props::SurfProps, Δt) where D
    @inbounds mass = species_data[species].mass
    clear_props!(surf_props)

    if has_outflow_bc(bc_list)
        n_deleted = 0

        @inbounds for cell in 1:grid.n_cells
            indexer = pia.indexer[cell, species]

            if indexer.n_group2 > 0
                for i in indexer.end2:-1:indexer.start2
                    p = particles[i]
                    if convect_single_particle!(rng, grid, bc_list, p, species, surf_props, mass, Δt)
                        # the cell index is swapped along with the particle index, so the particle
                        # moved into slot i keeps the cell index computed when it was convected
                        delete_particle_and_cell!(particles, pia, cell, species, i)
                        n_deleted += 1
                    else
                        particles.cell[i] = get_cell(grid, p.x)
                    end
                end
            end

            for i in indexer.end1:-1:indexer.start1
                p = particles[i]
                if convect_single_particle!(rng, grid, bc_list, p, species, surf_props, mass, Δt)
                    delete_particle_and_cell!(particles, pia, cell, species, i)
                    n_deleted += 1
                else
                    particles.cell[i] = get_cell(grid, p.x)
                end
            end
        end

        @inbounds if n_deleted > 0
            pia.contiguous[species] = false
        end

        surface_props_scale!(species, species_data, surf_props, Δt)
        return n_deleted
    end

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
    return 0
end

end