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
Whether the particle left the domain and is to be deleted.
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

        if bc_deletes_particles(bc_list, bc_id)
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
Whether the particle left the domain and is to be deleted.
A particle leaving the domain is counted as incident on the surface,
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

        if bc_deletes_particles(bc_list, bc_id)
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
The number of particles that left the domain and were deleted (0 if `bc_list` contains no [`Merzbild.ParticleDeletingBC`](@ref)).
"""
function convect_particles!(rng, grid::Grid1DUniform, bc_list, particles::ParticleVector{D}, pia, species, species_data, Δt) where D
    if has_particle_deleting_bc(bc_list)
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
0 (the number of deleted particles, to be consistent with other convection routines).
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
The number of particles that left the domain and were deleted (0 if `bc_list` contains no [`Merzbild.ParticleDeletingBC`](@ref)).
"""
function convect_particles!(rng, grid::Grid1DUniform, bc_list, particles::ParticleVector{D}, pia, species, species_data, surf_props::SurfProps, Δt) where D
    @inbounds mass = species_data[species].mass
    clear_props!(surf_props)

    if has_particle_deleting_bc(bc_list)
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
The number of particles that left the domain and were deleted (0 if `bc_list` contains no [`Merzbild.ParticleDeletingBC`](@ref)).
"""
function convect_particles_and_compute_cell!(rng, grid::Grid1DUniform, bc_list, particles::ParticleVector{D}, pia, species, species_data, Δt) where D
    if has_particle_deleting_bc(bc_list)
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
0 (the number of deleted particles, to be consistent with other convection routines).
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
The number of particles that left the domain and were deleted (0 if `bc_list` contains no [`Merzbild.ParticleDeletingBC`](@ref)).
"""
function convect_particles_and_compute_cell!(rng, grid::Grid1DUniform, bc_list, particles::ParticleVector{D}, pia, species, species_data, surf_props::SurfProps, Δt) where D
    @inbounds mass = species_data[species].mass
    clear_props!(surf_props)

    if has_particle_deleting_bc(bc_list)
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

"""
    convect_emitted_particle!(rng, grid::Grid1DUniform, bc_list, particle::Particle{D}, species, surf_props, mass, Fnum, surf_factor, bc_id, Δt) where D

Convect a particle emitted from boundary `bc_id` of a 1-D uniform grid for its residual flight time.
If `surf_props` is a `SurfProps` instance, the emitted particle is recorded as reflected from boundary `bc_id`,
and any impacts during the flight are recorded as well, with contributions already scaled
as by [`Merzbild.surface_props_scale!`](@ref). This scaling uses the area of the boundary `bc_id` also for impacts on the
other boundary, which is exact for equal areas (as is the case for `SurfProps(pia, grid::Grid1DUniform)`).
If `surf_props` is `nothing`, no surface properties are recorded.

# Positional arguments
* `rng`: the random number generator
* `grid`: the grid on which the convection is performed
* `bc_list`: the `Tuple` of boundary conditions (left and right wall)
* `particle`: the emitted particle
* `species`: the index of the species of the particle
* `surf_props`: the `SurfProps` struct where the computed surface properties will be stored, or `nothing`
* `mass`: the molecular mass of the species
* `Fnum`: the computational weight of the particle
* `surf_factor`: the scaling factor of the surface properties, `mass / (Δt_full * area)`,
    where `Δt_full` is the full timestep
* `bc_id`: the index of the boundary the particle is emitted from
* `Δt`: the residual flight time

# Returns
Whether the particle left the domain and is to be discarded.
"""
@inline function convect_emitted_particle!(rng, grid::Grid1DUniform, bc_list, particle::Particle{D}, species, surf_props::Nothing,
                                           mass, Fnum, surf_factor, bc_id, Δt) where D
    return convect_single_particle!(rng, grid, bc_list, particle, species, Δt)
end

@inline function convect_emitted_particle!(rng, grid::Grid1DUniform, bc_list, particle::Particle{D}, species, surf_props::SurfProps,
                                           mass, Fnum, surf_factor, bc_id, Δt) where D
    # the recorded quantities (apart from the impact count) are linear in the weight,
    # and the boundary conditions CURRENTLY do not depend on it
    # TODO: remove this weight hack in future versions
    particle.w = Fnum * surf_factor
    update_surface_reflected!(particle, species, surf_props, bc_id)
    left_domain = convect_single_particle!(rng, grid, bc_list, particle, species, surf_props, mass, Δt)
    particle.w = Fnum
    return left_domain
end

"""
    emit_particles_bc!(rng, grid::Grid1DUniform, bc_list, bc::AbstractBC, particles::ParticleVector{D}, pia, species, surf_props, mass, Fnum, bc_id, Δt) where D

Emit particles from boundary `bc_id` of a 1-D uniform grid over a timestep and convect them. Boundary conditions
other than [`InflowBC1D`](@ref) do not emit particles.

For an `InflowBC1D` with a uniform particle weight, the entries (particle entering the domain)
are a Poisson process with rate `bc.ndens_flux / Fnum` (per unit area),
and each particle is convected for the time remaining in the timestep after its entry; emitted particles leaving
the domain during that time are discarded. The particles are appended to the 2nd group of the cell
given by [`Merzbild.get_group2_append_cell`](@ref), and their actual cell indices are written to `particles.cell`,
so that the 4-argument version of [`sort_particles!`](@ref) can be used (if the convection routine
for the non-inflow particles also writes post-convection cell indices).

# Positional arguments
* `rng`: the random number generator
* `grid`: the grid into which the particles are emitted
* `bc_list`: the `Tuple` of boundary conditions (left and right wall)
* `bc`: the boundary condition of boundary `bc_id`
* `particles`: the `ParticleVector` of particles
* `pia`: the `ParticleIndexerArray` instance
* `species`: the index of the species being emitted
* `surf_props`: the `SurfProps` struct where the computed surface properties will be stored, or `nothing`
* `mass`: the molecular mass of the species
* `Fnum`: the computational weight of the emitted particles
* `bc_id`: the index of the boundary
* `Δt`: the timestep

# Returns
The number of particles added.
"""
@inline function emit_particles_bc!(rng, grid::Grid1DUniform, bc_list, bc::AbstractBC, particles::ParticleVector{D}, pia,
                                    species, surf_props, mass, Fnum, bc_id, Δt) where D
    return 0
end

function emit_particles_bc!(rng, grid::Grid1DUniform, bc_list, bc::InflowBC1D, particles::ParticleVector{D}, pia,
                            species, surf_props, mass, Fnum, bc_id, Δt) where D
    # @inbounds if bc.normal[1] != grid.surface_normals[bc_id][1]
    #     throw(ArgumentError("InflowBC1D with inward normal " * string(bc.normal) * " is used at boundary " *
    #                         string(bc_id) * " with surface normal " * string(grid.surface_normals[bc_id])))
    # end

    if bc.ndens_flux == 0.0
        return 0
    end

    if bc_id == 1
        x0 = set_x(zero(SVector{D,Float64}), grid.min_x)
        cell = get_group2_append_cell(pia, 1, species)
    else
        x0 = set_x(zero(SVector{D,Float64}), grid.max_x)
        cell = get_group2_append_cell(pia, grid.n_cells, species)
    end

    @inbounds surf_factor = surf_props === nothing ? 0.0 : mass / Δt * surf_props.inv_areas[bc_id]
    inv_rate = Fnum / bc.ndens_flux

    n_emitted = 0
    t = randexp(rng) * inv_rate

    while t < Δt
        @inbounds if length(particles) <= pia.index_last[species]
            resize!(particles, length(particles) + DELTA_PARTICLES)
        end

        # the slot that update_buffer_index_new_particle! takes from the buffer
        @inbounds particle = particles.particles[particles.buffer[particles.nbuffer]]
        particle.w = Fnum
        particle.v = sample_inflow_velocity(rng, bc)
        particle.x = x0

        if convect_emitted_particle!(rng, grid, bc_list, particle, species, surf_props, mass, Fnum, surf_factor, bc_id, Δt - t)
            particle.w = 0.0
        else
            update_buffer_index_new_particle!(particles, pia, cell, species)
            @inbounds particles.cell[pia.index_last[species]] = get_cell(grid, particle.x)
            n_emitted += 1
        end

        t += randexp(rng) * inv_rate
    end

    return n_emitted
end

"""
    emit_particles!(rng, grid::Grid1DUniform, bc_list, particles::ParticleVector{D}, pia, species, species_data, Fnum, Δt) where D

Emit particles across the [`InflowBC1D`](@ref) boundaries of a 1-D uniform grid over a timestep.
The number of particles entering across a boundary is Poisson-distributed with mean `Γ Δt / Fnum`,
where `Γ` is the number flux of the boundary condition (per unit area), and each particle is convected for the time
remaining in the timestep after its entry. The function is therefore to be called after the particles in the domain
have been convected, and before they are sorted.
Particles leaving the domain during their residual flight (either
through an [`OutflowBC`](@ref) or through an [`InflowBC1D`](@ref)) are discarded.

The new particles are appended to the 2nd group of a single cell (see [`Merzbild.get_group2_append_cell`](@ref)),
not necessarily the one they are located in, and their cell indices are written to `particles.cell`,
so the particles have to be sorted before any per-cell operation.

# Positional arguments
* `rng`: the random number generator
* `grid`: the grid into which the particles are emitted
* `bc_list`: the `Tuple` of boundary conditions (left and right wall)
* `particles`: the `ParticleVector` of particles
* `pia`: the `ParticleIndexerArray` instance
* `species`: the index of the species being emitted
* `species_data`: the vector of `Species` data
* `Fnum`: the computational weight of the emitted particles
* `Δt`: the timestep

# Returns
The number of particles added.
"""
function emit_particles!(rng, grid::Grid1DUniform, bc_list, particles::ParticleVector{D}, pia, species, species_data,
                         Fnum, Δt) where D
    return emit_particles!(rng, grid, bc_list, particles, pia, species, species_data, nothing, Fnum, Δt, 1:grid.n_cells)
end

"""
    emit_particles!(rng, grid::Grid1DUniform, bc_list, particles::ParticleVector{D}, pia, species, species_data, Fnum, Δt, cell_chunk) where D

Emit uniformly weighted particles across the [`InflowBC1D`](@ref) boundaries of a 1-D uniform grid over a timestep,
for a chunk of cells: particles are emitted across the left (right) boundary only if the first (last) cell
of the grid is the first (last) element of `cell_chunk`. See the version of the function without `cell_chunk` for details.

# Positional arguments
* `rng`: the random number generator
* `grid`: the grid into which the particles are emitted
* `bc_list`: the `Tuple` of boundary conditions (left and right wall)
* `particles`: the `ParticleVector` of particles
* `pia`: the `ParticleIndexerArray` instance
* `species`: the index of the species being emitted
* `species_data`: the vector of `Species` data
* `Fnum`: the computational weight of the emitted particles
* `Δt`: the timestep
* `cell_chunk`: the list of cell indices or range of the chunk of cells, should be ordered in increasing order

# Returns
The number of particles added.
"""
function emit_particles!(rng, grid::Grid1DUniform, bc_list, particles::ParticleVector{D}, pia, species, species_data,
                         Fnum, Δt, cell_chunk) where D
    return emit_particles!(rng, grid, bc_list, particles, pia, species, species_data, nothing, Fnum, Δt, cell_chunk)
end

"""
    emit_particles!(rng, grid::Grid1DUniform, bc_list, particles::ParticleVector{D}, pia, species, species_data, surf_props::SurfProps, Fnum, Δt) where D

Emit uniformly weighted particles across the [`InflowBC1D`](@ref) boundaries of a 1-D uniform grid over a timestep,
recording the emitted particles as reflected from the boundary they enter across, as well as any
impacts during their residual flight.
The contributions are added already scaled, so the function is to be called after
the convection routine that cleared and scaled `surf_props`, and before the surface properties are averaged.
See the version of the function without `surf_props` for details.

# Positional arguments
* `rng`: the random number generator
* `grid`: the grid into which the particles are emitted
* `bc_list`: the `Tuple` of boundary conditions (left and right wall)
* `particles`: the `ParticleVector` of particles
* `pia`: the `ParticleIndexerArray` instance
* `species`: the index of the species being emitted
* `species_data`: the vector of `Species` data
* `surf_props`: the `SurfProps` struct where the computed surface properties will be stored
* `Fnum`: the computational weight of the emitted particles
* `Δt`: the timestep

# Returns
The number of particles added.
"""
function emit_particles!(rng, grid::Grid1DUniform, bc_list, particles::ParticleVector{D}, pia, species, species_data,
                         surf_props::SurfProps, Fnum, Δt) where D
    return emit_particles!(rng, grid, bc_list, particles, pia, species, species_data, surf_props, Fnum, Δt, 1:grid.n_cells)
end

"""
    emit_particles!(rng, grid::Grid1DUniform, bc_list, particles::ParticleVector{D}, pia, species, species_data, surf_props::SurfProps, Fnum, Δt, cell_chunk) where D

Emit uniformly weighted particles across the [`InflowBC1D`](@ref) boundaries of a 1-D uniform grid over a timestep for a chunk of cells,
recording the emitted particles as reflected from the boundary they enter across, as well as any
impacts during their residual flight.
The contributions are added already scaled, so the function is to be called after
the convection routine that cleared and scaled `surf_props`, and before the surface properties are averaged.
See the version of the function without `surf_props` for details.

# Positional arguments
* `rng`: the random number generator
* `grid`: the grid into which the particles are emitted
* `bc_list`: the `Tuple` of boundary conditions (left and right wall)
* `particles`: the `ParticleVector` of particles
* `pia`: the `ParticleIndexerArray` instance
* `species`: the index of the species being emitted
* `species_data`: the vector of `Species` data
* `surf_props`: the `SurfProps` struct where the computed surface properties will be stored
* `Fnum`: the computational weight of the emitted particles
* `Δt`: the timestep
* `cell_chunk`: the list of cell indices or range of the chunk of cells, should be ordered in increasing order

# Returns
The number of particles added.
"""
function emit_particles!(rng, grid::Grid1DUniform, bc_list, particles::ParticleVector{D}, pia, species, species_data,
                         surf_props::Union{SurfProps,Nothing}, Fnum, Δt, cell_chunk) where D
    @inbounds mass = species_data[species].mass
    n_emitted = 0

    if first(cell_chunk) == 1
        @inbounds n_emitted += emit_particles_bc!(rng, grid, bc_list, bc_list[1], particles, pia, species,
                                                  surf_props, mass, Fnum, 1, Δt)
    end

    if last(cell_chunk) == grid.n_cells
        @inbounds n_emitted += emit_particles_bc!(rng, grid, bc_list, bc_list[2], particles, pia, species,
                                                  surf_props, mass, Fnum, 2, Δt)
    end

    return n_emitted
end

end