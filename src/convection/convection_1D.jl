function convect_particles!(rng, grid::Grid1DUniform, boundaries::MaxwellWalls, particles, pia, species, species_data, Δt)
    @inbounds @simd for i in 1:pia.n_total[species]
        t_rest = Δt
        x_old = particles[i].x[1]
        x_new = particles[i].x[1] + particles[i].v[1] * Δt

        while (x_new >= grid.L) || (x_new <= 0.0)
            
            if (x_new >= grid.L)
                t_rest -= abs((grid.L - x_old) / particles[i].v[1])
                bc_id = 2
                wall_normal = -1
                x_old = grid.L
            else
                t_rest -= abs(x_old / particles[i].v[1])
                bc_id = 1
                wall_normal = 1
                x_old = 0.0
            end

            reflect_particle_x!(rng, particles[i], boundaries.reflection_velocities_sq[bc_id, species],
                                wall_normal,
                                boundaries.boundaries[bc_id].v,
                                boundaries.boundaries[bc_id].accomodation)

            x_new = x_old + particles[i].v[1] * t_rest
        end

        # if a particle is too near a wall, we offset it a bit to avoid particles
        # that are exactly at a wall screwing up counters, etc.
        if x_new < grid.min_x
            x_new = grid.min_x
        elseif x_new > grid.max_x
            x_new = grid.max_x
        end

        particles[i].x = SVector{3,Float64}(x_new, particles[i].x[2], particles[i].x[3])
    end
end