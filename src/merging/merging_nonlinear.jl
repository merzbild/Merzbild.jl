function compute_nonlinear_lhs!(lhs, sol_vector, n_particles, v_moments, x_moments)
    lhs .= 0.0
    for i in 1:n_particles
        lhs[1] += sol_vector[i]  # weight constraint, should sum up to 1

        for j in 1:length(v_moments)
            lhs[1+j] += sol_vector[i] * (sol_vector[n_particles + i]^v_moments[j])
        end
        for j in 1:length(x_moments)
            lhs[1+length(v_moments)+j] += sol_vector[i] * (sol_vector[n_particles + length(v_moments) + i]^x_moments[j])
        end
    end
end

function compute_jacobian_sq!(J, sol_vector, n_particles, v_moments, x_moments)
    # js = 1 * n_particles + 3 * n_particles + 3 * n_particles
    # J_ij = ∂F_i/∂x_j
    J[1:end, 1:n_particles] .= 1.0  # everything is always linear in the weights
    nvm = length(v_moments)
    nxm = length(x_moments)

    # velocity moments equations, derivatives w.r.t. velocity
    for i in 1:nvm
        for j in 1:nvm
            # equation for conservation of moment v_moments[i]
            J[1+i, n_particles + j] = sol_vector[i] * v_moments[i] * sol_vector[n_particles + j]^v_moments[i]
        end
    end

    # position moments equations, derivatives w.r.t. position
    for i in 1:nxm
        for j in 1:nxm
            # equation for conservation of moment v_moments[i]
            J[1+nvm+i, n_particles + nvm + j] = sol_vector[i] * x_moments[i] * sol_vector[n_particles + nvm + j]^x_moments[i]
        end
    end
end