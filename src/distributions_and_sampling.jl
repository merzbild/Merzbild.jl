using Random

export tester

struct UnitDVGrid
    nx::Int16
    ny::Int16
    nz::Int16
    dx::Float64
    dy::Float64
    dz::Float64
    vx_grid::Array{Float64}
    vy_grid::Array{Float64}
    vz_grid::Array{Float64}
end

mutable struct DVGrid
    const base_grid::UnitDVGrid
    vx_max::Float64
    vy_max::Float64
    vz_max::Float64
    dx::Float64
    dy::Float64
    dz::Float64
    vx_grid::Array{Float64}
    vy_grid::Array{Float64}
    vz_grid::Array{Float64}
end

mutable struct VDF
    const nx::Int16
    const ny::Int16
    const nz::Int16
    w::Array{Float64,3}
end

# generate a grid with extent [-1.0,1.0]x[-1.0,1.0]x[-1.0,1.0]
function generate_unit_dvgrid(nx, ny, nz)
    nothing
end

# generate a grid [-vx_max, vx_max]x[-vy_max, vy_max]x[-vz_max, vz_max]
function generate_noiseless_dvgrid(nx, ny, nz, vx_max, vy_max, vz_max)
    nothing
end

function maxwellian(T, m, vx, vy, vz)
    return (m / (2.0 * π * k_B * T))^(1.5) * exp(-m * (vx^2 + vy^2 + vz^2) / (2.0 * k_B * T))
end

function bkw(T, m, vx, vy, vz)
    nothing
end

function evaluate_distribution(distribution_function, T, m, vx, vy, vz)
    return distribution_function(T, m, vx, vy, vz)
end

# compute thermal velocity: sqrt(2kT/m)
function compute_thermal_velocity(T, m)
    return sqrt(2 * k_B * T / m)
end

function sample_maxwellian!(rng, v, T, m)
    vscale = sqrt(2 * k_B * T / m)  # TODO: fix/check!
    vn = vscale * sqrt(-log(rand(rng, Float64)))
    vr = vscale * sqrt(-log(rand(rng, Float64)))
    theta1 = 2 * π * rand(rng, Float64)
    theta2 = 2 * π * rand(rng, Float64)

    v[1] = vn * cos(theta1)
    v[2] = vr * cos(theta2)
    v[3] = vr * sin(theta2)
end

function sample_particles_equal_weight!(rng, particles, nparticles, T, m, Fnum, xlo, xhi, ylo, yhi, zlo, zhi)
    for i in 1:nparticles
        particles[i] = Particle(Fnum , [0.0, 0.0, 0.0], [xlo + rand(rng, Float64) * (xhi - xlo),
                                                         ylo + rand(rng, Float64) * (yhi - ylo),
                                                         zlo + rand(rng, Float64) * (zhi - zlo)])
        sample_maxwellian!(rng, particles[i].v, T, m)
    end
end