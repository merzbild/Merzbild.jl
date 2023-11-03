using .Constants: k_B

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

function maxwellian(m, T, vx, vy, vz)
    return (m / (2.0 * π * k_B * T))^(1.5) * exp(-m * (vx^2 + vy^2 + vz^2) / (2.0 * k_B * T))
end

function bkw(m, T, vx, vy, vz)
    nothing
end

function evaluate_distribution(distribution_function, m, T, vx, vy, vz)
    return distribution_function(m, T, vx, vy, vz)
end

# compute thermal velocity: sqrt(2kT/m)
function compute_thermal_velocity(T, m)
    return sqrt(2 * k_B * T / m)
end