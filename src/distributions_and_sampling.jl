using Random
import Distributions

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
function create_unit_dvgrid(nx, ny, nz)
    vx = Vector(LinRange(-1.0, 1.0, nx))
    vy = Vector(LinRange(-1.0, 1.0, ny))
    vz = Vector(LinRange(-1.0, 1.0, nz))
    return UnitDVGrid(nx, ny, nz,
                      vx[2] - vx[1], vy[2] - vy[1], vz[2] - vz[1],
                      vx, vy, vz)
end

# generate a grid [-vx_max, vx_max]x[-vy_max, vy_max]x[-vz_max, vz_max]
function create_noiseless_dvgrid(nx, ny, nz, vx_max, vy_max, vz_max)
    unitgrid = create_unit_dvgrid(nx, ny, nz)
    return DVGrid(unitgrid, vx_max, vy_max, vz_max,
                  unitgrid.dx * vx_max, unitgrid.dy * vy_max, unitgrid.dz * vz_max,
                  unitgrid.vx_grid * vx_max, unitgrid.vy_grid * vy_max, unitgrid.vz_grid * vz_max)
end

function create_vdf(nx, ny, nz)
    return VDF(nx, ny, nz, zeros(nx, ny, nz))
end

function maxwellian(T, m, vx, vy, vz)
    return (m / (2.0 * π * k_B * T))^(1.5) * exp(-m * (vx^2 + vy^2 + vz^2) / (2.0 * k_B * T))
end

function sample_bkw!(rng, particles, nparticles, T, m, v0)
    vscale = sqrt(2 * k_B * T / m) * sqrt(0.3)  # 0.3 comes from some scaling of the Chi distribution

    v_distribution = Distributions.Chi(5) 
    v_abs = rand(v_distribution, nparticles)

    Θ = rand(rng, Float64, nparticles) * π
    ϕ = rand(rng, Float64, nparticles) * twopi
    sintheta = sin.(Θ)

    vx = v_abs .* sintheta .* cos.(ϕ)
    vy = v_abs .* sintheta .* sin.(ϕ)
    vz = v_abs .* cos.(Θ)

    for i in 1:nparticles
        particles[i].v = vscale * SVector{3,Float64}(vx[i], vy[i], vz[i]) .+ v0
    end
end

function evaluate_distribution_on_grid!(vdf, distribution_function, grid, w_total, cutoff_v; normalize=true)
    w = 0.0
    for k in 1:grid.base_grid.nz
        for j in 1:grid.base_grid.ny
            for i in 1:grid.base_grid.nx
                if (sqrt(grid.vx_grid[i]^2 + grid.vy_grid[j]^2 + grid.vz_grid[k]^2) <= cutoff_v)
                    vdf.w[i,j,k] = distribution_function(grid.vx_grid[i], grid.vy_grid[j], grid.vz_grid[k])
                    w += vdf.w[i,j,k]
                end
            end
        end
    end
    if normalize
        vdf.w = vdf.w * w_total / w
    end
end

function sample_maxwellian_on_grid!(rng, particles, nv, T, m, n_total,
                                    xlo, xhi, ylo, yhi, zlo, zhi; v_mult=3.5, cutoff_mult=3.5, noise=0.0,
                                    v_offset=[0.0, 0.0, 0.0])

    vdf = create_vdf(nv, nv, nv)
    v_thermal = compute_thermal_velocity(T, m)
    v_grid = create_noiseless_dvgrid(nv, nv, nv, v_thermal * v_mult, v_thermal * v_mult, v_thermal * v_mult)

    maxwell_df = (vx,vy,vz) -> maxwellian(T, m, vx, vy, vz)

    evaluate_distribution_on_grid!(vdf, maxwell_df, v_grid, n_total, v_thermal * cutoff_mult; normalize=true)

    pid = 0
    n_sampled = 0
    for k in 1:v_grid.base_grid.nz
        for j in 1:v_grid.base_grid.ny
            for i in 1:v_grid.base_grid.nx
                if vdf.w[i,j,k] > 0.0
                    pid += 1
                    n_sampled += 1
                    particles[pid] = Particle(vdf.w[i,j,k],
                                            # add random noise to velocity
                                            SVector{3}(v_grid.vx_grid[i] + noise * v_grid.dx * (0.5 - rand(rng, Float64)) + v_offset[1],
                                            v_grid.vy_grid[j] + noise * v_grid.dy * (0.5 - rand(rng, Float64)) + v_offset[2],
                                            v_grid.vz_grid[j] + noise * v_grid.dz * (0.5 - rand(rng, Float64)) + v_offset[3]),
                                            SVector{3}(xlo + rand(rng, Float64) * (xhi - xlo),
                                                    ylo + rand(rng, Float64) * (yhi - ylo),
                                                    zlo + rand(rng, Float64) * (zhi - zlo)))
                end
            end
        end
    end

    return n_sampled
end

# compute thermal velocity: sqrt(2kT/m)
function compute_thermal_velocity(T, m)
    return sqrt(2 * k_B * T / m)
end

function sample_maxwellian_single!(rng, v, T, m, v0)
    vscale = sqrt(2 * k_B * T / m)
    vn = vscale * sqrt(-log(rand(rng, Float64)))
    vr = vscale * sqrt(-log(rand(rng, Float64)))
    theta1 = twopi * rand(rng, Float64)
    theta2 = twopi * rand(rng, Float64)

    v[1] = vn * cos(theta1) + v0[1]
    v[2] = vr * cos(theta2) + v0[2]
    v[3] = vr * sin(theta2) + v0[3]
end

function sample_maxwellian!(rng, particles, nparticles, T, m, v0)
    vscale = sqrt(2 * k_B * T / m)

    for i in 1:nparticles
        vn = sqrt(-log(rand(rng, Float64)))
        vr = sqrt(-log(rand(rng, Float64)))
        theta1 = twopi * rand(rng, Float64)
        theta2 = twopi * rand(rng, Float64)

        particles[i].v = vscale * SVector{3,Float64}(vn * cos(theta1), vr * cos(theta2), vr * sin(theta2)) + v0
    end
end

function sample_particles_equal_weight!(rng, particles, nparticles, T, m, Fnum, xlo, xhi, ylo, yhi, zlo, zhi;
                                        distribution=:Maxwellian, vx0=0.0, vy0=0.0, vz0=0.0)
    # other options: Dirac Delta (T = 0: single point)
    # other options: Dirac Delta (T = N: two points)
    for i in 1:nparticles
        particles[i] = Particle(Fnum,  SVector{3}(0.0, 0.0, 0.0),  SVector{3}(xlo + rand(rng, Float64) * (xhi - xlo),
                                                                              ylo + rand(rng, Float64) * (yhi - ylo),
                                                                              zlo + rand(rng, Float64) * (zhi - zlo)))
    end

    v0 = SVector{3}(vx0, vy0, vz0)
    if distribution == :Maxwellian
        sample_maxwellian!(rng, particles, nparticles, T, m, v0)
    elseif distribution == :BKW
        sample_bkw!(rng, particles, nparticles, T, m, v0)
    end
end