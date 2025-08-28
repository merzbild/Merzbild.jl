using Random
import Distributions

@muladd begin

"""
    UnitDVGrid

Stores information about a uniform discrete 3-dimensional velocity grid with extent ``[-1.0,1.0]\\times[-1.0,1.0]\\times[-1.0,1.0]``.

# Fields
* `nx`: number of cells in x direction
* `ny`: number of cells in y direction
* `nz`: number of cells in z direction
* `dx`: grid spacing in x direction
* `dy`: grid spacing in y direction
* `dz`: grid spacing in z direction
* `vx_grid`: `Vector` of grid nodes in x direction
* `vy_grid`: `Vector` of grid nodes in y direction
* `vz_grid`: `Vector` of grid nodes in z direction
"""
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

"""
    DVGrid

Stores information about a uniform discrete 3-dimensional symmetric
velocity grid with extent ``[-v_{x,max},v_{x,max}]\\times[-v_{y,max},v_{y,max}]\\times[-v_{z,max},v_{z,max}]``.

# Fields
* `base_grid`: the underlying unit (non-scaled) `UnitDVGrid` uniform grid 
* `vx_max`: extent of the grid in the x direction
* `vy_max`: extent of the grid in the y direction
* `vz_max`: extent of the grid in the z direction
* `dy`: grid spacing in x direction
* `dy`: grid spacing in y direction
* `dz`: grid spacing in z direction
* `vx_grid`: `Vector` of grid nodes in x direction
* `vy_grid`: `Vector` of grid nodes in y direction
* `vz_grid`: `Vector` of grid nodes in z direction
"""
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

"""
    VDF

Stores the values of a function evaluated on 3-dimensional velocity grid.

# Fields
* `nx`: the number of grid elements in the x direction
* `ny`: the number of grid elements in the y direction
* `nz`: the number of grid elements in the z direction
* `w`: the 3-dimensional array of values
"""
mutable struct VDF
    const nx::Int16
    const ny::Int16
    const nz::Int16
    w::Array{Float64,3}
end

"""
    UnitDVGrid(nx, ny, nz)

Generate a velocity grid with extent `[-1.0,1.0]x[-1.0,1.0]x[-1.0,1.0]` of type `UnitDVGrid`.

# Positional arguments
* `nx`: the number of grid elements in the x direction
* `ny`: the number of grid elements in the y direction
* `nz`: the number of grid elements in the z direction
"""
function UnitDVGrid(nx, ny, nz)
    vx = Vector(LinRange(-1.0, 1.0, nx))
    vy = Vector(LinRange(-1.0, 1.0, ny))
    vz = Vector(LinRange(-1.0, 1.0, nz))
    return UnitDVGrid(nx, ny, nz,
                      vx[2] - vx[1], vy[2] - vy[1], vz[2] - vz[1],
                      vx, vy, vz)
end

"""
    DVGrid(nx, ny, nz, vx_max, vy_max, vz_max)

Generate a symmetric velocity grid with extent
``[-v_{x,max},v_{x,max}]\\times[-v_{y,max},v_{y,max}]\\times[-v_{z,max},v_{z,max}]``.

# Positional arguments
* `nx`: the number of grid elements in the x direction
* `ny`: the number of grid elements in the y direction
* `nz`: the number of grid elements in the z direction
* `vx_max`: extent of the grid in the x direction
* `vy_max`: extent of the grid in the y direction
* `vz_max`: extent of the grid in the z direction
"""
function DVGrid(nx, ny, nz, vx_max, vy_max, vz_max)
    unitgrid = UnitDVGrid(nx, ny, nz)
    return DVGrid(unitgrid, vx_max, vy_max, vz_max,
                  unitgrid.dx * vx_max, unitgrid.dy * vy_max, unitgrid.dz * vz_max,
                  unitgrid.vx_grid * vx_max, unitgrid.vy_grid * vy_max, unitgrid.vz_grid * vz_max)
end

"""
    VDF(nx, ny, nz)

Create an empty `VDF` instance of size `nx * ny * nz`.

# Positional arguments
* `nx`: the number of grid elements in the x direction
* `ny`: the number of grid elements in the y direction
* `nz`: the number of grid elements in the z direction
"""
function VDF(nx, ny, nz)
    return VDF(nx, ny, nz, zeros(nx, ny, nz))
end

"""
    maxwellian(vx, vy, vz, m, T)

Evaluate the Maxwell distribution with temperature `T` for a species with mass `m`
    at a velocity `(vx, vy, vz)`

# Positional arguments
* `vx`: x velocity
* `vy`: y velocity
* `vz`: z velocity
* `m`: species' mass
* `T`: temperature
"""
function maxwellian(vx, vy, vz, m, T)
    return (m / (2.0 * π * k_B * T))^(1.5) * exp(-m * (vx^2 + vy^2 + vz^2) / (2.0 * k_B * T))
end

"""
    bkw(vx, vy, vz, m, T, scaled_time)

Evaluate the Bobylev-Krook-Wu (BKW) distribution with temperature `T` for a species with mass `m`
    at a velocity `(vx, vy, vz)` and scaled time `scaled_time`

# Positional arguments
* `vx`: x velocity
* `vy`: y velocity
* `vz`: z velocity
* `m`: species' mass
* `T`: temperature
* `scaled_time`: the scaled_time
"""
function bkw(vx, vy, vz, m, T, scaled_time)
    xk = 1.0 - 0.4 * exp(-scaled_time / 6.0)

    # bkw = norm * ( 5.d0 * xk - 3.d0 + 2.d0 * ( 1.d0 - xk ) * Csq * mass / ( xk * temp ) ) * &
    # exp( -Csq * mass / ( xk * temp ) )


    Csq = (vx^2 + vy^2 + vz^2)
    return (5 * xk - 3 + 2 * (1.0 - xk) * Csq * m / (2 * k_B * xk * T)) * exp(-Csq * m / (2 * k_B * xk * T))
end

"""
    sample_bkw!(rng, particles, nparticles, offset, m, T, v0)

Sample particle velocities from the BKW distribution with temperature `T` for a species with mass `m`
at `t=0` and add a velocity offset.
Note: This does not update the particle weights, positions, or any indexing structures.

# Positional arguments
* `rng`: the random number generator
* `particles`: the `Vector`-like structure holding the particles
* `nparticles`: the number of particles to sample
* `offset`: offset the starting position for writing the sampled particles to the `particles` array 
* `m`: species' mass
* `T`: temperature
* `v0`: the 3-dimensional velocity to add to the sampled velocities
"""
function sample_bkw!(rng, particles, nparticles, offset, m, T, v0)
    # BKW at t=0
    vscale = sqrt(2 * k_B * T / m) * sqrt(0.3)  # 0.3 comes from some scaling of the Chi distribution

    v_distribution = Distributions.Chi(5) 
    v_abs = rand(rng, v_distribution, nparticles)

    Θ = rand(rng, Float64, nparticles) * π
    ϕ = rand(rng, Float64, nparticles) * twopi
    sintheta = sin.(Θ)

    vx = v_abs .* sintheta .* cos.(ϕ)
    vy = v_abs .* sintheta .* sin.(ϕ)
    vz = v_abs .* cos.(Θ)

    for i in 1:nparticles
        @inbounds particles[i+offset].v = vscale * SVector{3,Float64}(vx[i], vy[i], vz[i]) .+ v0
    end
end

"""
    sample_bkw!(rng, particles, nparticles, m, T, v0)

Sample particles' velocities from the BKW distribution with temperature `T` for a species with mass `m`
    at `t=0` and add a velocity offset. This does not update the particle weights, positions, or any indexing structures.
    The sampled particles are written to the start of the `particles` array.

# Positional arguments
* `rng`: the random number generator
* `particles`: the `Vector`-like structure holding the particles
* `nparticles`: the number of particles to sample
* `m`: species' mass
* `T`: temperature
* `v0`: the 3-dimensional velocity to add to the sampled velocities
"""
function sample_bkw!(rng, particles, nparticles, m, T, v0)
    sample_bkw!(rng, particles, nparticles, 0, m, T, v0)
end

"""
    evaluate_distribution_on_grid!(vdf, distribution_function, grid, w_total, cutoff_v; normalize=true)

Evaluate a distribution on a discrete velocity grid, considering only points inside a sphere of a given radius
(the value of the VDF at points outside of the sphere will be 0.0). The values can be re-normalized so that
the distribution has the prescribed computational weight/density.

# Positional arguments
* `vdf`: The `VDF` instance where the evaluated values will be stored
* `distribution_function`: the distribution function to be evaluated which takes the x, y, and z velocities as parameters
* `grid`: the `DVGrid` on which the distribution function is evaluated
* `w_total`: the total weight used in the re-normalization
* `cutoff_v`: the radius of the sphere used to cut off high velocities: on grid points outside the sphere
    the VDF will be 0

# Keyword arguments
* `normalize`: if `true`, the resulting values of the VDF will be renormalized so that their sum is equal to
    `w_total`
"""
function evaluate_distribution_on_grid!(vdf, distribution_function, grid, w_total, cutoff_v; normalize=true)
    w = 0.0
    for k in 1:grid.base_grid.nz
        for j in 1:grid.base_grid.ny
            for i in 1:grid.base_grid.nx
                @inbounds if (sqrt(grid.vx_grid[i]^2 + grid.vy_grid[j]^2 + grid.vz_grid[k]^2) <= cutoff_v)
                    @inbounds vdf.w[i,j,k] = distribution_function(grid.vx_grid[i], grid.vy_grid[j], grid.vz_grid[k])
                    @inbounds w += vdf.w[i,j,k]
                end
            end
        end
    end
    if normalize
        vdf.w = vdf.w * w_total / w
    end
end

"""
    sample_on_grid!(rng, vdf_func, particles, nv, m, T, n_total,
                    xlo, xhi, ylo, yhi, zlo, zhi; v_mult=3.5, cutoff_mult=3.5, noise=0.0,
                    v_offset=[0.0, 0.0, 0.0])

Sample particles by evaluating a distribution on a discrete velocity grid,
considering only points inside a sphere of a given radius
(the value of the VDF at points outside of the sphere will be 0.0). The values of the VDF at the grid points
will then be the computational weights of the particles, and the particles velocities are taken
to be the velocities of the corresponding grid nodes (with additional uniformly distributed noise).
Note: this can produce a large amount of particles for fine grids (as the number of grid nodes scales
as `nv^3`.)
The grid is assumed to have the same number of nodes `nv` in each direction, and the extent is computed
as `v_mult * v_thermal`, where `v_thermal` is the thermal velocity ``\\sqrt(2kT/m)``, and `v_mult` is a
user-defined parameter. The positions of the particles are assumed to be randomly distributed in a cuboid.

# Positional arguments
* `rng`: The random number generator
* `vdf_func`: the distribution function to be evaluated which takes the x, y, and z velocities as parameters
* `particles`: the `Vector`-like structure holding the particles
* `nv`: the number of grid nodes in each direction
* `m`: the molecular mass of the species
* `T`: the temperature used to compute the thermal velocity
* `n_total`: the total computational weight (number of physical particles) to be sampled
* `xlo`: the lower bound of the x coordinates of the particles
* `xhi`: the upper bound of the x coordinates of the particles
* `ylo`: the lower bound of the y coordinates of the particles
* `yhi`: the upper bound of the y coordinates of the particles
* `zlo`: the lower bound of the z coordinates of the particles
* `zhi`: the upper bound of the z coordinates of the particles

# Keyword arguments
* `v_mult`: the value by which the thermal velocity is multiplied to compute the extent of the velocity grid
* `cutoff_mult`: the value by which the thermal velocity is multiplied to compute the radius for the sphere 
    used to cut-off the higher velocities
* `noise`: controls the amount of noise added to the particle velocities (the noise is uniformly distributed
    on the interval `[-noise*dv, noise*dv]`, where `dv` is the grid spacing)
* `v_offset`: the streaming velocity vector to be added to the particle velocities

# Returns
* The number of particles created
"""
function sample_on_grid!(rng, vdf_func, particles, nv, m, T, n_total,
                         xlo, xhi, ylo, yhi, zlo, zhi; v_mult=3.5, cutoff_mult=3.5, noise=0.0,
                         v_offset=[0.0, 0.0, 0.0])

    vdf = VDF(nv, nv, nv)
    v_thermal = compute_thermal_velocity(m, T)
    v_grid = DVGrid(nv, nv, nv, v_thermal * v_mult, v_thermal * v_mult, v_thermal * v_mult)

    # maxwell_df = (vx,vy,vz) -> vdf_func(T, m, vx, vy, vz)

    evaluate_distribution_on_grid!(vdf, vdf_func, v_grid, n_total, v_thermal * cutoff_mult; normalize=true)

    pid = 0
    n_sampled = 0
    for k in 1:v_grid.base_grid.nz
        for j in 1:v_grid.base_grid.ny
            for i in 1:v_grid.base_grid.nx
                @inbounds if vdf.w[i,j,k] > 0.0
                    pid += 1
                    n_sampled += 1

                    @inbounds add_particle!(particles, pid, vdf.w[i,j,k],
                                  SVector{3}(v_grid.vx_grid[i] + noise * v_grid.dx * (0.5 - rand(rng, Float64)) + v_offset[1],
                                             v_grid.vy_grid[j] + noise * v_grid.dy * (0.5 - rand(rng, Float64)) + v_offset[2],
                                             v_grid.vz_grid[k] + noise * v_grid.dz * (0.5 - rand(rng, Float64)) + v_offset[3]),
                                  SVector{3}(xlo + rand(rng, Float64) * (xhi - xlo),
                                             ylo + rand(rng, Float64) * (yhi - ylo),
                                             zlo + rand(rng, Float64) * (zhi - zlo)))
                end
            end
        end
    end

    return n_sampled
end


"""
    sample_maxwellian_on_grid!(rng, particles, nv, m, T, n_total,
                               xlo, xhi, ylo, yhi, zlo, zhi; v_mult=3.5, cutoff_mult=3.5, noise=0.0,
                               v_offset=[0.0, 0.0, 0.0])

Sample particles by evaluating a Maxwellian on a discrete velocity grid,
considering only points inside a sphere of a given radius
(the value of the VDF at points outside of the sphere will be 0.0). The values of the VDF at the grid points
will then be the computational weights of the particles, and the particles velocities are taken
to be the velocities of the corresponding grid nodes (with additional uniformly distributed noise).
Note: this can produce a large amount of particles for fine grids (as the number of grid nodes scales
as `nv^3`.)
The grid is assumed to have the same number of nodes `nv` in each direction, and the extent is computed
as `v_mult * v_thermal`, where `v_thermal` is the thermal velocity ``\\sqrt(2kT/m)``, and `v_mult` is a
user-defined parameter. The positions of the particles are assumed to be randomly distributed in a cuboid.

# Positional arguments
* `rng`: The random number generator
* `particles`: the `Vector`-like structure holding the particles
* `nv`: the number of grid nodes in each direction
* `m`: the molecular mass of the species
* `T`: the temperature used to compute the thermal velocity
* `n_total`: the total computational weight (number of physical particles) to be sampled
* `xlo`: the lower bound of the x coordinates of the particles
* `xhi`: the upper bound of the x coordinates of the particles
* `ylo`: the lower bound of the y coordinates of the particles
* `yhi`: the upper bound of the y coordinates of the particles
* `zlo`: the lower bound of the z coordinates of the particles
* `zhi`: the upper bound of the z coordinates of the particles

# Keyword arguments
* `v_mult`: the value by which the thermal velocity is multiplied to compute the extent of the velocity grid
* `cutoff_mult`: the value by which the thermal velocity is multiplied to compute the radius for the sphere 
    used to cut-off the higher velocities
* `noise`: controls the amount of noise added to the particle velocities (the noise is uniformly distributed
    on the interval `[-noise*dv, noise*dv]`, where `dv` is the grid spacing)
* `v_offset`: the streaming velocity vector to be added to the particle velocities

# Returns
The function returns the number of particles created
"""
function sample_maxwellian_on_grid!(rng, particles, nv, m, T, n_total,
                                    xlo, xhi, ylo, yhi, zlo, zhi; v_mult=3.5, cutoff_mult=3.5, noise=0.0,
                                    v_offset=[0.0, 0.0, 0.0])

    maxwell_df = (vx,vy,vz) -> maxwellian(vx, vy, vz, m, T)
    return sample_on_grid!(rng, maxwell_df, particles, nv, m, T, n_total,
                           xlo, xhi, ylo, yhi, zlo, zhi; v_mult=v_mult, cutoff_mult=cutoff_mult, noise=noise,
                           v_offset=v_offset)
end

"""
    compute_thermal_velocity(m, T)

Compute the thermal velocity ``\\sqrt(2kT/m)``.

# Positional arguments
* `m`: the molecular mass of the species
* `T`: the temperature

# Returns
The thermal velocity
"""
function compute_thermal_velocity(m, T)
    return sqrt(2 * k_B * T / m)
end

"""
    sample_maxwellian!(rng, particles, nparticles, offset, m, T, v0)

Sample `nparticles` particles from a Maxwellian with temperature T for a species with mass m
and add a velocity offset.
Note: This does not update the particle weights, positions, or any indexing structures.

# Positional arguments
* `rng`: the random number generator
* `particles`: the `Vector`-like structure holding the particles
* `nparticles`: the number of particles to sample
* `offset`: offset the starting position for writing the sampled particles to the `particles` array 
* `m`: species' mass
* `T`: temperature
* `v0`: the 3-dimensional velocity to add to the sampled velocities
"""
function sample_maxwellian!(rng, particles, nparticles, offset, m, T, v0)
    vscale = compute_thermal_velocity(m, T)

    for i in 1:nparticles
        vn = sqrt(-log(rand(rng, Float64)))
        vr = sqrt(-log(rand(rng, Float64)))
        theta1 = twopi * rand(rng, Float64)
        theta2 = twopi * rand(rng, Float64)

        @inbounds particles[i+offset].v = vscale * SVector{3,Float64}(vn * cos(theta1), vr * cos(theta2), vr * sin(theta2)) + v0
    end
end

"""
    sample_particles_equal_weight!(rng, particles, pia, cell, species,
                                        nparticles, m, T, Fnum, xlo, xhi, ylo, yhi, zlo, zhi;
                                        distribution=:Maxwellian, vx0=0.0, vy0=0.0, vz0=0.0)

Sample equal-weight particles of a specific species in a specific cell from a distribution.
The positions of the particles are assumed to be randomly distributed in a cuboid.
Note: this does not work if applied twice in a row to the same cell.

# Positional arguments
* `rng`: the random number generator
* `particles`: the `Vector`-like structure holding the particles
* `pia`: the `ParticleIndexerArray`
* `cell`: the cell index
* `species`: the species index
* `nparticles`: the number of particles to sample
* `m`: species' mass
* `T`: temperature
* `Fnum`: the computational weight of the particles
* `xlo`: the lower bound of the x coordinates of the particles
* `xhi`: the upper bound of the x coordinates of the particles
* `ylo`: the lower bound of the y coordinates of the particles
* `yhi`: the upper bound of the y coordinates of the particles
* `zlo`: the lower bound of the z coordinates of the particles
* `zhi`: the upper bound of the z coordinates of the particles

# Keyword arguments
* `distribution`: the distribution to sample from (either `:Maxwellian` or `:BKW`)
* `vx0`: the x-velocity offset to add to the particle velocities
* `vy0`: the y-velocity offset to add to the particle velocities
* `vz0`: the z-velocity offset to add to the particle velocities
"""
function sample_particles_equal_weight!(rng, particles, pia, cell, species,
                                        nparticles, m, T, Fnum, xlo, xhi, ylo, yhi, zlo, zhi;
                                        distribution=:Maxwellian, vx0=0.0, vy0=0.0, vz0=0.0)

    @inbounds start = pia.n_total[species] + 1                                    
    @inbounds pia.indexer[cell, species].n_local = nparticles
    @inbounds pia.n_total[species] += nparticles

    @inbounds pia.indexer[cell, species].start1 = start
    @inbounds pia.indexer[cell, species].end1 = start - 1 + nparticles
    @inbounds pia.indexer[cell, species].n_group1 = nparticles

    @inbounds pia.indexer[cell, species].start2 = 0
    @inbounds pia.indexer[cell, species].end2 = -1
    @inbounds pia.indexer[cell, species].n_group2 = 0

    offset = start - 1

    for i in 1:nparticles
        add_particle!(particles, i+offset, Fnum,  SVector{3}(0.0, 0.0, 0.0),
                      SVector{3}(xlo + rand(rng, Float64) * (xhi - xlo),
                                 ylo + rand(rng, Float64) * (yhi - ylo),
                                 zlo + rand(rng, Float64) * (zhi - zlo)))
        @inbounds particles.cell[i+offset] = cell
    end

    v0 = SVector{3}(vx0, vy0, vz0)
    if distribution == :Maxwellian
        sample_maxwellian!(rng, particles, nparticles, offset, m, T, v0)
    elseif distribution == :BKW
        sample_bkw!(rng, particles, nparticles, m, T, v0)
    end
end

end