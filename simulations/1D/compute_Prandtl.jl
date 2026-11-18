# include("../../src/Merzbild.jl")

using Merzbild
using Printf
using Random

isdir("scratch/data") || error("Output directory `scratch/data` does not exist (the path is relative to the working directory). Create it, e.g. with `mkdir -p scratch/data`.")

"""
    run_flow(seed, T_bg0, T_wall1, T_wall2, v_wall, L, ndens, nx,
             ppc,
             Δt, n_timesteps, avg_start,
             name; reorder_freq = 10)

Run a single fixed-weight 1-D flow between two fully diffuse walls and time-average the physical
and flux properties. Setting `v_wall` to 0 with `T_wall1 != T_wall2` gives a Fourier flow; setting
`T_wall1 == T_wall2` with a non-zero `v_wall` gives a Couette flow.
Time-averaged output is written to `scratch/data/Prandtl_[name]_[nx]_[n_timesteps].nc` and
`scratch/data/Prandtl_[name]_[nx]_[n_timesteps]_flux.nc`.

Positional arguments:
* `seed`: random seed value
* `T_bg0`: initial temperature of gas
* `T_wall1`: temperature of left wall
* `T_wall2`: temperature of right wall
* `v_wall`: left wall will have y-velocity of `-v_wall`, right wall will have y-velocity of `v_wall`
* `L`: length of domain in meters
* `ndens`: initial number density
* `nx`: number of grid cells
* `ppc`: number of particles sampled in each cell
* `Δt`: timestep size
* `n_timesteps`: total number of timesteps to run for
* `avg_start`: timestep after which averaging begins
* `name`: name to append to the output files

Keyword arguments:
* `reorder_freq`: how frequently the ordering of particle indices is restored

Returns:
A named tuple `(grid, species_data, interaction_data, phys_props_avg, flux_props_avg)`.
"""
function run_flow(seed, T_bg0, T_wall1, T_wall2, v_wall, L, ndens, nx,
                  ppc,
                  Δt, n_timesteps, avg_start,
                  name; reorder_freq = 10)
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    particles_data_path = joinpath(MERZBILD_DATA_PATH, "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")
    interaction_data_path = joinpath(MERZBILD_DATA_PATH, "vhs.toml")
    interaction_data::Array{Interaction, 2} = load_interaction_data(interaction_data_path, species_data)

    grid = Grid1DUniform(L, nx)
    bc_list = (FullyDiffuseBC1D(1, species_data, T_wall1, [0.0, -v_wall, 0.0]),
               FullyDiffuseBC1D(1, species_data, T_wall2, [0.0, v_wall, 0.0]))

    T_mean = 0.5 * (T_wall1 + T_wall2)
    mass = species_data[1].mass
    λ = mean_free_path(interaction_data, 1, ndens, T_mean)
    ν = mean_collision_frequency(interaction_data, 1, species_data, ndens, T_mean)
    v_thermal = sqrt(2 * k_B * T_mean / mass)
    μ_mean = interaction_data[1, 1].vhs_muref * (T_mean / interaction_data[1, 1].vhs_Tref)^interaction_data[1, 1].vhs_o
    t_diffusion = L^2 * ndens * mass / μ_mean

    println("--- $(name) flow: T_wall = ($(T_wall1), $(T_wall2)) K, v_wall = $(v_wall) m/s ---")
    @printf("  λ = %.4e m, Kn = λ/L = %.4e, λ/Δx = %.3f\n", λ, λ / L, λ / grid.Δx)
    @printf("  Δt/τ_coll = %.4f, Δt/(Δx/v_thermal) = %.4f\n", Δt * ν, Δt * v_thermal / grid.Δx)
    @printf("  diffusion time L²ρ/μ = %.4e s, transient = %.4e s (%.2f diffusion times)\n",
            t_diffusion, avg_start * Δt, avg_start * Δt / t_diffusion)

    n_particles = ppc * nx
    particles = [ParticleVector{1}(n_particles)]
    pia = ParticleIndexerArray(grid.n_cells, 1)
    gridsorter = GridSortInPlace(grid, n_particles)

    # Fnum * ppc = Np in cell = ndens * V_cell
    Fnum = grid.cells[1].V * ndens / ppc
    sample_particles_equal_weight!(rng, grid, particles[1], pia, 1,
                                   species_data, ndens, T_bg0, Fnum)

    collision_data = CollisionData()
    collision_factors = create_collision_factors_array(pia, interaction_data, species_data, T_bg0, Fnum)

    phys_props = PhysProps(pia; ndens_not_Np=true)
    phys_props_avg = PhysProps(pia; ndens_not_Np=true)
    flux_props = FluxProps(pia)
    flux_props_avg = FluxProps(pia)

    ds_avg = NCDataHolder("scratch/data/Prandtl_$(name)_$(nx)_$(n_timesteps).nc",
                          species_data, phys_props_avg)
    ds_flux_avg = NCDataHolderFlux("scratch/data/Prandtl_$(name)_$(nx)_$(n_timesteps)_flux.nc",
                                   species_data, flux_props_avg)
    write_grid("scratch/data/Prandtl_$(nx)_grid.nc", grid)

    n_avg = n_timesteps - avg_start + 1
    index_inv_map = zeros(Int64, n_particles)

    for t in 1:n_timesteps
        if t % 20000 == 0
            println("  t = $t")
        end

        for cell in 1:grid.n_cells
            ntc_equal_weight!(rng, collision_factors[1, 1, cell],
                              collision_data, interaction_data, particles[1],
                              pia, cell, 1, Δt, grid.cells[cell].V)
        end

        convect_particles_and_compute_cell!(rng, grid, bc_list, particles[1], pia, 1, species_data, Δt)
        sort_particles!(gridsorter, particles[1], pia, 1)

        if t % reorder_freq == 0
            restore_particle_ordering!(particles[1], index_inv_map)
        end

        if (t >= avg_start)
            compute_props_sorted!(particles, pia, species_data, phys_props, grid)
            compute_flux_props_sorted!(particles, pia, species_data, phys_props, flux_props, grid)
            avg_props!(phys_props_avg, phys_props, n_avg)
            avg_props!(flux_props_avg, flux_props, n_avg)
        end
    end

    write_netcdf(ds_avg, phys_props_avg, n_timesteps)
    write_netcdf(ds_flux_avg, flux_props_avg, n_timesteps)

    close_netcdf(ds_avg)
    close_netcdf(ds_flux_avg)

    return (grid=grid, species_data=species_data, interaction_data=interaction_data,
            phys_props_avg=phys_props_avg, flux_props_avg=flux_props_avg)
end

"""
    interior_cell_range(grid, interaction_data, ndens, T_mean, n_knudsen_layers)

Compute the range of cells lying further than `n_knudsen_layers` mean free paths away from either
wall, i.e. the cells in which the Navier-Stokes-level constitutive relations are expected to hold.

# Positional arguments
* `grid`: the physical grid
* `interaction_data`: the 2-dimensional array of `Interaction` instances
* `ndens`: the number density used to estimate the mean free path
* `T_mean`: the temperature used to estimate the mean free path
* `n_knudsen_layers`: number of mean free paths to skip at each wall

# Returns
The range of interior cell indices.
"""
function interior_cell_range(grid, interaction_data, ndens, T_mean, n_knudsen_layers)
    λ = mean_free_path(interaction_data, 1, ndens, T_mean)
    n_skip = max(1, ceil(Int64, n_knudsen_layers * λ / grid.Δx))

    if grid.n_cells - 2 * n_skip < 5
        error("Too few cells outside the Knudsen layers ($(grid.n_cells - 2 * n_skip) cells for " *
              "$(grid.n_cells) cells in total, skipping $(n_skip) cells at each wall); " *
              "increase the number of cells or reduce the Knudsen number")
    end

    return (n_skip + 1):(grid.n_cells - n_skip)
end

"""
    transport_coefficient(profile, flux, T, Δx, cell_range, T_ref, vhs_o)

Extract a transport coefficient from a time-averaged flux and the gradient of the corresponding
time-averaged profile, ``\\zeta_i = -F_i / (\\partial_x f)_i``, where the gradient is computed using
central differences. As the VHS transport coefficients scale as ``T^\\omega``, each per-cell value is
scaled to the VHS reference temperature by ``(T_{ref}/T_i)^\\omega`` before averaging, which removes
the influence of the temperature variation across the domain.

# Positional arguments
* `profile`: the cell-centered profile the gradient of which drives the flux (temperature or velocity)
* `flux`: the cell-centered flux (heat flux or shear stress)
* `T`: the cell-centered temperature
* `Δx`: the cell size
* `cell_range`: the range of cells over which to average
* `T_ref`: the VHS reference temperature
* `vhs_o`: the VHS exponent (omega)

# Returns
Tuple of the mean transport coefficient at `T_ref` and its standard deviation across `cell_range`.
"""
function transport_coefficient(profile, flux, T, Δx, cell_range, T_ref, vhs_o)
    inv_2Δx = 0.5 / Δx
    mean = 0.0
    mean_sq = 0.0
    n = 0

    @inbounds for cell in cell_range
        gradient = (profile[cell + 1] - profile[cell - 1]) * inv_2Δx
        ζ = -flux[cell] * (T_ref / T[cell])^vhs_o / gradient
        mean += ζ
        mean_sq += ζ * ζ
        n += 1
    end

    mean /= n
    return mean, sqrt(max(0.0, mean_sq / n - mean * mean))
end

"""
    relative_flux_variation(flux, cell_range)

Compute the standard deviation of a flux across a range of cells relative to its mean. At steady
state both the heat flux and the shear stress are constant across the domain, so a large value
indicates that the flow has not converged or that the statistics are insufficient.

# Positional arguments
* `flux`: the cell-centered flux
* `cell_range`: the range of cells over which to compute the variation

# Returns
The standard deviation of the flux divided by its mean.
"""
function relative_flux_variation(flux, cell_range)
    mean = 0.0
    mean_sq = 0.0
    n = 0

    @inbounds for cell in cell_range
        mean += flux[cell]
        mean_sq += flux[cell] * flux[cell]
        n += 1
    end

    mean /= n
    return sqrt(max(0.0, mean_sq / n - mean * mean)) / abs(mean)
end

"""
    run(nx, n_timesteps; kwargs...)

Compute the Prandtl number ``Pr = c_p \\mu / \\kappa`` of a monatomic VHS gas from two DSMC
simulations: the viscosity is extracted from the shear stress and the velocity gradient of a Couette
flow, and the thermal conductivity from the heat flux and the temperature gradient of a Fourier flow,
both taken from the `FluxProps` and `PhysProps` computed outside of the Knudsen layers.
Both coefficients are scaled to the VHS reference temperature, and are compared to the
Chapman-Enskog values. The first approximation of the Chapman-Enskog theory gives ``Pr = 2/3``
for a monatomic gas.

# Positional arguments
* `nx`: number of grid cells
* `n_timesteps`: total number of timesteps to run each of the two flows for

# Keyword arguments
* `seed`: random seed value
* `ndens`: number density
* `L`: length of domain in meters
* `ppc`: number of particles sampled in each cell
* `Δt`: timestep size
* `avg_start`: timestep after which averaging begins
* `T_cold`: temperature of the cold wall in the Fourier flow
* `T_hot`: temperature of the hot wall in the Fourier flow
* `T_couette`: temperature of both walls in the Couette flow
* `v_wall`: wall velocity in the Couette flow
* `n_knudsen_layers`: number of mean free paths skipped at each wall
* `reorder_freq`: how frequently the ordering of particle indices is restored
"""
function run(nx, n_timesteps;
             seed = 1234, ndens = 2.5e23, L = 5e-4, ppc = 200, Δt = 1e-9,
             avg_start = max(1, n_timesteps ÷ 3),
             T_cold = 300.0, T_hot = 600.0, T_couette = 300.0, v_wall = 500.0,
             n_knudsen_layers = 3, reorder_freq = 10)

    fourier = run_flow(seed, 0.5 * (T_cold + T_hot), T_cold, T_hot, 0.0, L, ndens, nx,
                       ppc, Δt, n_timesteps, avg_start, "fourier"; reorder_freq=reorder_freq)
    couette = run_flow(seed, T_couette, T_couette, T_couette, v_wall, L, ndens, nx,
                       ppc, Δt, n_timesteps, avg_start, "couette"; reorder_freq=reorder_freq)

    grid = fourier.grid
    interaction_data = fourier.interaction_data
    mass = fourier.species_data[1].mass
    T_ref = interaction_data[1, 1].vhs_Tref
    vhs_o = interaction_data[1, 1].vhs_o

    T_fourier = @view fourier.phys_props_avg.T[:, 1]
    heat_flux = @view fourier.flux_props_avg.kinetic_energy_flux[1, :, 1]
    cells_fourier = interior_cell_range(grid, interaction_data, ndens, 0.5 * (T_cold + T_hot), n_knudsen_layers)
    κ, κ_std = transport_coefficient(T_fourier, heat_flux, T_fourier, grid.Δx, cells_fourier, T_ref, vhs_o)

    T_couette_profile = @view couette.phys_props_avg.T[:, 1]
    vy = @view couette.phys_props_avg.v[2, :, 1]
    shear_stress = @view couette.flux_props_avg.off_diagonal_momentum_flux[1, :, 1]
    cells_couette = interior_cell_range(grid, interaction_data, ndens, T_couette, n_knudsen_layers)
    μ, μ_std = transport_coefficient(vy, shear_stress, T_couette_profile, grid.Δx, cells_couette, T_ref, vhs_o)

    c_p = 2.5 * k_B / mass
    Pr = c_p * μ / κ

    μ_theory = interaction_data[1, 1].vhs_muref
    κ_theory = 3.75 * (k_B / mass) * μ_theory
    Pr_theory = 2.0 / 3.0

    println("--- Prandtl number at T_ref = $(T_ref) K ---")
    @printf("  interior cells: Fourier %s, Couette %s (of %d cells)\n",
            cells_fourier, cells_couette, grid.n_cells)
    @printf("  relative variation of q_x  = %.4f\n", relative_flux_variation(heat_flux, cells_fourier))
    @printf("  relative variation of P_xy = %.4f\n", relative_flux_variation(shear_stress, cells_couette))
    @printf("  μ = %.6e ± %.2e Pa s (VHS: %.6e, %+.2f%%)\n",
            μ, μ_std, μ_theory, 100.0 * (μ / μ_theory - 1.0))
    @printf("  κ = %.6e ± %.2e W/(m K) (Chapman-Enskog: %.6e, %+.2f%%)\n",
            κ, κ_std, κ_theory, 100.0 * (κ / κ_theory - 1.0))
    @printf("  Pr = %.6f (Chapman-Enskog: %.6f, %+.2f%%)\n",
            Pr, Pr_theory, 100.0 * (Pr / Pr_theory - 1.0))

    return Pr
end

const n_t = 1000000
const n_cells = 200
run(n_cells, n_t)
# this setup gives
# μ = 2.212562e-05 ± 3.60e-07 Pa s (VHS: 2.177627e-05, +1.60%)
# κ = 1.733018e-02 ± 1.13e-03 W/(m K) (Chapman-Enskog: 1.700531e-02, +1.91%)
# Pr = 0.664664 (Chapman-Enskog: 0.666667, -0.30%)

# const n_t = 4500000
# const n_cells = 600
# run(n_cells, n_t; ppc=300)
# this setup gives
# μ = 2.217103e-05 ± 2.43e-07 Pa s (VHS: 2.177627e-05, +1.81%)
# κ = 1.728503e-02 ± 7.95e-04 W/(m K) (Chapman-Enskog: 1.700531e-02, +1.64%)
# Pr = 0.667768 (Chapman-Enskog: 0.666667, +0.17%)