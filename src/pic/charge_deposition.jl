@muladd begin

"""
    deposit_charge!(grid::Grid1DUniform, particles::ParticleVector, pia, species, species_data, field_props)

Deposit the charge of the particles of a single species on the nodes of a 1-D uniform grid using
first-order (cloud-in-cell) weighting: a particle located at a distance ``\\xi \\Delta x`` from the
left node of its cell contributes a fraction ``1 - \\xi`` of its charge to the left node and a
fraction ``\\xi`` to the right node.

The deposited values are accumulated, so the charge density has to be cleared by a call to
[`clear_charge_density!`](@ref) or [`clear_props!`](@ref) before the first species is deposited,
and [`normalize_charge_density!`](@ref) has to be called exactly once after the last species has
been deposited in order to turn the deposited charge into a charge density.
Neutral species are skipped.

The particles are assumed to be sorted on the grid.
This routine is serial; for a multi-threaded simulation, see the `cell_chunk` version of the
routine and [`reduce_field_props!`](@ref).

# Positional arguments
* `grid`: the `Grid1DUniform` grid
* `particles`: the `ParticleVector` of the particles of the species being deposited
* `pia`: the `ParticleIndexerArray` instance
* `species`: the index of the species being deposited
* `species_data`: the `Vector` of `Species` data
* `field_props`: the `ElectrostaticFieldProps` instance in which the deposited charge is stored
"""
@inline function deposit_charge!(grid::Grid1DUniform, particles::ParticleVector, pia, species, species_data,
                                 field_props)
    deposit_charge!(grid, particles, pia, species, species_data, field_props, 1:grid.n_cells)
end

"""
    deposit_charge!(grid::Grid1DUniform, particles::ParticleVector, pia, species, species_data, field_props, cell_chunk)

Deposit the charge of the particles of a single species located in a subset of the cells of a 1-D
uniform grid on the nodes of the grid.

In a multi-threaded simulation, each thread should deposit the particles of the cells it owns into its
own `ElectrostaticFieldProps` instance; the per-thread instances are then summed into the global one
with [`reduce_field_props!`](@ref).

# Positional arguments
* `grid`: the `Grid1DUniform` grid
* `particles`: the `ParticleVector` of the particles of the species being deposited
* `pia`: the `ParticleIndexerArray` instance
* `species`: the index of the species being deposited
* `species_data`: the `Vector` of `Species` data
* `field_props`: the `ElectrostaticFieldProps` instance in which the deposited charge is stored
* `cell_chunk`: the list of cell indices or range in which the charge of the particles is deposited
"""
function deposit_charge!(grid::Grid1DUniform, particles::ParticleVector, pia, species, species_data,
                         field_props, cell_chunk)
    @inbounds q_s = q_e * species_data[species].charge

    if q_s == 0.0
        return
    end

    inv_Δx = grid.inv_Δx

    @inbounds for cell in cell_chunk
        shift = cell - 1

        s1 = pia.indexer[cell,species].start1
        e1 = pia.indexer[cell,species].end1

        for i in s1:e1
            p = particles[i]
            ξ = p.x[1] * inv_Δx - shift
            q_w = q_s * p.w

            field_props.charge_density[cell] += (1.0 - ξ) * q_w
            field_props.charge_density[cell+1] += ξ * q_w
        end

        s2 = pia.indexer[cell,species].start2
        if s2 > 0
            e2 = pia.indexer[cell,species].end2

            for i in s2:e2
                p = particles[i]
                ξ = p.x[1] * inv_Δx - shift
                q_w = q_s * p.w

                field_props.charge_density[cell] += (1.0 - ξ) * q_w
                field_props.charge_density[cell+1] += ξ * q_w
            end
        end
    end
end

"""
    normalize_charge_density!(poisson_solver::PoissonSolver1DUniform, field_props)

Turn the charge deposited on the nodes of a 1D uniform grid by [`deposit_charge!`](@ref) into a charge density,
by dividing it by the volume of the dual cell of each node.
The dual cell of a boundary node is ``\\Delta x / 2``, so the deposited
charge in the boundary nodes is multiplied by 2.

Note that this routine has to be
called **exactly once** per timestep, after all the species have been deposited.
The solver is used only to determine the types of the boundary conditions and the cell size.

# Positional arguments
* `poisson_solver`: the `PoissonSolver1DUniform` instance
* `field_props`: the `ElectrostaticFieldProps` instance holding the deposited charge
"""
function normalize_charge_density!(poisson_solver, field_props)
    N = poisson_solver.n_nodes

    @inbounds field_props.charge_density[1] *= 2.0
    @inbounds field_props.charge_density[N] *= 2.0

    @inbounds for j in 1:N
        field_props.charge_density[j] *= poisson_solver.inv_Δx
    end
end


"""
    normalize_charge_density!(poisson_solver::PoissonSolver1DUniform{PeriodicFieldBC1D, PeriodicFieldBC1D},
                              field_props)

Turn the charge deposited on the nodes of a 1D uniform grid by [`deposit_charge!`](@ref) into a charge density,
by dividing it by the volume of the dual cell of each node; periodic boundary conditions.
The charge deposited in the
last node is folded into the first node and mirrored back, and the dual
cell of every node is ``\\Delta x``.

Note that this routine has to be
called **exactly once** per timestep, after all the species have been deposited.

# Positional arguments
* `poisson_solver`: the `PoissonSolver1DUniform` instance with periodic BCs
* `field_props`: the `ElectrostaticFieldProps` instance holding the deposited charge
"""
function normalize_charge_density!(poisson_solver::PoissonSolver1DUniform{PeriodicFieldBC1D, PeriodicFieldBC1D},
                                   field_props)
    N = poisson_solver.n_nodes

    @inbounds field_props.charge_density[1] += field_props.charge_density[N]
    @inbounds field_props.charge_density[N] = field_props.charge_density[1]

    @inbounds for j in 1:N
        field_props.charge_density[j] *= poisson_solver.inv_Δx
    end
end

"""
    deposit_charge!(poisson_solver, grid::Grid1DUniform, particles, pia, species_data, field_props)

Clear the charge density, deposit the charge of the particles of all the species on the nodes of a
1-D uniform grid, and normalize the result to a charge density. This is the recommended way of
computing the charge density, as it performs the clearing, the deposition, and the normalization
in the correct order.

The particles are assumed to be sorted on the grid. This routine is serial.

# Positional arguments
* `poisson_solver`: the `PoissonSolver1DUniform` instance
* `grid`: the `Grid1DUniform` grid
* `particles`: the `Vector` of `ParticleVector`s containing all the particles in a simulation
* `pia`: the `ParticleIndexerArray` instance
* `species_data`: the `Vector` of `Species` data
* `field_props`: the `ElectrostaticFieldProps` instance in which the charge density is stored
"""
function deposit_charge!(poisson_solver::PoissonSolver1DUniform, grid::Grid1DUniform, particles, pia,
                         species_data, field_props)
    clear_charge_density!(field_props)

    for species in 1:pia.n_species
        @inbounds deposit_charge!(grid, particles[species], pia, species, species_data, field_props)
    end

    normalize_charge_density!(poisson_solver, field_props)
end

"""
    reduce_field_props!(field_props_target, field_props_chunks)

Sum the charge deposited in the `ElectrostaticFieldProps` instances of a `field_props_chunks` list
into `field_props_target`, which is cleared first. This is used in a multi-threaded simulation,
where each thread deposits the particles of the cells it owns into its own
`ElectrostaticFieldProps` instance, since depositing into a shared instance is not thread-safe.

Only the charge density is reduced; the potential and the electric field of the target are left
untouched, as they are computed by a single [`solve_poisson!`](@ref) call afterwards. The charge in
the per-thread instances is expected to be the raw deposited charge, so
[`normalize_charge_density!`](@ref) has to be called exactly once, on the target, after the
reduction.

# Positional arguments
* `field_props_target`: the `ElectrostaticFieldProps` instance which will hold the reduced charge density
* `field_props_chunks`: the list of `ElectrostaticFieldProps` instances to use for the reduction operation
"""
function reduce_field_props!(field_props_target, field_props_chunks)
    clear_charge_density!(field_props_target)

    for field_props in field_props_chunks
        @inbounds @simd for j in 1:field_props_target.n_nodes
            field_props_target.charge_density[j] += field_props.charge_density[j]
        end
    end
end

end
