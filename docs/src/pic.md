# Particle-in-Cell

Merzbild.jl currently implements an electrostatic Particle-in-Cell (PIC) capability on a 1-D uniform grid.
The routines are serial, but the deposition can be threaded over chunks of
cells, see [the section on multithreading](@ref "Multithreading") below; the other
routines are thread-safe for chunks of cells.

Jump to the ["The timestep"](@ref "The timestep") section to see immediately an example
of the function calls for a coupled PIC simulation.

## Node layout

The grid ([`Grid1DUniform`](@ref)) is cell-based, whereas the field quantities are node-centered:

```
node j        1       2       3            n_cells  n_cells+1
              |-------|-------|--- ... ----|---------|
cell i            1       2         ...      n_cells
```

with ``x_j = (j-1) \Delta x``, ``j = 1 \dots n_\textrm{nodes}``, ``n_\textrm{nodes} = n_\textrm{cells} + 1``.
The charge density ``\rho``, the potential ``\phi``, and the x-component of the electric field
``E_x`` are stored in an [`ElectrostaticFieldProps`](@ref) instance as vectors of length
`n_cells + 1`. In a periodic simulation, node
``n_\textrm{nodes}`` is the periodic image of node 1, and the values stored in it are mirrored from
node 1, so that the deposition and the gather are free of any boundary-specific branching — only the
solution of the Poisson equation knows about the boundary conditions.

[`ElectrostaticFieldProps`](@ref) stores no reference to the grid or to the solver, and its
constructor takes the number of **nodes**, in contrast to `PhysProps(n_cells, n_species)`; the
[`ElectrostaticFieldProps(grid::Grid1DUniform)`](@ref) constructor is the recommended way of
creating it.

## Charge deposition

The charge is deposited using first-order (cloud-in-cell) weighting. For a particle of species
``s`` with a computational weight ``w_p`` located in cell ``i`` at position ``x``, with
``\xi = x / \Delta x - (i-1) \in [0,1)``:

```math
Q_i \mathrel{+}= (1 - \xi) q_s w_p, \qquad Q_{i+1} \mathrel{+}= \xi q_s w_p,
```

where ``q_s = q_e Z_s`` is the charge of the species in Coulomb (`Species.charge` stores the charge
in units of the elementary charge). Summing over all the particles gives
``\sum_j Q_j = \sum_p q_s w_p`` exactly.

[`deposit_charge!`](@ref) accumulates the deposited charge of a single species, so the charge
density has to be cleared beforehand with [`clear_charge_density!`](@ref) (which leaves the
potential and the field of the previous timestep untouched) or with [`clear_props!`](@ref).
Once all the species have been deposited, [`normalize_charge_density!`](@ref) turns the deposited
charge into a charge density by dividing it by the volume of the dual cell of each node:

* non-periodic: the dual cell of a boundary node is ``\Delta x / 2``, so
  ``\rho_1 \mathrel{*}= 2``, ``\rho_N \mathrel{*}= 2``, followed by a division of all the nodal
  values by ``\Delta x``;
* periodic: the charge deposited in the last node is folded into the first node and mirrored back
  (``\rho_1 \mathrel{+}= \rho_N``, ``\rho_N = \rho_1``), and every node has a full ``\Delta x``
  dual cell.

The normalization has to be performed exactly once per timestep; the
multi-species convenience method
`deposit_charge!(poisson_solver, grid, particles, pia, species_data, field_props)` performs the
clearing, the deposition of all the species, and the normalization in the correct order and is the
recommended way of computing the charge density. The particles have to be sorted on the grid, as
the cell indices are taken from the `ParticleIndexerArray`. This also means that no particles should be pointed
to by `group2` of a `ParticleIndexer`, and that even if that is the case, merging might still cause issues,
as some merging algorithms displace particles outside of their cell; therefore, particles should be re-sorted
if charge deposition is performed immediately after collisions and merging, even if no convection
occurred.

## The discrete Poisson equation

The Poisson equation is discretized via second-order finite differences:

```math
-\phi'' = \rho / \varepsilon_0 \quad \Rightarrow \quad
\frac{-\phi_{j-1} + 2 \phi_j - \phi_{j+1}}{\Delta x^2} = \frac{\rho_j}{\varepsilon_0}.
```

### Boundary conditions

The field boundary conditions subtype [`Merzbild.AbstractFieldBC1D`](@ref) and are separate from the particle-surface
boundary conditions (which subtype [`Merzbild.AbstractBC`](@ref)):

| Boundary condition | Description | Type |
| --- | --- | --- |
| Dirichlet | Prescribes the potential ``\phi`` at the boundary, in V | [`DirichletFieldBC1D`](@ref) |
| Neumann | Prescribes the x-component of the electric field ``E_x`` at the boundary, in V/m | [`NeumannFieldBC1D`](@ref) |
| Periodic | Periodic domain, has to be used on both sides | [`PeriodicFieldBC1D`](@ref) |

[`DirichletFieldBC1D`](@ref) and [`NeumannFieldBC1D`](@ref) are mutable, so that the prescribed
values can be changed from within a user's time loop, for example to drive an RF electrode
(`bc_left.ϕ = V0 * sin(2π * f * t)`) or to model a dielectric wall accumulating a surface charge
``\sigma`` (`bc_right.E_x = σ / eps_0`). The boundary condition *types* are type parameters of
[`PoissonSolver1DUniform`](@ref) and cannot change.

Note that the Neumann boundary condition prescribes ``E_x`` and not ``d\phi/dx = -E_x``.
Example Neumann BCs: ``E_x = 0`` for a symmetry plane or a floating wall, and
``E_x = \sigma / \varepsilon_0`` for an accumulated surface charge on a dielectric.

Depending on the boundary conditions, the unknowns of the tridiagonal system are the values of the
potential in the following nodes (``N = n_\textrm{nodes}``):

| left BC | right BC | unknown nodes | `n_unknowns` | `node_offset` |
|---|---|---|---|---|
| Dirichlet | Dirichlet | ``2 \dots N-1`` | ``n_\textrm{cells} - 1`` | 1 |
| Neumann | Dirichlet | ``1 \dots N-1`` | ``n_\textrm{cells}`` | 0 |
| Dirichlet | Neumann | ``2 \dots N`` | ``n_\textrm{cells}`` | 1 |
| Neumann | Neumann | — | — | rejected at construction |
| periodic | periodic | ``1 \dots n_\textrm{cells}-1`` | ``n_\textrm{cells} - 1`` | 0 |

with ``\phi_{\textrm{node\_offset} + k}`` being the ``k``-th unknown. Neumann boundary conditions on
both sides are rejected by the constructor of [`PoissonSolver1DUniform`](@ref), as the resulting
system is singular and no gauge is available to fix the constant; a periodic boundary condition
paired with a non-periodic one is rejected as well.

For periodic boundaries, the mean charge density over the ``n_\textrm{cells}`` unique nodes is subtracted,
so that ``\sum_j (\rho_j - \bar{\rho}) = 0``. The subtracted value is
stored in `field_props.net_charge_density` as a diagnostic: a non-zero value means the simulation
carries a net charge, which is unphysical in a periodic domain
In the solver, the gauge is chosen to correspond to a zero mean of the cell-averaged potential.

### The solver

The tridiagonal system is solved by the Thomas algorithm. As the matrix never changes during a
simulation (assuming the PIC grid is constant), the forward elimination of the matrix is performed once in the constructor of
[`PoissonSolver1DUniform`](@ref) (in `factorize_poisson!`) and stored (the diagonals of the matrix are stored as well,
but are not used by the solve).

The matrix is diagonally dominant in every admissible configuration, so the Thomas algorithm is stable
without pivoting.

### The electric field

The electric field ``E_x = -d\phi/dx`` is computed in the nodes by central differences in the
interior,

```math
E_j = -\frac{\phi_{j+1} - \phi_{j-1}}{2 \Delta x}, \qquad j = 2 \dots N-1,
```

whereas in the boundary nodes a one-sided difference is used for a Dirichlet boundary condition, the
prescribed value is copied for a Neumann boundary condition, and the stencil wraps around the domain
for a periodic one.

[`solve_poisson!`](@ref) performs the whole sequence
after the charge density has been computed: the assembly of the right-hand side, the
tridiagonal solve, the write-back of the known values (and, in the periodic case, of the pinned node
and the gauge), and the computation of the electric field. It does not modify the computed charge density.

## Gather and push

The electric field is interpolated to the position of a particle with the same shape function as
the one used for the deposition:

```math
E_p = (1 - \xi) E_i + \xi E_{i+1}, \qquad v_x \mathrel{+}= \frac{q_s}{m_s} E_p \Delta t.
```
The push is performed by
[`accelerate_electric_field_x!`](@ref), either for a single cell or for all the cells of a species.
Note that it needs neither the solver nor the charge density — only the grid (for ``1/\Delta x``),
the particle indexing, and the field.

To use leapfrog time-stepping, one should initialize the simulation
by a single backward half-kick after the first solution of the Poisson equation:

```julia
deposit_charge!(poisson_solver, grid, particles, pia, species_data, field_props)
solve_poisson!(poisson_solver, field_props)
accelerate_electric_field_x!(grid, particles[1], pia, 1, species_data, field_props, -0.5 * Δt)
```

## The timestep

A PIC timestep is (assuming a leapfrog push with an already existing time-offset, see above)

```julia
deposit_charge!(poisson_solver, grid, particles, pia, species_data, field_props)
solve_poisson!(poisson_solver, field_props)

for species in 1:pia.n_species
    accelerate_electric_field_x!(grid, particles[species], pia, species, species_data, field_props, Δt)
    convect_particles!(rng, grid, bc_list, particles[species], pia, species, species_data, Δt)
    sort_particles!(gridsorter, grid, particles[species], pia, species)
end

# collisions, merging, computation of the macroscopic properties, I/O
```

with `convect_particles_periodic!` used instead of `convect_particles!` in a periodic domain.
The deposition has to run on sorted particles, so the sorting closes the timestep.

**Important note**: if merging is performed, depending on the algorithm used, particles may end up
outside of their pre-merge cells: one needs to re-sort the particles if that is the case.

## Multithreading

The deposition writes to nodes that are shared between neighbouring cells, so depositing the cells
owned by different threads into a **shared** [`ElectrostaticFieldProps`](@ref) instance is not
thread-safe: the two cells meeting at a node both update it.

Each thread therefore deposits the cells it owns into its own [`ElectrostaticFieldProps`](@ref)
instance, using the `cell_chunk` version of [`deposit_charge!`](@ref); the per-thread instances are
then summed into the global one with [`reduce_field_props!`](@ref), which is followed by a single
[`normalize_charge_density!`](@ref) call (the per-thread instances hold the raw deposited charge)
and a single [`solve_poisson!`](@ref):

```julia
Threads.@threads for chunk_id in 1:n_chunks
    clear_charge_density!(field_props_chunks[chunk_id])
    for species in 1:pia.n_species
        deposit_charge!(grid, particles_chunks[chunk_id][species], pia_chunks[chunk_id], species,
                        species_data, field_props_chunks[chunk_id], cell_chunks[chunk_id])
    end
end

reduce_field_props!(field_props, field_props_chunks)
normalize_charge_density!(poisson_solver, field_props)
solve_poisson!(poisson_solver, field_props)
```

[`reduce_field_props!`](@ref) reduces only the charge density and leaves the potential and the
electric field of the target untouched, as those are computed by the solve that follows.

The Poisson solve itself is inherently serial, but is ``O(n_\textrm{cells})`` against the
``O(n_\textrm{particles})`` deposition and push, so leaving it serial is not a bottleneck.
The gather and push ([`accelerate_electric_field_x!`](@ref)) only read the electric field and write
per-particle velocities, so they are safe to run per cell chunk on the shared instance.

## Stability constraints

An electrostatic PIC simulation resolves the electron plasma oscillation and the Debye shielding,
which imposes

```math
\Delta x \lesssim \lambda_D, \qquad \omega_p \Delta t \lesssim 0.2,
```

where the Debye length and the plasma frequency are computed by [`debye_length`](@ref) and
[`plasma_frequency`](@ref). Violating the first constraint leads to the finite-grid instability
(a numerical heating of the plasma), violating the second one to an unstable particle push.

## I/O

Output of the electrostatic field data to NetCDF format is done via [`ElectrostaticFieldProps`](@ref)
and [`write_netcdf`](@ref), with optional omitting of output fields via [`IOSkipListField`](@ref):

```julia
field_props = ElectrostaticFieldProps(grid)
ds_field = NCDataHolderField("plasma_fields.nc", field_props)

for t in 1:n_timesteps
  # compute fields, move particles, etc.

  write_netcdf(ds_field, field_props, t)
end

close_netcdf(ds_field)
```

The field data in the output is independent of the underlying grid and contains only the data
on the electric field, potential, and charge density. The values of the boundary conditions
are not written out but support will be added in the future.

## Examples

Two example simulations are provided in `simulations/1D`: `plasma_oscillation.jl`, a cold electron
plasma oscillation in a periodic domain on top of a fixed neutralizing ion background, and
`rf_electrode.jl`, a plasma between an RF-driven electrode (a Dirichlet boundary condition updated
from within the time loop) and a floating wall (a Neumann boundary condition with ``E_x = 0``).
