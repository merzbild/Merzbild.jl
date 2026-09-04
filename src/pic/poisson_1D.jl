@muladd begin

"""
    PoissonSolver1DUniform

Structure holding data used to discretize and solve
the 1-D Poisson equation ``-\\phi'' = \\rho / \\varepsilon_0`` via second-order finite differences
using the Thomas algorithm for tridiagonal matrices.

The field quantities themselves are not stored here, but in an [`ElectrostaticFieldProps`](@ref)
instance. The boundary conditions are type parameters, so that the assembly of the right-hand side
is resolved at compile time.

The unknowns are the values of the potential in the nodes which are not prescribed by the boundary
conditions; the potential in node `node_offset + k` is the `k`-th unknown. The number of unknowns
and the node offset depend on the boundary conditions (`N = n_nodes`, `n_cells = N - 1`):

| left BC | right BC | unknown nodes | `n_unknowns` | `node_offset` |
|---|---|---|---|---|
| Dirichlet | Dirichlet | `2 … N-1` | `n_cells - 1` | 1 |
| Neumann | Dirichlet | `1 … N-1` | `n_cells` | 0 |
| Dirichlet | Neumann | `2 … N` | `n_cells` | 1 |
| Neumann | Neumann | — | — | rejected at construction |
| periodic | periodic | `1 … n_cells-1` | `n_cells - 1` | 0 |

# Fields
* `n_nodes`: number of nodes
* `n_unknowns`: number of unknowns of the tridiagonal system
* `node_offset`: offset between the index of an unknown and the index of the corresponding node
* `Δx`: cell size
* `inv_Δx`: inverse of the cell size
* `inv_Δx2`: inverse of the squared cell size, used for the Dirichlet right-hand side correction
* `inv_eps_0`: inverse of the vacuum permittivity
* `bc_left`: the `AbstractFieldBC1D` boundary condition on the left boundary
* `bc_right`: the `AbstractFieldBC1D` boundary condition on the right boundary
* `a`: sub-diagonal of the matrix, with `a[1] = 0` as it lies outside of the matrix
    (not used in the solve)
* `b`: main diagonal of the matrix (not used in the solve)
* `c`: super-diagonal of the matrix, with `c[n_unknowns] = 0` as it lies outside of the matrix
    (not used in the solve)
* `cp`: pre-computed `c[i] / m[i]` factorization coefficients
* `inv_m`: pre-computed `1 / (b[i] - cp[i-1] * a[i])` factorization coefficients
* `a_inv_m`: pre-computed `-a[i] / m[i]` factorization coefficients
* `rhs`: right-hand side of the tridiagonal system
* `dp`: workspace of the forward sweep of the Thomas algorithm
* `x`: solution of the tridiagonal system (values of the potential in the unknown nodes)
"""
mutable struct PoissonSolver1DUniform{BCL<:AbstractFieldBC1D, BCR<:AbstractFieldBC1D}
    n_nodes::Int64
    n_unknowns::Int64
    node_offset::Int64
    Δx::Float64
    inv_Δx::Float64
    inv_Δx2::Float64
    inv_eps_0::Float64
    bc_left::BCL
    bc_right::BCR
    a::Vector{Float64}
    b::Vector{Float64}
    c::Vector{Float64}
    cp::Vector{Float64}
    inv_m::Vector{Float64}
    a_inv_m::Vector{Float64}
    rhs::Vector{Float64}
    dp::Vector{Float64}
    x::Vector{Float64}
end

"""
    poisson_unknowns(n_cells, bc_left, bc_right)

Compute the number of unknowns of the discrete 1-D Poisson system and the offset between
the index of an unknown and the index of the corresponding node, for a given pair of boundary
conditions. Throws an error for the boundary condition combinations which are not admissible
(Neumann on both sides, which is singular and leaves no way of fixing the gauge, and a periodic
boundary condition paired with a non-periodic one).

# Positional arguments
* `n_cells`: number of cells of the grid
* `bc_left`: the `AbstractFieldBC1D` boundary condition on the left boundary
* `bc_right`: the `AbstractFieldBC1D` boundary condition on the right boundary

# Returns
Tuple of the number of unknowns and the node offset.
"""
poisson_unknowns(n_cells, bc_left::DirichletFieldBC1D, bc_right::DirichletFieldBC1D) = (n_cells - 1, Int64(1))
poisson_unknowns(n_cells, bc_left::NeumannFieldBC1D, bc_right::DirichletFieldBC1D) = (n_cells, Int64(0))
poisson_unknowns(n_cells, bc_left::DirichletFieldBC1D, bc_right::NeumannFieldBC1D) = (n_cells, Int64(1))
poisson_unknowns(n_cells, bc_left::PeriodicFieldBC1D, bc_right::PeriodicFieldBC1D) = (n_cells - 1, Int64(0))

function poisson_unknowns(n_cells, bc_left::NeumannFieldBC1D, bc_right::NeumannFieldBC1D)
    throw(ErrorException("Neumann boundary conditions on both sides of the domain lead to a singular"
                         * " Poisson system with no way of fixing the gauge;"
                         * " prescribe the potential on at least one boundary"))
end

function poisson_unknowns(n_cells, bc_left::AbstractFieldBC1D, bc_right::AbstractFieldBC1D)
    throw(ErrorException("A periodic field boundary condition has to be used on both sides"
                         * " of the domain, got $(typeof(bc_left)) and $(typeof(bc_right))"))
end

"""
    set_neumann_diagonal_left!(poisson_solver, bc_left)

Replace the first row of the discrete Poisson matrix by the finite-volume integration over the
left half-cell if the left boundary condition is a Neumann one, otherwise do nothing.
The diagonal entry is halved with respect to an interior row, as the control volume of the
boundary node is ``\\Delta x / 2``.

# Positional arguments
* `poisson_solver`: the `PoissonSolver1DUniform` instance
* `bc_left`: the `AbstractFieldBC1D` boundary condition on the left boundary
"""
@inline set_neumann_diagonal_left!(poisson_solver, bc_left::AbstractFieldBC1D) = nothing

@inline function set_neumann_diagonal_left!(poisson_solver, bc_left::NeumannFieldBC1D)
    poisson_solver.b[1] = poisson_solver.inv_Δx2
    return nothing
end

"""
    set_neumann_diagonal_right!(poisson_solver, bc_right)

Replace the last row of the discrete Poisson matrix by the finite-volume integration over the
right half-cell if the right boundary condition is a Neumann one, otherwise do nothing.
The diagonal entry is halved with respect to an interior row, as the control volume of the
boundary node is ``\\Delta x / 2``.

# Positional arguments
* `poisson_solver`: the `PoissonSolver1DUniform` instance
* `bc_right`: the `AbstractFieldBC1D` boundary condition on the right boundary
"""
@inline set_neumann_diagonal_right!(poisson_solver, bc_right::AbstractFieldBC1D) = nothing

@inline function set_neumann_diagonal_right!(poisson_solver, bc_right::NeumannFieldBC1D)
    poisson_solver.b[poisson_solver.n_unknowns] = poisson_solver.inv_Δx2
    return nothing
end

"""
    PoissonSolver1DUniform(grid::Grid1DUniform, bc_left::AbstractFieldBC1D, bc_right::AbstractFieldBC1D)

Construct a solver of the 1-D Poisson equation on a uniform grid for a given pair of boundary
conditions and factorize the resulting tridiagonal matrix.

Neumann boundary conditions on both sides are rejected (the system is singular and no gauge is
available), as is a periodic boundary condition paired with a non-periodic one.

# Positional arguments
* `grid`: the `Grid1DUniform` grid
* `bc_left`: the `AbstractFieldBC1D` boundary condition on the left boundary
* `bc_right`: the `AbstractFieldBC1D` boundary condition on the right boundary
"""
function PoissonSolver1DUniform(grid::Grid1DUniform, bc_left::AbstractFieldBC1D, bc_right::AbstractFieldBC1D)
    n_unknowns, node_offset = poisson_unknowns(grid.n_cells, bc_left, bc_right)

    if n_unknowns < 1
        throw(ErrorException("A grid with $(grid.n_cells) cell(s) is too coarse for the chosen field"
                             * " boundary conditions, the Poisson system has $(n_unknowns) unknowns"))
    end

    inv_Δx2 = grid.inv_Δx * grid.inv_Δx

    a = fill(-inv_Δx2, n_unknowns)
    b = fill(2.0 * inv_Δx2, n_unknowns)
    c = fill(-inv_Δx2, n_unknowns)
    a[1] = 0.0
    c[n_unknowns] = 0.0

    poisson_solver = PoissonSolver1DUniform(grid.n_cells + 1, n_unknowns, node_offset,
                                            grid.Δx, grid.inv_Δx, inv_Δx2, 1.0 / eps_0,
                                            bc_left, bc_right,
                                            a, b, c,
                                            zeros(n_unknowns), zeros(n_unknowns), zeros(n_unknowns),
                                            zeros(n_unknowns), zeros(n_unknowns), zeros(n_unknowns))

    set_neumann_diagonal_left!(poisson_solver, bc_left)
    set_neumann_diagonal_right!(poisson_solver, bc_right)
    factorize_poisson!(poisson_solver)

    return poisson_solver
end

"""
    factorize_poisson!(poisson_solver)

Pre-compute the forward elimination of the tridiagonal matrix of the Thomas algorithm and store
the `inv_m`, `cp`, and `a_inv_m` coefficients. These depend only on the matrix and not on the
right-hand side, so with a fixed matrix only the right-hand side recurrence and the back
substitution have to be performed per solve.

This assumes that the matrix never changes during a simulations, i.e. PIC grid and BC **types**
are fixed (the prescribed values of the potential/electric fields at Dirichlet/Neumann BCs can vary).

# Positional arguments
* `poisson_solver`: the `PoissonSolver1DUniform` instance to be factorized
"""
function factorize_poisson!(poisson_solver)
    n = poisson_solver.n_unknowns

    @inbounds poisson_solver.inv_m[1] = 1.0 / poisson_solver.b[1]
    @inbounds poisson_solver.cp[1] = poisson_solver.c[1] * poisson_solver.inv_m[1]
    @inbounds poisson_solver.a_inv_m[1] = 0.0

    @inbounds for i in 2:n
        poisson_solver.inv_m[i] = 1.0 / (poisson_solver.b[i] - poisson_solver.cp[i-1] * poisson_solver.a[i])
        poisson_solver.cp[i] = poisson_solver.c[i] * poisson_solver.inv_m[i]
        poisson_solver.a_inv_m[i] = -poisson_solver.a[i] * poisson_solver.inv_m[i]
    end
end

"""
    poisson_rhs_left!(poisson_solver, field_props, bc_left)

Apply the correction of the first row of the right-hand side of the discrete Poisson system
stemming from the left boundary condition: for a Dirichlet boundary condition, the eliminated
known value of the potential in the boundary node is added to the right-hand side; for a Neumann
boundary condition, the row is that of the finite-volume integration over the left half-cell,
with a halved charge density and the prescribed value of the electric field.

# Positional arguments
* `poisson_solver`: the `PoissonSolver1DUniform` instance
* `field_props`: the `ElectrostaticFieldProps` instance holding the charge density
* `bc_left`: the `AbstractFieldBC1D` boundary condition on the left boundary
"""
@inline function poisson_rhs_left!(poisson_solver, field_props, bc_left::DirichletFieldBC1D)
    @inbounds poisson_solver.rhs[1] += bc_left.ϕ * poisson_solver.inv_Δx2
    return nothing
end

@inline function poisson_rhs_left!(poisson_solver, field_props, bc_left::NeumannFieldBC1D)
    @inbounds poisson_solver.rhs[1] = 0.5 * field_props.charge_density[1] * poisson_solver.inv_eps_0 +
                                      bc_left.E_x * poisson_solver.inv_Δx
    return nothing
end

"""
    poisson_rhs_right!(poisson_solver, field_props, bc_right)

Apply the correction of the last row of the right-hand side of the discrete Poisson system
stemming from the right boundary condition, mirroring [`poisson_rhs_left!`](@ref).

# Positional arguments
* `poisson_solver`: the `PoissonSolver1DUniform` instance
* `field_props`: the `ElectrostaticFieldProps` instance holding the charge density
* `bc_right`: the `AbstractFieldBC1D` boundary condition on the right boundary
"""
@inline function poisson_rhs_right!(poisson_solver, field_props, bc_right::DirichletFieldBC1D)
    @inbounds poisson_solver.rhs[poisson_solver.n_unknowns] += bc_right.ϕ * poisson_solver.inv_Δx2
    return nothing
end

@inline function poisson_rhs_right!(poisson_solver, field_props, bc_right::NeumannFieldBC1D)
    @inbounds poisson_solver.rhs[poisson_solver.n_unknowns] =
        0.5 * field_props.charge_density[poisson_solver.n_nodes] * poisson_solver.inv_eps_0 -
        bc_right.E_x * poisson_solver.inv_Δx
    return nothing
end

"""
    assemble_poisson_rhs!(poisson_solver, field_props)

Assemble the right-hand side of the discrete Poisson system from the charge density stored in
`field_props`, applying the boundary corrections. In the periodic case, the mean charge density
over the unique nodes is subtracted (otherwise the singular periodic system has no solution at all)
and stored in `field_props.net_charge_density`; a non-zero value of the latter means that the
simulation carries a net charge, which is unphysical in a periodic domain.

The charge density is not modified.

# Positional arguments
* `poisson_solver`: the `PoissonSolver1DUniform` instance
* `field_props`: the `ElectrostaticFieldProps` instance holding the charge density
"""
function assemble_poisson_rhs!(poisson_solver::PoissonSolver1DUniform, field_props)
    @inbounds for k in 1:poisson_solver.n_unknowns
        poisson_solver.rhs[k] = field_props.charge_density[poisson_solver.node_offset + k] * poisson_solver.inv_eps_0
    end

    poisson_rhs_left!(poisson_solver, field_props, poisson_solver.bc_left)
    poisson_rhs_right!(poisson_solver, field_props, poisson_solver.bc_right)

    field_props.net_charge_density = 0.0
end

function assemble_poisson_rhs!(poisson_solver::PoissonSolver1DUniform{PeriodicFieldBC1D, PeriodicFieldBC1D},
                               field_props)
    n_cells = poisson_solver.n_nodes - 1

    ρ_mean = 0.0
    @inbounds for j in 1:n_cells
        ρ_mean += field_props.charge_density[j]
    end
    ρ_mean /= n_cells

    @inbounds for k in 1:poisson_solver.n_unknowns
        poisson_solver.rhs[k] = (field_props.charge_density[k] - ρ_mean) * poisson_solver.inv_eps_0
    end

    field_props.net_charge_density = ρ_mean
end

"""
    solve_tridiagonal!(poisson_solver)

Solve the tridiagonal system using the pre-computed factorization of the matrix
(see [`factorize_poisson!`](@ref)) and the right-hand side stored in `poisson_solver.rhs`,
storing the result in `poisson_solver.x`.

# Positional arguments
* `poisson_solver`: the `PoissonSolver1DUniform` instance
"""
function solve_tridiagonal!(poisson_solver)
    n = poisson_solver.n_unknowns

    @inbounds poisson_solver.dp[1] = poisson_solver.rhs[1] * poisson_solver.inv_m[1]
    @inbounds for i in 2:n
        poisson_solver.dp[i] = poisson_solver.rhs[i] * poisson_solver.inv_m[i] +
                               poisson_solver.dp[i-1] * poisson_solver.a_inv_m[i]
    end

    @inbounds poisson_solver.x[n] = poisson_solver.dp[n]
    @inbounds for i in n-1:-1:1
        poisson_solver.x[i] = poisson_solver.dp[i] - poisson_solver.cp[i] * poisson_solver.x[i+1]
    end
end

"""
    apply_gauge!(poisson_solver, field_props)

Shift the potential in the unique nodes `1:n_cells` of a periodic domain so that its mean is zero.
On a periodic uniform grid this gauge is identical to requiring a zero mean of the cell-averaged
potential, as each unique node contributes to exactly two cells.

# Positional arguments
* `poisson_solver`: the `PoissonSolver1DUniform` instance
* `field_props`: the `ElectrostaticFieldProps` instance holding the potential
"""
function apply_gauge!(poisson_solver, field_props)
    n_cells = poisson_solver.n_nodes - 1

    ϕ_mean = 0.0
    @inbounds for j in 1:n_cells
        ϕ_mean += field_props.potential[j]
    end
    ϕ_mean /= n_cells

    @inbounds for j in 1:n_cells
        field_props.potential[j] -= ϕ_mean
    end
end

"""
    finalize_potential_left!(field_props, bc_left)

Write the known value of the potential into the left boundary node for a Dirichlet boundary
condition; for a Neumann boundary condition the boundary node is an unknown of the system and
nothing is done.

# Positional arguments
* `field_props`: the `ElectrostaticFieldProps` instance holding the potential
* `bc_left`: the `AbstractFieldBC1D` boundary condition on the left boundary
"""
@inline finalize_potential_left!(field_props, bc_left::AbstractFieldBC1D) = nothing

@inline function finalize_potential_left!(field_props, bc_left::DirichletFieldBC1D)
    @inbounds field_props.potential[1] = bc_left.ϕ
    return nothing
end

"""
    finalize_potential_right!(field_props, bc_right)

Write the known value of the potential into the right boundary node for a Dirichlet boundary
condition, mirroring [`finalize_potential_left!`](@ref).

# Positional arguments
* `field_props`: the `ElectrostaticFieldProps` instance holding the potential
* `bc_right`: the `AbstractFieldBC1D` boundary condition on the right boundary
"""
@inline finalize_potential_right!(field_props, bc_right::AbstractFieldBC1D) = nothing

@inline function finalize_potential_right!(field_props, bc_right::DirichletFieldBC1D)
    @inbounds field_props.potential[field_props.n_nodes] = bc_right.ϕ
    return nothing
end

"""
    finalize_potential!(poisson_solver, field_props)

Write the values of the potential which are not part of the tridiagonal solution into the nodes:
the prescribed values for Dirichlet boundary conditions, and, in the periodic case, the pinned
value ``\\phi = 0`` in node `n_cells` (which has to be written every timestep, as a stale value
left over from the previous timestep would corrupt both the gauge and the electric field in the
neighbourhood of that node), followed by the gauge shift and the mirroring of node 1 into the
last node.

# Positional arguments
* `poisson_solver`: the `PoissonSolver1DUniform` instance
* `field_props`: the `ElectrostaticFieldProps` instance holding the potential
"""
function finalize_potential!(poisson_solver::PoissonSolver1DUniform, field_props)
    finalize_potential_left!(field_props, poisson_solver.bc_left)
    finalize_potential_right!(field_props, poisson_solver.bc_right)
end

function finalize_potential!(poisson_solver::PoissonSolver1DUniform{PeriodicFieldBC1D, PeriodicFieldBC1D},
                             field_props)
    n_cells = poisson_solver.n_nodes - 1
    @inbounds field_props.potential[n_cells] = 0.0
    apply_gauge!(poisson_solver, field_props)
    @inbounds field_props.potential[poisson_solver.n_nodes] = field_props.potential[1]
end

"""
    electric_field_left!(poisson_solver, field_props, bc_left)

Compute the x-component of the electric field in the left boundary node: a one-sided difference
of the potential for a Dirichlet boundary condition, and the prescribed value for a Neumann
boundary condition.

# Positional arguments
* `poisson_solver`: the `PoissonSolver1DUniform` instance
* `field_props`: the `ElectrostaticFieldProps` instance holding the potential and the electric field
* `bc_left`: the `AbstractFieldBC1D` boundary condition on the left boundary
"""
@inline function electric_field_left!(poisson_solver, field_props, bc_left::DirichletFieldBC1D)
    @inbounds field_props.electric_field[1] = -(field_props.potential[2] - field_props.potential[1]) * poisson_solver.inv_Δx
    return nothing
end

@inline function electric_field_left!(poisson_solver, field_props, bc_left::NeumannFieldBC1D)
    @inbounds field_props.electric_field[1] = bc_left.E_x
    return nothing
end

"""
    electric_field_right!(poisson_solver, field_props, bc_right)

Compute the x-component of the electric field in the right boundary node: a one-sided difference
of the potential for a Dirichlet boundary condition, and the prescribed value for a Neumann
boundary condition.

# Positional arguments
* `poisson_solver`: the `PoissonSolver1DUniform` instance
* `field_props`: the `ElectrostaticFieldProps` instance holding the potential and the electric field
* `bc_right`: the `AbstractFieldBC1D` boundary condition on the right boundary
"""
@inline function electric_field_right!(poisson_solver, field_props, bc_right::DirichletFieldBC1D)
    N = poisson_solver.n_nodes
    @inbounds field_props.electric_field[N] = -(field_props.potential[N] - field_props.potential[N-1]) * poisson_solver.inv_Δx
    return nothing
end

@inline function electric_field_right!(poisson_solver, field_props, bc_right::NeumannFieldBC1D)
    @inbounds field_props.electric_field[poisson_solver.n_nodes] = bc_right.E_x
    return nothing
end

"""
    compute_electric_field!(poisson_solver, field_props)

Compute the x-component of the electric field ``E_x = -d\\phi/dx`` in the nodes from the potential,
using central differences in the interior nodes and the boundary condition-specific treatment of
the boundary nodes.

# Positional arguments
* `poisson_solver`: the `PoissonSolver1DUniform` instance
* `field_props`: the `ElectrostaticFieldProps` instance holding the potential and the electric field
"""
function compute_electric_field!(poisson_solver::PoissonSolver1DUniform, field_props)
    half_inv_Δx = 0.5 * poisson_solver.inv_Δx

    @inbounds for j in 2:poisson_solver.n_nodes-1
        field_props.electric_field[j] = -(field_props.potential[j+1] - field_props.potential[j-1]) * half_inv_Δx
    end

    electric_field_left!(poisson_solver, field_props, poisson_solver.bc_left)
    electric_field_right!(poisson_solver, field_props, poisson_solver.bc_right)
end

function compute_electric_field!(poisson_solver::PoissonSolver1DUniform{PeriodicFieldBC1D, PeriodicFieldBC1D},
                                 field_props)
    N = poisson_solver.n_nodes
    n_cells = N - 1
    half_inv_Δx = 0.5 * poisson_solver.inv_Δx

    @inbounds for j in 2:N-1
        field_props.electric_field[j] = -(field_props.potential[j+1] - field_props.potential[j-1]) * half_inv_Δx
    end

    @inbounds field_props.electric_field[1] = -(field_props.potential[2] - field_props.potential[n_cells]) * half_inv_Δx
    @inbounds field_props.electric_field[N] = field_props.electric_field[1]
end

"""
    solve_poisson!(poisson_solver, field_props)

Compute the electrostatic potential and the electric field in the nodes from the charge density
stored in `field_props`, by assembling and solving the discrete Poisson system.
The charge density is not modified, so repeated calls do not change the result; it is however assumed
that the charge density has already been normalized by a call to
[`normalize_charge_density!`](@ref).

# Positional arguments
* `poisson_solver`: the `PoissonSolver1DUniform` instance
* `field_props`: the `ElectrostaticFieldProps` instance holding the charge density, the potential,
    and the electric field
"""
function solve_poisson!(poisson_solver, field_props)
    assemble_poisson_rhs!(poisson_solver, field_props)
    solve_tridiagonal!(poisson_solver)

    @inbounds for k in 1:poisson_solver.n_unknowns
        field_props.potential[poisson_solver.node_offset + k] = poisson_solver.x[k]
    end

    finalize_potential!(poisson_solver, field_props)
    compute_electric_field!(poisson_solver, field_props)
end

end
