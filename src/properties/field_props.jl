@muladd begin

"""
    ElectrostaticFieldProps

Structure to store the electrostatic field quantities on the nodes of a grid.
The quantities are node-centered: for a 1D uniform grid with `n_cells` cells, the nodes are located at
``x_j = (j-1) \\Delta x``, ``j = 1 \\dots n_{\\textrm{cells}}+1``, so all the vectors are of
length `n_cells + 1`. For a periodic simulation, the last node is the periodic image of the first
one, and the values stored in it are mirrored from the first node.

The structure is grid-independent: it stores no reference to a grid or to a Poisson solver.

# Fields
* `n_nodes`: number of nodes
* `charge_density`: the charge density ``\\rho`` in the nodes, C/m³ (before a call to
    [`normalize_charge_density!`](@ref) this holds the deposited charge in C and not a density)
* `potential`: the electrostatic potential ``\\phi`` in the nodes, V
* `electric_field`: the x-component of the electric field ``E_x`` in the nodes, V/m
* `net_charge_density`: the mean charge density subtracted in the periodic case, `0.0` otherwise, C/m³
"""
mutable struct ElectrostaticFieldProps
    n_nodes::Int64
    charge_density::Vector{Float64}
    potential::Vector{Float64}
    electric_field::Vector{Float64}
    net_charge_density::Float64
end

"""
    ElectrostaticFieldProps(n_nodes::Integer)

Construct an `ElectrostaticFieldProps` instance given the number of **nodes**
(and not the number of cells, in contrast to `PhysProps(n_cells, n_species)`).
For a 1-D uniform grid with `n_cells` cells, the number of nodes is `n_cells + 1`.

# Positional arguments
* `n_nodes`: number of nodes
"""
ElectrostaticFieldProps(n_nodes::Integer) = ElectrostaticFieldProps(n_nodes, zeros(n_nodes), zeros(n_nodes),
                                                                   zeros(n_nodes), 0.0)

"""
    ElectrostaticFieldProps(grid::Grid1DUniform)

Construct an `ElectrostaticFieldProps` instance for a 1-D uniform grid, using
`n_nodes = grid.n_cells + 1` nodes.

# Positional arguments
* `grid`: the `Grid1DUniform` grid
"""
ElectrostaticFieldProps(grid::Grid1DUniform) = ElectrostaticFieldProps(grid.n_cells + 1)

"""
    clear_props!(field_props::ElectrostaticFieldProps)

Clear all data from an `ElectrostaticFieldProps` instance (the charge density, the potential,
and the electric field). For clearing only the charge density in a time loop,
use [`clear_charge_density!`](@ref).

# Positional arguments
* `field_props`: the `ElectrostaticFieldProps` instance to be cleared
"""
function clear_props!(field_props::ElectrostaticFieldProps)
    fill!(field_props.charge_density, 0.0)
    fill!(field_props.potential, 0.0)
    fill!(field_props.electric_field, 0.0)
    field_props.net_charge_density = 0.0
    return nothing
end

"""
    clear_charge_density!(field_props::ElectrostaticFieldProps)

Clear the charge density stored in an `ElectrostaticFieldProps` instance, leaving the potential
and the electric field computed at the previous timestep untouched. This is to be called before
the charge of the particles is deposited on the nodes at a new timestep.

# Positional arguments
* `field_props`: the `ElectrostaticFieldProps` instance for which the charge density is cleared
"""
function clear_charge_density!(field_props::ElectrostaticFieldProps)
    fill!(field_props.charge_density, 0.0)
    field_props.net_charge_density = 0.0
    return nothing
end

end
