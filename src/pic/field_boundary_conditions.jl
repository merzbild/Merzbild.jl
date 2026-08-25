@muladd begin

"""
    DirichletFieldBC1D <: AbstractFieldBC1D

A boundary condition prescribing the value of the electrostatic potential at a boundary of a 1-D domain.
The structure is mutable, so that the value can be changed from within a time loop
(for example, to drive an RF electrode: `bc_left.ϕ = V0 * sin(2π * f * t)`).

# Fields
* `ϕ`: the prescribed value of the potential, V
"""
mutable struct DirichletFieldBC1D <: AbstractFieldBC1D
    ϕ::Float64
end

"""
    NeumannFieldBC1D <: AbstractFieldBC1D

A boundary condition prescribing the x-component of the electric field at a boundary of a 1-D domain
(and not the derivative of the potential, given by ``d\\phi/dx = -E_x``).
The structure is mutable, so that the value can be changed from within a time loop
(for example, for a dielectric wall accumulating a surface charge ``\\sigma``:
`bc_right.E_x = σ / eps_0`).

# Fields
* `E_x`: the prescribed value of the x-component of the electric field, V/m
"""
mutable struct NeumannFieldBC1D <: AbstractFieldBC1D
    E_x::Float64
end

"""
    PeriodicFieldBC1D <: AbstractFieldBC1D

A periodic boundary condition for the electrostatic field in a 1-D domain. It has to be used
on both sides of the domain; as no values are prescribed, no data is stored in the struct.
"""
struct PeriodicFieldBC1D <: AbstractFieldBC1D
end

end
