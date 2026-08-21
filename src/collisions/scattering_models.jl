"""
    AbstractScatteringModel

Abstract type for the elastic scattering models used by the DSMC/SWPM collision routines.

Each concrete model is a zero-field singleton type ([`VHS`](@ref), [`VSS`](@ref))
that is used as a compile-time tag: the total cross-section
[`Merzbild.sigma`](@ref) and the computation of the post-collision velocities
[`Merzbild.scatter!`](@ref) are dispatched on it. As the tag is a zero-size value of a concrete
type, it is not stored anywhere and the inner collision loops stay free of branching,
dynamic dispatch, and allocations.

The model of a species pair is stored in the corresponding `Interaction` instance as a
[`Merzbild.ScatteringModel`](@ref) enum value (so that the array of `Interaction` instances stays
concretely typed); the collision drivers convert it to the singleton tag exactly once per call,
see [`Merzbild.@scattering_barrier`](@ref).
"""
abstract type AbstractScatteringModel end

"""
    VHS

Singleton type tag for the Variable Hard Sphere model:
the total cross-section follows the power law ``\\sigma = C g^{1 - 2\\omega}`` and the scattering
is isotropic.
"""
struct VHS <: AbstractScatteringModel end

"""
    VSS

Singleton type tag for the Variable Soft Sphere model: the total cross-section is the same
as for the [`VHS`](@ref) model, but the scattering is anisotropic, with the cosine of the deflection
angle sampled as ``\\cos\\chi = 2 R^{1/\\alpha} - 1``, where ``R`` is a uniformly distributed random
number and ``\\alpha`` is the VSS exponent (``\\alpha = 1`` recovers isotropic scattering).

# References
* K. Koura, H. Matsumoto, Variable soft sphere molecular model for inverse-power-law or Lennard-Jones
    potential. [Phys. Fluids A, 1991](https://doi.org/10.1063/1.857792).
"""
struct VSS <: AbstractScatteringModel end

"""
    ScatteringModel ScatteringVHS=1 ScatteringVSS=2

Enum of the elastic scattering models, stored in an `Interaction` instance to define
the model used for the species pair in question. `ScatteringVHS` corresponds to [`VHS`](@ref),
`ScatteringVSS` to [`VSS`](@ref).
"""
@enum ScatteringModel ScatteringVHS=1 ScatteringVSS=2

"""
Tuple of all `(enum value, singleton tag)` pairs of the implemented elastic scattering models,
used by [`Merzbild.@scattering_barrier`](@ref) to generate the enum-to-tag conversion.
Any newly added scattering model has to be listed here.
"""
const SCATTERING_MODEL_TAGS = ((ScatteringVHS, VHS()), (ScatteringVSS, VSS()))

"""
    @scattering_barrier model call

Turn the [`Merzbild.ScatteringModel`](@ref) enum value `model` into the corresponding singleton
scattering model tag and insert it into the function call `call` as the second positional argument
(that is, directly after `rng`).

The macro expands into an `if`/`elseif` chain with a literal singleton tag in each arm, so that
each arm is a separate call with a concrete tag type. This acts as a function barrier: the branch
is taken once per call of a collision driver, and everything below it is compiled for one fixed
model, with the cross-section and scattering functions inlined and no dispatch left in the
collision loop.

# Example
```julia
@scattering_barrier interaction_l.model ntc!(rng, collision_factors, collision_data, interaction,
                                             particles, pia, cell, species, Δt, V; dw_tol=dw_tol)
```
"""
macro scattering_barrier(model, call)
    Meta.isexpr(call, :call) || throw(ArgumentError("@scattering_barrier expects a function call as its second argument"))

    insert_pos = Meta.isexpr(call.args[2], :parameters) ? 4 : 3

    escaped_model = esc(model)
    expr = :(throw(ArgumentError("unknown scattering model: " * string($escaped_model))))

    for (model_enum, model_tag) in reverse(SCATTERING_MODEL_TAGS)
        branch = copy(call)
        insert!(branch.args, insert_pos, model_tag)
        expr = Expr(:if, :($escaped_model == $model_enum), esc(branch), expr)
    end

    return expr
end

"""
    parse_scattering_model(name)

Convert the name of an elastic scattering model, as written in an interaction data TOML file,
to the corresponding [`Merzbild.ScatteringModel`](@ref) enum value. The comparison is
case-insensitive; the recognized names are `"VHS"` and `"VSS"`.

# Positional arguments
* `name`: the name of the scattering model

# Returns
The `ScatteringModel` enum value.

# Throws
`ArgumentError` if the model name is not recognized.
"""
function parse_scattering_model(name)
    lowercase_name = lowercase(name)

    if lowercase_name == "vhs"
        return ScatteringVHS
    elseif lowercase_name == "vss"
        return ScatteringVSS
    else
        throw(ArgumentError("unknown scattering model: " * name))
    end
end

"""
    sigma(model, interaction, g)

Compute the total elastic collision cross-section for a given scattering `model`.
The [`VSS`](@ref) model uses the same power law as the [`VHS`](@ref) model.

# Positional arguments
* `model`: the `AbstractScatteringModel` singleton tag of the scattering model
* `interaction`: the `Interaction` instance
* `g`: the relative velocity of the collision

# Returns
The value of the computed cross-section.
"""
@inline sigma(::VHS, interaction, g) = sigma_vhs(interaction, g)
@inline sigma(::VSS, interaction, g) = sigma_vhs(interaction, g)
