
"""
    AbstractGrid

Abstract grid type
"""
abstract type AbstractGrid end

"""
    AbstractNCDataHolder

Abstract type that holds NetCDF-output related data for I/O
"""
abstract type AbstractNCDataHolder end

"""
    AbstractBC

Abstract type for boundary conditions
"""
abstract type AbstractBC end

"""
    AbstractFieldBC1D

Abstract type for boundary conditions imposed on the electrostatic field in 1-D simulations.
Not to be confused with [`AbstractBC`](@ref), which describes particle-surface interaction.
"""
abstract type AbstractFieldBC1D end