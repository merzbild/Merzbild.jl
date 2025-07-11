"""
    Cell1D

Cell element of a 1-D grid

# Fields
* `xlo`: coordinate of the left end of the element
* `xhi`: coordinate of the right end of the element 
* `V`: cell volume
* `inv_V`: inverse of cell volume
"""
struct Cell1D
    xlo::Float64
    xhi::Float64
    V::Float64
    inv_V::Float64

    @doc """
        Cell1D(xlo, xhi, V)

    Create cell of 1-D uniform grid

    # Positional arguments
    * `xlo`: coordinate of the left end of the element
    * `xhi`: coordinate of the right end of the element 
    * `V`: cell volume
    """
    function Cell1D(xlo, xhi, V)
        return new(xlo, xhi, V, 1.0/V)
    end
end

"""
    Grid1DUniform

1-D Uniform grid for a domain ``[0,L]``

# Fields
* `L`: length of the domain
* `nx`: number of cells
* `Δx`: cell size
* `inv_Δx`: inverse of cell size
* `cells`: `Vector` of `Cell1D` elements
* `min_x`: minimum allowed `x` coordinate for particles (slightly larger than ``0``)
* `max_x`: maximum allowed `x` coordinate for particles (slightly smaller than ``L``)
"""
struct Grid1DUniform
    L::Float64
    n_cells::Int64
    Δx::Float64
    inv_Δx::Float64
    cells::Vector{Cell1D}
    min_x::Float64  # so that we don't get particles stuck exactly at the wall
    max_x::Float64  # so that we don't get particles stuck exactly at the wall

    @doc """
        Grid1DUniform(L, nx; wall_offset=1e-12)

    Create 1-D uniform grid for a domain ``[0, L]`` with `nx` cells

    # Positional arguments
    * `L`: length of the domain
    * `nx`: number of cells

    # Keyword arguments
    * `wall_offset`: specifies a small relative offset from the walls so that particles
        never end up with a coordinate of exactly ``0`` or ``L``, otherwise some
        sorting routines may produce cell indices outside of the `1:nx` range.
        The offset is computed as `Δx * wall_offset`, where `Δx` is the cell size.
    """
    function Grid1DUniform(L, nx; wall_offset=1e-12)
        cells = Vector{Cell1D}(undef, nx)
        dx = L / nx

        for i in 1:nx
            xlo = (i-1) * dx
            xhi = i * dx
            V = dx
            cells[i] = Cell1D(xlo, xhi, V)
        end

        return new(L, nx, dx, 1.0 / dx, cells, dx * wall_offset, L - dx * wall_offset)
    end
end

"""
    get_cell(grid1duniform::Grid1DUniform, x_pos)

Find in which cell of 1-D uniform grid the coordinate is located

# Positional arguments
* `grid1duniform`: the 1-D uniform grid
* `x_pos`: the 3-D coordinate vector for the which the cell index is to be determined (only the first component is used)
"""
@inline function get_cell(grid1duniform::Grid1DUniform, x_pos)
    return floor(Int64, x_pos[1] * grid1duniform.inv_Δx) + 1
end

"""
    sample_particles_equal_weight!(rng, grid1duniform, particles, pia, species, species_data, ppc::Integer, T, Fnum)

Sample particles from a Maxwellian distribution in each cell of 1-D uniform grid given the number of particles per cell

# Positional arguments
* `rng`: the random number generator
* `grid1duniform`: the 1-D uniform grid
* `particles`: the `ParticleVector` of particles
* `pia`: the `ParticleIndexerArray`
* `species`: the index of the species to be sampled for
* `species_data`: `Vector` of `Species` data
* `ppc`: number of particles per cell to be sampled
* `T`: the temperature
* `Fnum`: the computational weight of the particles
"""
function sample_particles_equal_weight!(rng, grid1duniform, particles, pia, species,
                                        species_data, ppc::Integer, T, Fnum)

    for cell in 1:grid1duniform.n_cells
        sample_particles_equal_weight!(rng, particles, pia, cell, species,
                                       ppc, species_data[species].mass, T, Fnum,
                                       grid1duniform.cells[cell].xlo, grid1duniform.cells[cell].xhi,
                                       0.0, 1.0,
                                       0.0, 1.0;
                                       distribution=:Maxwellian, vx0=0.0, vy0=0.0, vz0=0.0)
    end
end

"""
    sample_particles_equal_weight!(rng, grid1duniform, particles, pia, species, species_data, ndens::Float64, T, Fnum)

Sample particles from a Maxwellian distribution in each cell of 1-D uniform grid given the target number density.
If the computed number of particles is not an integer value, the fractional remainder is used to
probabilistically sample an extra particle, so that on average, the expected number density is achieved.

# Positional arguments
* `rng`: the random number generator
* `grid1duniform`: the 1-D uniform grid
* `particles`: the `ParticleVector` of particles
* `pia`: the `ParticleIndexerArray`
* `species`: the index of the species to be sampled for
* `species_data`: `Vector` of `Species` data
* `ndens`: target number density
* `T`: the temperature
* `Fnum`: the computational weight of the particles
"""
function sample_particles_equal_weight!(rng, grid1duniform, particles, pia, species,
                                        species_data, ndens::Float64, T, Fnum)
    # ndens is target ndensity per cell
    for cell in 1:grid1duniform.n_cells
        n_in_cell = ndens * grid1duniform.cells[cell].V

        ppc = n_in_cell / Fnum
        ppc_int = floor(Int64, ppc)
        remainder = ppc - ppc_int

        R = rand(rng)
        if R < remainder
            ppc_int += 1
        end

        sample_particles_equal_weight!(rng, particles, pia, cell, species,
                                       ppc_int, species_data[species].mass, T, Fnum,
                                       grid1duniform.cells[cell].xlo, grid1duniform.cells[cell].xhi,
                                       0.0, 1.0,
                                       0.0, 1.0;
                                       distribution=:Maxwellian, vx0=0.0, vy0=0.0, vz0=0.0)
    end
end