struct Cell1D
    xlo::Float64
    xhi::Float64
    V::Float64
end

struct Grid1DUniform
    L::Float64
    n_cells::Int64
    Δx::Float64
    inv_Δx::Float64
    cells::Vector{Cell1D}
    min_x::Float64  # so that we don't get particles stuck exactly at the wall
    max_x::Float64  # so that we don't get particles stuck exactly at the wall
end

function create_grid1D_uniform(L, nx)
    cells = Vector{Cell1D}(undef, nx)
    dx = L / nx

    for i in 1:nx
        xlo = (i-1) * dx
        xhi = i * dx
        V = dx
        cells[i] = Cell1D(xlo, xhi, V)
    end

    return Grid1DUniform(L, nx, dx, 1.0 / dx, cells, dx * 1e-12, (1.0 - 1e-12) * L)
end

function get_cell(grid1duniform::Grid1DUniform, x_pos)
    return floor(Int64, x_pos[1] * grid1duniform.inv_Δx) + 1
end

function sample_particles_equal_weight!(rng, grid1duniform, particles, pia, species,
                                        species_data, ppc::Int64, T, Fnum)

    for cell in 1:grid1duniform.n_cells
        sample_particles_equal_weight!(rng, particles, pia, cell, species,
                                       ppc, species_data[species].mass, T, Fnum,
                                       grid1duniform.cells[cell].xlo, grid1duniform.cells[cell].xhi,
                                       0.0, 1.0,
                                       0.0, 1.0;
                                       distribution=:Maxwellian, vx0=0.0, vy0=0.0, vz0=0.0)
    end
end

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