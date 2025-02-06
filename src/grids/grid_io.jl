@doc """
    write_grid(nc_filename, grid::Grid1DUniform; global_attributes=Dict{Any,Any}())

Write grid info to a NetCDF file

# Positional arguments
* `nc_filename`: path to the NetCDF file
* `grid`: the grid to be written out to a file

# Keyword arguments
* `global_attributes`: dictionary of additional global attributes to be written to the netCDF file
"""
function write_grid(nc_filename, grid::Grid1DUniform; global_attributes=Dict{Any,Any}())
    cells_dim = NcDim("cells", grid.n_cells, unlimited=false)

    v_xlo = NcVar("xlo", [cells_dim], t=Float64, compress=-1)
    v_xhi = NcVar("xhi", [cells_dim], t=Float64, compress=-1)
    v_V = NcVar("cell_volume", [cells_dim], t=Float64, compress=-1)
    varlist = [v_xlo, v_xhi, v_V]

    gatts = deepcopy(global_attributes)
    gatts["dx"] = grid.Δx
    gatts["L"] = grid.L
    gatts["dim"] = 1

    filehandle = NetCDF.create(nc_filename, varlist, gatts=gatts, mode=NC_NETCDF4)

    filehandle["xlo"][:] = [grid.cells[i].xlo for i in 1:grid.n_cells]
    filehandle["xhi"][:] = [grid.cells[i].xhi for i in 1:grid.n_cells]
    filehandle["cell_volume"][:] = [grid.cells[i].V for i in 1:grid.n_cells]

    NetCDF.sync(filehandle)
    finalize(filehandle)
end