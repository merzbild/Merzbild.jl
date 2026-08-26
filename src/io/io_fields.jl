using NCDatasets
using NetCDF

"""
    IOSkipListField

Struct that holds track of which variables are not to be written to NetCDF file for the
electrostatic field quantities computed on the nodes of a grid.
If the field value is `true`, the corresponding field quantity will not be output to the file.

# Fields
* `skip_charge_density`: whether the output of the charge density should be skipped
* `skip_potential`: whether the output of the electrostatic potential should be skipped
* `skip_electric_field`: whether the output of the electric field should be skipped
"""
struct IOSkipListField
    skip_charge_density::Bool
    skip_potential::Bool
    skip_electric_field::Bool

    @doc """
        IOSkipListField(list_of_variables_to_skip)

    Construct an `IOSkipListField` from a list of variable names.
    The possible names are: `rho` or `charge_density`, `phi` or `potential`,
    `E` or `electric_field`.

    # Positional arguments
    * `list_of_variables_to_skip`: list of variable names to skip
    """
    function IOSkipListField(list_of_variables_to_skip)

        skip_charge_density = false
        skip_potential = false
        skip_electric_field = false

        if ("rho" in list_of_variables_to_skip) || ("charge_density" in list_of_variables_to_skip)
            skip_charge_density = true
        end

        if ("phi" in list_of_variables_to_skip) || ("potential" in list_of_variables_to_skip)
            skip_potential = true
        end

        if ("E" in list_of_variables_to_skip) || ("electric_field" in list_of_variables_to_skip)
            skip_electric_field = true
        end

        return new(skip_charge_density, skip_potential, skip_electric_field)
    end

    @doc """
        IOSkipListField()

    Construct an empty `IOSkipListField`.
    """
    function IOSkipListField()
        return IOSkipListField([])
    end
end

"""
    NCDataHolderField <: AbstractNCDataHolder

Struct that holds NetCDF-output related data for the I/O of the electrostatic field quantities
stored in the nodes of a grid.

# Fields
* `filehandle`: handle to the open NetCDF file
* `timestep_dim`: timestep dimension that used to keep track of the number of output steps
* `v_timestep`: variable to hold the simulation timestep number (dimension `time`)
* `v_charge_density`: variable to hold the charge density (dimension `n_nodes x time`)
* `v_potential`: variable to hold the electrostatic potential (dimension `n_nodes x time`)
* `v_electric_field`: variable to hold the x-component of the electric field (dimension `n_nodes x time`)
* `n_nodes`: number of nodes
* `n_nodes_1`: constant vector `[n_nodes, 1]` (used for offsets during I/O)
* `currtimesteps`: vector `[n_t_output]`, where `n_t_output` is the current output timestep (i.e. how many times the
    properties have already been output, not the simulation timestep) (used for offsets during I/O)
* `currtimesteps_1`: vector `[1, n_t_output]`, where `n_t_output` is the current output timestep (i.e. how many times the
    properties have already been output, not the simulation timestep) (used for offsets during I/O)
* `timestep`: vector storing the current simulation timestep
* `skip_list`: `IOSkipListField` instance of variables to skip during output
"""
mutable struct NCDataHolderField <: AbstractNCDataHolder
    filehandle::NcFile
    timestep_dim::NcDim  # timestep dimension, used to keep track of where we are in the file
    v_timestep::NcVar  # timestep

    v_charge_density::NcVar  # charge density: "n_nodes" x "time"
    v_potential::NcVar  # potential: "n_nodes" x "time"
    v_electric_field::NcVar  # electric field: "n_nodes" x "time"

    n_nodes::Int64

    # some constant offsets of ones (to count number of written elements)
    n_nodes_1::Vector{Int64}
    currtimesteps::Vector{Int64}
    currtimesteps_1::Vector{Int64}
    timestep::Vector{Float64}

    skip_list::IOSkipListField

    @doc """
        NCDataHolderField(nc_filename, names_skip_list, field_props; global_attributes=Dict{Any,Any}(), mode=NC_64BIT_OFFSET)

    Construct a `NCDataHolderField` instance with a list of variables to skip.

    # Positional arguments
    * `nc_filename`: filename to write output to
    * `names_skip_list`: list of variable names to skip, see [`IOSkipListField`](@ref) for more details
    * `field_props`: the `ElectrostaticFieldProps` instance which will be used for the output of the field quantities

    # Keyword arguments
    * `global_attributes`: dictionary of any additional attributes to write to the netCDF file as a global attribute
    * `mode`: NetCDF file format mode (default: `NC_64BIT_OFFSET` for older and faster format, `NC_NETCDF4` for NetCDF4 format)
    """
    function NCDataHolderField(nc_filename, names_skip_list, field_props; global_attributes=Dict{Any,Any}(), mode=NC_64BIT_OFFSET)
        skip_list = IOSkipListField(names_skip_list)

        gatts = deepcopy(global_attributes)

        nodes_dim = NcDim("n_nodes", field_props.n_nodes, unlimited=false)
        timestep_dim = NcDim("timestep", 0, unlimited=true)

        v_timestep = NcVar("timestep", [timestep_dim], t=Float64, compress=-1)

        v_charge_density = NcVar("charge_density", [nodes_dim, timestep_dim], t=Float64, compress=-1)
        v_potential = NcVar("potential", [nodes_dim, timestep_dim], t=Float64, compress=-1)
        v_electric_field = NcVar("electric_field", [nodes_dim, timestep_dim], t=Float64, compress=-1)

        varlist = NetCDF.NcVar[v_timestep]

        if !skip_list.skip_charge_density
            push!(varlist, v_charge_density)
        end

        if !skip_list.skip_potential
            push!(varlist, v_potential)
        end

        if !skip_list.skip_electric_field
            push!(varlist, v_electric_field)
        end

        filehandle = NetCDF.create(nc_filename, varlist, gatts=gatts, mode=mode)

        return new(filehandle,
                   timestep_dim, v_timestep,
                   v_charge_density, v_potential, v_electric_field,
                   field_props.n_nodes,
                   [field_props.n_nodes, 1], [1], [1, 1], [0.0],
                   skip_list)
    end

    @doc """
        NCDataHolderField(nc_filename, field_props; global_attributes=Dict{Any,Any}(), mode=NC_64BIT_OFFSET)

    Construct a `NCDataHolderField` instance with an empty list of variables to skip.

    # Positional arguments
    * `nc_filename`: filename to write output to
    * `field_props`: the `ElectrostaticFieldProps` instance which will be used for the output of the field quantities

    # Keyword arguments
    * `global_attributes`: dictionary of any additional attributes to write to the netCDF file as a global attribute
    * `mode`: NetCDF file format mode (default: `NC_64BIT_OFFSET` for older and faster format, `NC_NETCDF4` for NetCDF4 format)
    """
    function NCDataHolderField(nc_filename, field_props; global_attributes=Dict{Any,Any}(), mode=NC_64BIT_OFFSET)
        return NCDataHolderField(nc_filename, [], field_props; global_attributes=global_attributes, mode=mode)
    end
end

"""
    write_netcdf(ds, field_props::ElectrostaticFieldProps, timestep; sync_freq=0)

Write the electrostatic field quantities stored in an `ElectrostaticFieldProps` instance to a
NetCDF file and synchronize file to disk if necessary.

# Positional arguments
* `ds`: the `NCDataHolderField` for the file to which the output will be written
* `field_props`: the `ElectrostaticFieldProps` instance containing the field quantities
* `timestep`: the simulation timestep

# Keyword arguments
* `sync_freq`: if larger than 0 and if the number of timesteps output is proportional to `sync_freq`,
    the data will be synchronized to disk. If set to 1, will sync data to disk at every timestep at which
    data is written to the file.
"""
function write_netcdf(ds, field_props::ElectrostaticFieldProps, timestep; sync_freq=0)
    currtimesteps = ds.timestep_dim.dimlen + 1

    @inbounds ds.currtimesteps[1] = currtimesteps
    @inbounds ds.currtimesteps_1[2] = currtimesteps
    @inbounds ds.timestep[1] = timestep

    NetCDF.putvar(ds.v_timestep, ds.timestep, start=ds.currtimesteps)

    if !ds.skip_list.skip_charge_density
        NetCDF.putvar(ds.v_charge_density, field_props.charge_density, start=ds.currtimesteps_1, count=ds.n_nodes_1)
    end

    if !ds.skip_list.skip_potential
        NetCDF.putvar(ds.v_potential, field_props.potential, start=ds.currtimesteps_1, count=ds.n_nodes_1)
    end

    if !ds.skip_list.skip_electric_field
        NetCDF.putvar(ds.v_electric_field, field_props.electric_field, start=ds.currtimesteps_1, count=ds.n_nodes_1)
    end

    if (sync_freq > 0) && (currtimesteps % sync_freq == 0)
        NetCDF.sync(ds.filehandle)
    end
end
