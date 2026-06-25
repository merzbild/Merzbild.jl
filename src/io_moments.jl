using NCDatasets
using NetCDF

"""
    NCDataHolderMoments <: AbstractNCDataHolder

Struct that holds NetCDF-output related data for moments I/O.

# Fields
* `filehandle`: handle to the open NetCDF file
* `timestep_dim`: timestep dimension that used to keep track of the number of output steps
* `v_timestep`: variable to hold the simulation timestep number (dimension `time`)
* `v_mompows`: variable to hold list of total moment powers (dimension `n_moments`)
* `v_moments`: variable to hold total moments (dimension `n_moments x n_cells x n_species x time`)
* `n_moments`: number of moments
* `n_cells`: number of cells
* `n_species`: number of species
* `n_moments_n_cells_n_species_1`: constant vector `[n_moments, n_cells, n_species, 1]` (used for offsets during I/O)
* `currtimesteps`: vector `[n_t_output]`, where `n_t_output` is the current output timestep (i.e. how many times the
    properties have already been output, not the simulation timestep) (used for offsets during I/O)
* `currtimesteps_1_1_1`: vector `[1, 1, 1, n_t_output]`, where `n_t_output` is the current output timestep (i.e. how many times the
    properties have already been output, not the simulation timestep) (used for offsets during I/O)
* `timestep`: vector storing the current simulation timestep
"""
mutable struct NCDataHolderMoments <: AbstractNCDataHolder
    filehandle::NcFile
    timestep_dim::NcDim  # timestep dimension, used to keep track of where we are in the file
    v_timestep::NcVar  # timestep
    v_mompows::NcVar  # moment powers: "n_moments"
    v_moments::NcVar  # moment values: "n_moments" x "n_cells" x "n_species" x "time"

    n_moments::Int64
    n_cells::Int64
    n_species::Int64

    # some constant offsets of ones (to count number of written elements)
    n_moments_n_cells_n_species_1::Vector{Int64}
    currtimesteps::Vector{Int64}
    currtimesteps_1_1_1::Vector{Int64}
    timestep::Vector{Float64}

    @doc """
        NCDataHolderMoments(nc_filename, species_data, n_cells, n_species, moment_powers; global_attributes=Dict{Any,Any}(), mode=NC_64BIT_OFFSET)

    Construct a `NCDataHolderMoments` instance.

    # Positional arguments
    * `nc_filename`: filename to write output to
    * `species_data`: the vector of `Species` data for the species in the simulation
    * `n_cells`: number of cells
    * `n_species`: number of species
    * `moment_powers`: vector of moment powers (e.g., `Int8[2, 4, 6]`)
    
    # Keyword arguments
    * `global_attributes`: dictionary of any additional attributes to write to the netCDF file as a global attribute
    * `mode`: NetCDF file format mode (default: `NC_64BIT_OFFSET` for older and faster format, can use `NC_NETCDF4` for NetCDF4 format)
    """
    function NCDataHolderMoments(nc_filename, species_data, n_cells, n_species, moment_powers; global_attributes=Dict{Any,Any}(), mode=NC_64BIT_OFFSET)
        gatts = deepcopy(global_attributes)
        gatts["species_names"] = join([species.name for species in species_data], ",")

        n_moments = length(moment_powers)

        moments_dim = NcDim("n_moments", n_moments, unlimited=false)
        cells_dim = NcDim("n_cells", n_cells, unlimited=false)
        species_dim = NcDim("n_species", n_species, unlimited=false)
        timestep_dim = NcDim("timestep", 0, unlimited=true)

        v_timestep = NcVar("timestep", [timestep_dim], t=Float64)
        v_mompows = NcVar("moment_powers", [moments_dim], t=Int32)
        v_moments = NcVar("moments", [moments_dim, cells_dim, species_dim, timestep_dim], t=Float64)

        varlist::Vector{NetCDF.NcVar} = [v_timestep, v_mompows, v_moments]

        filehandle = NetCDF.create(nc_filename, varlist, gatts=gatts, mode=mode)

        NetCDF.putvar(v_mompows, moment_powers)

        return new(filehandle, 
                   timestep_dim, v_timestep, v_mompows, v_moments,
                   n_moments, n_cells, n_species,
                   [n_moments, n_cells, n_species, 1],
                   [1], [1, 1, 1, 1], [0.0])
    end
end

"""
    write_netcdf(ds::NCDataHolderMoments, moments, timestep; sync_freq=0)
    
Write 'moments' to a NetCDF file and synchronize file to disk if necessary.

# Positional arguments
* `ds`: the `NCDataHolderMoments` for the file to which the output will be written
* `moments`: the moments array (dimension `n_moments x n_cells x n_species`)
* `timestep`: the simulation timestep

# Keyword arguments
* `sync_freq`: if larger than 0 and if the number of timesteps output is proportional to `sync_freq`,
    the data will be synchronized to disk. If set to 1, will sync data to disk at every timestep at which
    data is written to the file.
"""
function write_netcdf(ds::NCDataHolderMoments, moments, timestep; sync_freq=0)
    currtimesteps = ds.timestep_dim.dimlen + 1

    @inbounds ds.currtimesteps[1] = currtimesteps
    @inbounds ds.currtimesteps_1_1_1[4] = currtimesteps
    @inbounds ds.timestep[1] = timestep

    NetCDF.putvar(ds.v_timestep, ds.timestep, start=ds.currtimesteps)
    NetCDF.putvar(ds.v_moments, moments, start=ds.currtimesteps_1_1_1, count=ds.n_moments_n_cells_n_species_1)

    if (sync_freq > 0) && (currtimesteps % sync_freq == 0)
        NetCDF.sync(ds.filehandle)
    end
end
