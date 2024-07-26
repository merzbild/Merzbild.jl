using NCDatasets
using NetCDF

mutable struct NCDataHolder
    filehandle::NcFile
    timestep_dim::NcDim  # timestep dimension, used to keep track of where we are in the file
    v_spn::NcVar  # species names: "n_species"
    v_timestep::NcVar  # timestep
    v_lpa::NcVar  # length of particle array: "n_species" x "time"
    v_mompows::NcVar  # moment powers: "n_moments"
    v_moments::NcVar  # moment values: "n_moments" x "n_cells" x "n_species" x "time"
    v_np::NcVar  # number of particles: "n_cells" x "n_species" x "time"
    v_ndens::NcVar  # number density: "n_cells" x "n_species" x "time"
    v_v::NcVar  # velocity: 3 x "n_cells" x "n_species" x "time"
    v_T::NcVar  # temperature: "n_cells" x "n_species" x "time"
    n_species_1::Vector{Int64}
    n_cells_n_species_1::Vector{Int64}
    n_v_n_cells_n_species_1::Vector{Int64}
    currtimesteps::Vector{Int64}
    currtimesteps_1::Vector{Int64}
    currtimesteps_1_1::Vector{Int64}
    currtimesteps_1_1_1::Vector{Int64}
    timestep::Vector{Float64}
end

struct IOSkipList
    skip_length_particle_array::Bool
    skip_moments::Bool
    skip_number_of_particles::Bool
    skip_number_density::Bool
    skip_velocity::Bool
    skip_temperature::Bool
end

function create_IO_skip_list(list_of_variables_to_skip)

    skip_length_particle_array = false
    skip_moments = false
    skip_number_of_particles = false
    skip_number_density = false
    skip_velocity = false
    skip_temperature = false

    if "length_particle_array" in list_of_variables_to_skip
        skip_length_particle_array = true
    end

    if "moments" in list_of_variables_to_skip
        skip_moments = true
    end

    if "np" in list_of_variables_to_skip
        skip_number_of_particles = true
    end

    if "ndens" in list_of_variables_to_skip
        skip_number_density = true
    end

    if "v" in list_of_variables_to_skip
        skip_velocity = true
    end

    if "T" in list_of_variables_to_skip
        skip_temperature = true
    end

    return IOSkipList(skip_length_particle_array, skip_moments, skip_number_of_particles,
                      skip_number_density, skip_velocity, skip_temperature)
end

function close_netcdf(ds::NCDataHolder)
    finalize(ds.filehandle)
end

function create_netcdf_phys_props(nc_filename, species_data, phys_props; global_attributes=Dict{Any,Any}())
    v_dim = NcDim("vel_components", 3, unlimited=false)
    cells_dim = NcDim("n_cells", phys_props.n_cells, unlimited=false)
    species_dim = NcDim("n_species", phys_props.n_species, unlimited=false)
    moments_dim = NcDim("n_moments", phys_props.n_moments, unlimited=false)
    timestep_dim = NcDim("timestep", 0, unlimited=true)

    v_spn = NcVar("species_names", [species_dim], t=String, compress=-1)
    v_timestep = NcVar("timestep", [timestep_dim], t=Float64, compress=-1)
    v_lpa = NcVar("length_particle_array", [species_dim, timestep_dim], t=Float64, compress=-1)
    v_mompows = NcVar("moment_powers", [moments_dim], t=Int32, compress=-1)
    v_moments = NcVar("moments", [moments_dim, cells_dim, species_dim, timestep_dim], t=Float64, compress=-1)
    v_np = NcVar("np", [cells_dim, species_dim, timestep_dim], t=Float64, compress=-1)
    v_ndens = NcVar("ndens", [cells_dim, species_dim, timestep_dim], t=Float64, compress=-1)
    v_v = NcVar("v", [v_dim, cells_dim, species_dim, timestep_dim], t=Float64, compress=-1)
    v_T = NcVar("T", [cells_dim, species_dim, timestep_dim], t=Float64, compress=-1)

    varlist = [v_spn, v_timestep, v_lpa, v_mompows, v_moments, v_np, v_ndens, v_v, v_T]

    filehandle = NetCDF.create(nc_filename, varlist, gatts=global_attributes, mode=NC_NETCDF4)

    NetCDF.putvar(v_spn, [species.name for species in species_data])
    NetCDF.putvar(v_mompows, phys_props.moment_powers)

    return NCDataHolder(filehandle, timestep_dim, v_spn, v_timestep, v_lpa, v_mompows, v_moments, v_np, v_ndens, v_v, v_T,
                        [phys_props.n_species, 1], [phys_props.n_cells, phys_props.n_species, 1],
                        [3, phys_props.n_cells, phys_props.n_species, 1], [1], [1, 1], [1, 1, 1], [1, 1, 1, 1], [0.0])
end

function write_netcdf_phys_props(ds, phys_props, timestep; sync_freq=0)
    currtimesteps = ds.timestep_dim.dimlen + 1

    ds.currtimesteps[1] = currtimesteps
    ds.currtimesteps_1[2] = currtimesteps
    ds.currtimesteps_1_1[3] = currtimesteps
    ds.currtimesteps_1_1_1[4] = currtimesteps
    ds.timestep[1] = timestep

    NetCDF.putvar(ds.v_timestep, ds.timestep, start=ds.currtimesteps)

    NetCDF.putvar(ds.v_lpa, phys_props.lpa, start=ds.currtimesteps_1, count=ds.n_species_1)
    NetCDF.putvar(ds.v_np, phys_props.np, start=ds.currtimesteps_1_1, count=ds.n_cells_n_species_1)
    NetCDF.putvar(ds.v_ndens, phys_props.n, start=ds.currtimesteps_1_1, count=ds.n_cells_n_species_1)
    NetCDF.putvar(ds.v_v, phys_props.v, start=ds.currtimesteps_1_1_1, count=ds.n_v_n_cells_n_species_1)
    NetCDF.putvar(ds.v_T, phys_props.T, start=ds.currtimesteps_1_1, count=ds.n_cells_n_species_1)

    if (phys_props.n_moments > 0)
        NetCDF.putvar(ds.v_moments, phys_props.moments, start=[1, 1, 1, currtimesteps], count=[phys_props.n_moments, phys_props.n_cells, phys_props.n_species, 1])
    end

    if (sync_freq > 0) && (currtimesteps % sync_freq == 0)
        NetCDF.sync(ds.filehandle)
    end
end

function write_netcdf_phys_props_nov(ds, phys_props, timestep; sync_freq=0)
    currtimesteps = ds.timestep_dim.dimlen + 1

    ds.currtimesteps[1] = currtimesteps
    ds.currtimesteps_1[2] = currtimesteps
    ds.currtimesteps_1_1[3] = currtimesteps
    ds.currtimesteps_1_1_1[4] = currtimesteps
    ds.timestep[1] = timestep

    NetCDF.putvar(ds.v_timestep, ds.timestep, start=ds.currtimesteps)

    NetCDF.putvar(ds.v_lpa, phys_props.lpa, start=ds.currtimesteps_1, count=ds.n_species_1)
    NetCDF.putvar(ds.v_np, phys_props.np, start=ds.currtimesteps_1_1, count=ds.n_cells_n_species_1)
    NetCDF.putvar(ds.v_ndens, phys_props.n, start=ds.currtimesteps_1_1, count=ds.n_cells_n_species_1)

    if (sync_freq > 0) && (currtimesteps % sync_freq == 0)
        NetCDF.sync(ds.filehandle)
    end
end


function write_netcdf_phys_props(ds, skip_list, phys_props, timestep; sync_freq=0)
    currtimesteps = ds.timestep_dim.dimlen + 1

    ds.currtimesteps[1] = currtimesteps
    ds.currtimesteps_1[2] = currtimesteps
    ds.currtimesteps_1_1[3] = currtimesteps
    ds.currtimesteps_1_1_1[4] = currtimesteps
    ds.timestep[1] = timestep

    NetCDF.putvar(ds.v_timestep, ds.timestep, start=ds.currtimesteps)

    if !skip_list.skip_length_particle_array
        NetCDF.putvar(ds.v_lpa, phys_props.lpa, start=ds.currtimesteps_1, count=ds.n_species_1)
    end

    if !skip_list.skip_number_of_particles
        NetCDF.putvar(ds.v_np, phys_props.np, start=ds.currtimesteps_1_1, count=ds.n_cells_n_species_1)
    end

    if !skip_list.skip_number_density
        NetCDF.putvar(ds.v_ndens, phys_props.n, start=ds.currtimesteps_1_1, count=ds.n_cells_n_species_1)
    end

    if !skip_list.skip_velocity
        NetCDF.putvar(ds.v_v, phys_props.v, start=ds.currtimesteps_1_1_1, count=ds.n_v_n_cells_n_species_1)
    end

    if !skip_list.skip_temperature
        NetCDF.putvar(ds.v_T, phys_props.T, start=ds.currtimesteps_1_1, count=ds.n_cells_n_species_1)
    end

    if !skip_list.skip_moments
        if (phys_props.n_moments > 0)
            NetCDF.putvar(ds.v_moments, phys_props.moments, start=[1, 1, 1, currtimesteps], count=[phys_props.n_moments, phys_props.n_cells, phys_props.n_species, 1])
        end
    end

    if (sync_freq > 0) && (currtimesteps % sync_freq == 0)
        NetCDF.sync(ds.filehandle)
    end
end
