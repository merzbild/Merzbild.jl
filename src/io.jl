using NCDatasets
using NetCDF


struct IOSkipList
    skip_length_particle_array::Bool
    skip_moments::Bool
    skip_number_of_particles::Bool
    skip_number_density::Bool
    skip_velocity::Bool
    skip_temperature::Bool

    function IOSkipList(list_of_variables_to_skip)

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
    
        return new(skip_length_particle_array, skip_moments, skip_number_of_particles,
                   skip_number_density, skip_velocity, skip_temperature)
    end


    function IOSkipList()
        return IOSkipList([])
    end
end

mutable struct NCDataHolder
    filehandle::NcFile
    ndens_not_Np::Bool
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

    # some constant offsets of ones (to count number of written elements)
    n_species_1::Vector{Int64}
    n_cells_n_species_1::Vector{Int64}
    n_v_n_cells_n_species_1::Vector{Int64}
    currtimesteps::Vector{Int64}
    currtimesteps_1::Vector{Int64}
    currtimesteps_1_1::Vector{Int64}
    currtimesteps_1_1_1::Vector{Int64}
    timestep::Vector{Float64}

    skip_list::IOSkipList

    function NCDataHolder(nc_filename, names_skip_list, species_data, phys_props; global_attributes=Dict{Any,Any}())
        skip_list = IOSkipList(names_skip_list)

        gatts = deepcopy(global_attributes)

        if phys_props.ndens_not_Np
            gatts["COMMENT ndens"] = "In this file, the ndens variable stores the number density, NOT # of physical particles in cell"
        else
            gatts["COMMENT ndens"] = "In this file, the ndens variable stores # of physical particles in cell, NOT number density"
        end

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

        varlist::Vector{NetCDF.NcVar} = [v_spn, v_timestep]

        if !skip_list.skip_length_particle_array
            push!(varlist, v_lpa)
        end

        if !skip_list.skip_moments
            push!(varlist, v_mompows)
            push!(varlist, v_moments)
        end

        if !skip_list.skip_number_of_particles
            push!(varlist, v_np)
        end

        if !skip_list.skip_number_density
            push!(varlist, v_ndens)
        end

        if !skip_list.skip_velocity
            push!(varlist, v_v)
        end

        if !skip_list.skip_temperature
            push!(varlist, v_T)
        end

        filehandle = NetCDF.create(nc_filename, varlist, gatts=gatts, mode=NC_NETCDF4)

        NetCDF.putvar(v_spn, [species.name for species in species_data])

        if !skip_list.skip_moments
            NetCDF.putvar(v_mompows, phys_props.moment_powers)
        end

        return new(filehandle, phys_props.ndens_not_Np,
                   timestep_dim, v_spn, v_timestep, v_lpa, v_mompows, v_moments, v_np, v_ndens, v_v, v_T,
                   [phys_props.n_species, 1], [phys_props.n_cells, phys_props.n_species, 1],
                   [3, phys_props.n_cells, phys_props.n_species, 1], [1], [1, 1], [1, 1, 1], [1, 1, 1, 1], [0.0],
                   skip_list)
    end

    function NCDataHolder(nc_filename, species_data, phys_props; global_attributes=Dict{Any,Any}())
        return NCDataHolder(nc_filename, [], species_data, phys_props; global_attributes=global_attributes)
    end
end

function close_netcdf(ds::NCDataHolder)
    finalize(ds.filehandle)
end

# function write_netcdf_phys_props(ds, phys_props, timestep; sync_freq=0)
#     if phys_props.ndens_not_Np != ds.ndens_not_Np
#         throw(ErrorException("Inconsistent computation of ndens/number of physical particles in cell"))
#     end

#     currtimesteps = ds.timestep_dim.dimlen + 1

#     ds.currtimesteps[1] = currtimesteps
#     ds.currtimesteps_1[2] = currtimesteps
#     ds.currtimesteps_1_1[3] = currtimesteps
#     ds.currtimesteps_1_1_1[4] = currtimesteps
#     ds.timestep[1] = timestep

#     NetCDF.putvar(ds.v_timestep, ds.timestep, start=ds.currtimesteps)

#     NetCDF.putvar(ds.v_lpa, phys_props.lpa, start=ds.currtimesteps_1, count=ds.n_species_1)
#     NetCDF.putvar(ds.v_np, phys_props.np, start=ds.currtimesteps_1_1, count=ds.n_cells_n_species_1)
#     NetCDF.putvar(ds.v_ndens, phys_props.n, start=ds.currtimesteps_1_1, count=ds.n_cells_n_species_1)
#     NetCDF.putvar(ds.v_v, phys_props.v, start=ds.currtimesteps_1_1_1, count=ds.n_v_n_cells_n_species_1)
#     NetCDF.putvar(ds.v_T, phys_props.T, start=ds.currtimesteps_1_1, count=ds.n_cells_n_species_1)

#     if (phys_props.n_moments > 0)
#         NetCDF.putvar(ds.v_moments, phys_props.moments, start=[1, 1, 1, currtimesteps], count=[phys_props.n_moments, phys_props.n_cells, phys_props.n_species, 1])
#     end

#     if (sync_freq > 0) && (currtimesteps % sync_freq == 0)
#         NetCDF.sync(ds.filehandle)
#     end
# end

function write_netcdf_phys_props(ds, phys_props, timestep; sync_freq=0)
    currtimesteps = ds.timestep_dim.dimlen + 1

    ds.currtimesteps[1] = currtimesteps
    ds.currtimesteps_1[2] = currtimesteps
    ds.currtimesteps_1_1[3] = currtimesteps
    ds.currtimesteps_1_1_1[4] = currtimesteps
    ds.timestep[1] = timestep

    NetCDF.putvar(ds.v_timestep, ds.timestep, start=ds.currtimesteps)

    if !ds.skip_list.skip_length_particle_array
        NetCDF.putvar(ds.v_lpa, phys_props.lpa, start=ds.currtimesteps_1, count=ds.n_species_1)
    end

    if !ds.skip_list.skip_number_of_particles
        NetCDF.putvar(ds.v_np, phys_props.np, start=ds.currtimesteps_1_1, count=ds.n_cells_n_species_1)
    end

    if !ds.skip_list.skip_number_density
        if phys_props.ndens_not_Np != ds.ndens_not_Np
            throw(ErrorException("Inconsistent computation of ndens/number of physical particles in cell"))
        end
        NetCDF.putvar(ds.v_ndens, phys_props.n, start=ds.currtimesteps_1_1, count=ds.n_cells_n_species_1)
    end

    if !ds.skip_list.skip_velocity
        NetCDF.putvar(ds.v_v, phys_props.v, start=ds.currtimesteps_1_1_1, count=ds.n_v_n_cells_n_species_1)
    end

    if !ds.skip_list.skip_temperature
        NetCDF.putvar(ds.v_T, phys_props.T, start=ds.currtimesteps_1_1, count=ds.n_cells_n_species_1)
    end

    if !ds.skip_list.skip_moments
        if (phys_props.n_moments > 0)
            NetCDF.putvar(ds.v_moments, phys_props.moments, start=[1, 1, 1, currtimesteps], count=[phys_props.n_moments, phys_props.n_cells, phys_props.n_species, 1])
        end
    end

    if (sync_freq > 0) && (currtimesteps % sync_freq == 0)
        NetCDF.sync(ds.filehandle)
    end
end
