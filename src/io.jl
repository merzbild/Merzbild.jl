using NCDatasets
using NetCDF

"""
    IOSkipList

Struct that holds track of which variables are not to be written to NetCDF file for physical properties
computed on a grid.
If the field value is `true`, the corresponding physical grid property will not be output to the file.

# Fields
* `skip_length_particle_array`: whether the length of the particle array should be skipped
* `skip_moments`: whether the output of the total moments should be skipped
* `skip_number_of_particles`: whether the output of the number of particles should be skipped
* `skip_number_density`: whether the output of the number density/number of physical particles should be skipped
* `skip_velocity`: whether the output of the velocity should be skipped
* `skip_temperature`: whether the output of the temperature should be skipped
"""
struct IOSkipList
    skip_length_particle_array::Bool
    skip_moments::Bool
    skip_number_of_particles::Bool
    skip_number_density::Bool
    skip_velocity::Bool
    skip_temperature::Bool

    @doc """
        IOSkipList(list_of_variables_to_skip)
    
    Construct an `IOSkipList` from a list of variable names.
    The possible names are: `length_particle_array`, `moments`,
    `np` or `nparticles`, `ndens`, `v`, `T`.

    # Positional arguments
    * `list_of_variables_to_skip`: list of variable names to skip
    """
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
    
        if ("np" in list_of_variables_to_skip) || ("nparticles" in list_of_variables_to_skip)
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

    @doc """
        IOSkipList()
    
    Construct an empty `IOSkipList`.
    """
    function IOSkipList()
        return IOSkipList([])
    end
end

"""
    AbstractNCDataHolder

Abstract type that holds NetCDF-output related data for I/O
"""
abstract type AbstractNCDataHolder end

"""
    NCDataHolder <: AbstractNCDataHolder

Struct that holds NetCDF-output related data for physical properties (grid properties) I/O.

# Fields
* `filehandle`: handle to the open NetCDF file
* `ndens_not_Np`: whether the number density or the number of physical particles is being output
* `timestep_dim`: timestep dimension that used to keep track of the number of output steps
* `v_spn`: variable to hold species' names (dimension `n_species`)
* `v_timestep`: variable to hold the simulation timestep number (dimension `time`)
* `v_lpa`: variable to hold lengths of particle arrays (dimension `n_species x time`)
* `v_mompows`: variable to hold list of total moment powers (dimension `n_moments`)
* `v_moments`: variable to hold total moments (dimension `n_moments x n_cells x n_species x time`)
* `v_np`:  variable to hold number of particles (dimension ` n_cells x n_species x time`)
* `v_ndens`: variable to hold number density or the number of physical particles (dimension ` n_cells x n_species x time`)
* `v_v`: variable to hold velocity (dimension `3 x n_cells x n_species x time`)
* `v_T`: variable to hold temperature (dimension `n_cells x n_species x time`)
* `n_species_1`: constant vector `[n_species, 1]` (used for offsets during I/O)
* `n_cells_n_species_1`: constant vector `[n_cells, n_species, 1]` (used for offsets during I/O)
* `n_v_n_cells_n_species_1`: constant vector `[3, n_cells, n_species, 1]` (used for offsets during I/O)
* `currtimesteps`: vector `[n_t_output]`, where `n_t_output` is the current output timestep (i.e. how many times the
    properties have already been output, not the simulation timestep) (used for offsets during I/O)
* `currtimesteps_1`: vector `[1, n_t_output]`, where `n_t_output` is the current output timestep (i.e. how many times the
    properties have already been output, not the simulation timestep) (used for offsets during I/O)
* `currtimesteps_1_1`: vector `[1, 1, n_t_output]`, where `n_t_output` is the current output timestep (i.e. how many times the
    properties have already been output, not the simulation timestep) (used for offsets during I/O)
* `currtimesteps_1_1_1`: vector `[1, 1, 1, n_t_output]`, where `n_t_output` is the current output timestep (i.e. how many times the
    properties have already been output, not the simulation timestep) (used for offsets during I/O)
* `timestep`: vector storing the current simulation timestep
* `skip_list`: `IOSkipList` instance of variables to skip during output
"""
mutable struct NCDataHolder <: AbstractNCDataHolder
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

    @doc """
        NCDataHolder(nc_filename, names_skip_list, species_data, phys_props; global_attributes=Dict{Any,Any}())

    Construct a `NCDataHolder` instance with a list of variables to skip.

    # Positional arguments
    * `nc_filename`: filename to write output to
    * `names_skip_list`: list of variable names to skip, see [`IOSkipList`](@ref) for more details
    * `species_data`: the vector of `Species` data for the species in the simulation
    * `phys_props`: the `PhysProps` instance which will be used for the computation and output of physical properties
    
    # Keyword arguments
    * `global_attributes`: dictionary of any additional attributes to write to the netCDF file as a global attribute
    """
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

    @doc """
        NCDataHolder(nc_filename, species_data, phys_props; global_attributes=Dict{Any,Any}())

    Construct a `NCDataHolder` instance with an empty list of variable to skip,

    # Positional arguments
    * `nc_filename`: filename to write output to
    * `species_data`: the vector of `Species` data for the species in the simulation
    * `phys_props`: the `PhysProps` instance which will be used for the computation and output of physical properties
    
    # Keyword arguments
    * `global_attributes`: dictionary of any additional attributes to write to the netCDF file as a global attribute
    """
    function NCDataHolder(nc_filename, species_data, phys_props; global_attributes=Dict{Any,Any}())
        return NCDataHolder(nc_filename, [], species_data, phys_props; global_attributes=global_attributes)
    end
end

"""
    close_netcdf(ds::AbstractNCDataHolder)

Close NetCDF file.

# Positional arguments
* `ds`: an `AbstractNCDataHolder` instance to close.
"""
function close_netcdf(ds::AbstractNCDataHolder)
    finalize(ds.filehandle)
end

"""
     write_netcdf_phys_props(ds, phys_props, timestep; sync_freq=0)

Write computed `PhysProps` to NetCDF file and synchronize file to disk if necessary.

# Positional arguments
* `ds`: the `NCDataHolder` for the file to which the output will be written
* `phys_props`: the `PhysProps` instance containing the computed properties
* `timestep`: the simulation timestep

# Keyword arguments
* `sync_freq`: if larger than 0 and if the number of timesteps output is proportional to `sync_freq`,
    the data will be synchronized to disk. If set to 1, will sync data to disk at every timestep at which
    data is written to the file.

# Throws
`ErrorException` if the `NCDataHolder` expects number density and the `phys_props` holds the number of physical
particles, or vice versa.
"""
function write_netcdf_phys_props(ds, phys_props, timestep; sync_freq=0)
    currtimesteps = ds.timestep_dim.dimlen + 1

    @inbounds ds.currtimesteps[1] = currtimesteps
    @inbounds ds.currtimesteps_1[2] = currtimesteps
    @inbounds ds.currtimesteps_1_1[3] = currtimesteps
    @inbounds ds.currtimesteps_1_1_1[4] = currtimesteps
    @inbounds ds.timestep[1] = timestep

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


"""
    IOSkipListSurf

Struct that holds track of which variables are not to be written to NetCDF file
for computed surface properties.
If the field value is `true`, the corresponding physical grid property will not be output to the file.

# Fields
* `skip_number_of_particles`: whether the output of the number of particles should be skipped
* `skip_fluxes`: whether the output of the incident/reflected fluxes should be skipped
* `skip_force`: whether the output of the force should be skipped
* `skip_normal_pressure`: whether the output of the normal pressure should be skipped
* `skip_shear_pressure`: whether the output of the shear pressure should be skipped
* `skip_kinetic_energy_flux`: whether the output of the kinetic energy flux should be skipped
"""
struct IOSkipListSurf
    skip_number_of_particles::Bool
    skip_fluxes::Bool
    skip_force::Bool
    skip_normal_pressure::Bool
    skip_shear_pressure::Bool
    skip_kinetic_energy_flux::Bool

    @doc """
        IOSkipListSurf(list_of_variables_to_skip)
    
    Construct an `IOSkipListSurf` from a list of variable names.
    The possible names are:
    `np` or `nparticles`, `fluxes`, `force`, `normal_pressure`, `shear_pressure`, "kinetic_energy_flux".

    # Positional arguments
    * `list_of_variables_to_skip`: list of variable names to skip
    """
    function IOSkipListSurf(list_of_variables_to_skip)

        skip_number_of_particles = false
        skip_fluxes = false
        skip_force = false
        skip_normal_pressure = false
        skip_shear_pressure = false
        skip_kinetic_energy_flux = false
    
        if ("np" in list_of_variables_to_skip) || ("nparticles" in list_of_variables_to_skip)
            skip_number_of_particles = true
        end
    
        if "fluxes" in list_of_variables_to_skip
            skip_fluxes = true
        end
    
        if "force" in list_of_variables_to_skip
            skip_force = true
        end
    
        if "normal_pressure" in list_of_variables_to_skip
            skip_normal_pressure = true
        end
    
        if "shear_pressure" in list_of_variables_to_skip
            skip_shear_pressure = true
        end
    
        if "kinetic_energy_flux" in list_of_variables_to_skip
            skip_kinetic_energy_flux = true
        end
    
        return new(skip_number_of_particles, skip_fluxes, skip_force,
                   skip_normal_pressure, skip_shear_pressure, skip_kinetic_energy_flux)
    end

    @doc """
        IOSkipListSurf()
    
    Construct an empty `IOSkipListSurf`
    """
    function IOSkipListSurf()
        return IOSkipListSurf([])
    end
end

"""
    NCDataHolderSurf

Struct that holds NetCDF-output related data for surface properties I/O.

# Fields
* `filehandle`: handle to the open NetCDF file
* `timestep_dim`: timestep dimension that used to keep track of the number of output steps
* `v_spn`: variable to hold species' names (dimension `n_species`)
* `v_timestep`: variable to hold the simulation timestep number (dimension `time`)
* `v_np`:  variable to hold number of particles that hit the surface (dimension ` n_elements x n_species x time`)
* `v_flux_incident`: variable to hold incident mass flux (dimension `n_elements x n_species x time`)
* `v_flux_reflected`: variable to hold reflected mass flux (dimension `n_elements x n_species x time`)
* `v_force`: variable to hold force (dimension `3 x n_elements x n_species x time`)
* `v_normal_pressure`: variable to hold normal pressure (dimension `n_elements x n_species x time`)
* `v_shear_pressure`: variable to hold shear pressure (dimension `3 x n_elements x n_species x time`)
* `v_kinetic_energy_flux`: variable to hold kinetic energy flux (dimension `n_elements x n_species x time`)
* `n_species_1`: constant vector `[n_species, 1]` (used for offsets during I/O)
* `n_elements_n_species_1`: constant vector `[n_elements, n_species, 1]` (used for offsets during I/O)
* `n_v_n_elements_n_species_1`: constant vector `[3, n_elements, n_species, 1]` (used for offsets during I/O)
* `currtimesteps`: vector `[n_t_output]`, where `n_t_output` is the current output timestep (i.e. how many times the
    properties have already been output, not the simulation timestep) (used for offsets during I/O)
* `currtimesteps_1`: vector `[1, n_t_output]`, where `n_t_output` is the current output timestep (i.e. how many times the
    properties have already been output, not the simulation timestep) (used for offsets during I/O)
* `currtimesteps_1_1`: vector `[1, 1, n_t_output]`, where `n_t_output` is the current output timestep (i.e. how many times the
    properties have already been output, not the simulation timestep) (used for offsets during I/O)
* `currtimesteps_1_1_1`: vector `[1, 1, 1, n_t_output]`, where `n_t_output` is the current output timestep (i.e. how many times the
    properties have already been output, not the simulation timestep) (used for offsets during I/O)
* `timestep`: vector storing the current simulation timestep
* `skip_list`: `IOSkipListSurf` instance of variables to skip during output
"""
mutable struct NCDataHolderSurf <: AbstractNCDataHolder
    filehandle::NcFile
    timestep_dim::NcDim  # timestep dimension, used to keep track of where we are in the file
    v_spn::NcVar  # species names: "n_species"
    v_timestep::NcVar  # timestep

    v_np::NcVar  # number of particles: "n_elements" x "n_species" x "time"
    v_flux_incident::NcVar  # incident flux: "n_elements" x "n_species" x "time"
    v_flux_reflected::NcVar  # reflected flux: "n_elements" x "n_species" x "time"
    v_force::NcVar  # velocity: 3 x "n_elements" x "n_species" x "time"
    v_normal_pressure::NcVar  # normal pressure: "n_elements" x "n_species" x "time"
    v_shear_pressure::NcVar  # shear pressure: 3 x "n_elements" x "n_species" x "time"
    v_kinetic_energy_flux::NcVar  # kinetic energy flux: "n_elements" x "n_species" x "time"

    # some constant offsets of ones (to count number of written elements)
    n_species_1::Vector{Int64}
    n_elements_n_species_1::Vector{Int64}
    n_v_n_elements_n_species_1::Vector{Int64}
    currtimesteps::Vector{Int64}
    currtimesteps_1::Vector{Int64}
    currtimesteps_1_1::Vector{Int64}
    currtimesteps_1_1_1::Vector{Int64}
    timestep::Vector{Float64}

    skip_list::IOSkipListSurf

    @doc """
        NCDataHolderSurf(nc_filename, names_skip_list, species_data, surf_props; global_attributes=Dict{Any,Any}())

    Construct a `NCDataHolderSurf` instance with a list of variables to skip.

    # Positional arguments
    * `nc_filename`: filename to write output to
    * `names_skip_list`: list of variable names to skip, see [`IOSkipListSurf`](@ref) for more details
    * `species_data`: the vector of `Species` data for the species in the simulation
    * `surf_props`: the `SurfProps` instance which will be used for the computation and output of surface properties
    
    # Keyword arguments
    * `global_attributes`: dictionary of any additional attributes to write to the netCDF file as a global attribute
    """
    function NCDataHolderSurf(nc_filename, names_skip_list, species_data, surf_props; global_attributes=Dict{Any,Any}())
        skip_list = IOSkipListSurf(names_skip_list)

        gatts = deepcopy(global_attributes)

        v_dim = NcDim("vector_components", 3, unlimited=false)
        elements_dim = NcDim("n_elements", surf_props.n_elements, unlimited=false)
        species_dim = NcDim("n_species", surf_props.n_species, unlimited=false)
        timestep_dim = NcDim("timestep", 0, unlimited=true)

        v_spn = NcVar("species_names", [species_dim], t=String, compress=-1)
        v_timestep = NcVar("timestep", [timestep_dim], t=Float64, compress=-1)

        v_np = NcVar("np", [elements_dim, species_dim, timestep_dim], t=Float64, compress=-1)
        v_flux_incident = NcVar("flux_incident", [elements_dim, species_dim, timestep_dim], t=Float64, compress=-1)
        v_flux_reflected = NcVar("flux_reflected", [elements_dim, species_dim, timestep_dim], t=Float64, compress=-1)
        v_force = NcVar("force", [v_dim, elements_dim, species_dim, timestep_dim], t=Float64, compress=-1)
        v_normal_pressure = NcVar("normal_pressure", [elements_dim, species_dim, timestep_dim], t=Float64, compress=-1)
        v_shear_pressure = NcVar("shear_pressure", [v_dim, elements_dim, species_dim, timestep_dim], t=Float64, compress=-1)
        v_kinetic_energy_flux = NcVar("kinetic_energy_flux", [elements_dim, species_dim, timestep_dim], t=Float64, compress=-1)

        varlist::Vector{NetCDF.NcVar} = [v_spn, v_timestep]

        if !skip_list.skip_number_of_particles
            push!(varlist, v_np)
        end

        if !skip_list.skip_fluxes
            push!(varlist, v_flux_incident)
            push!(varlist, v_flux_reflected)
        end

        if !skip_list.skip_force
            push!(varlist, v_force)
        end

        if !skip_list.skip_normal_pressure
            push!(varlist, v_normal_pressure)
        end

        if !skip_list.skip_shear_pressure
            push!(varlist, v_shear_pressure)
        end

        if !skip_list.skip_kinetic_energy_flux
            push!(varlist, v_kinetic_energy_flux)
        end

        filehandle = NetCDF.create(nc_filename, varlist, gatts=gatts, mode=NC_NETCDF4)

        NetCDF.putvar(v_spn, [species.name for species in species_data])

        return new(filehandle, 
                   timestep_dim, v_spn, v_timestep, v_np, v_flux_incident, v_flux_reflected, v_force,
                   v_normal_pressure, v_shear_pressure, v_kinetic_energy_flux,
                   [surf_props.n_species, 1], [surf_props.n_elements, surf_props.n_species, 1],
                   [3, surf_props.n_elements, surf_props.n_species, 1],
                   [1], [1, 1], [1, 1, 1], [1, 1, 1, 1], [0.0],
                   skip_list)
    end

    @doc """
        NCDataHolderSurf(nc_filename, species_data, surf_props; global_attributes=Dict{Any,Any}())

    Construct a `NCDataHolderSurf` instance with an empty list of variable to skip.

    # Positional arguments

    # Positional arguments
    * `nc_filename`: filename to write output to
    * `species_data`: the vector of `Species` data for the species in the simulation
    * `surf_props`: the `SurfProps` instance which will be used for the computation and output of surface properties
    
    # Keyword arguments
    * `global_attributes`: dictionary of any additional attributes to write to the netCDF file as a global attribute
    """
    function NCDataHolderSurf(nc_filename, species_data, surf_props; global_attributes=Dict{Any,Any}())
        return NCDataHolderSurf(nc_filename, [], species_data, surf_props; global_attributes=global_attributes)
    end
end


"""
    write_netcdf_surf_props(ds, surf_props, timestep; sync_freq=0)
    
Write SurfProps to NetCDF file and synchronize file to disk if necessary.

# Positional arguments
* `ds`: the `NCDataHolderSurf` for the file to which the output will be written
* `phys_props`: the `SurfProps` instance containing the computed properties
* `timestep`: the simulation timestep

# Keyword arguments
* `sync_freq`: if larger than 0 and if the number of timesteps output is proportional to `sync_freq`,
    the data will be synchronized to disk. If set to 1, will sync data to disk at every timestep at which
    data is written to the file.
"""
function write_netcdf_surf_props(ds, surf_props, timestep; sync_freq=0)
    currtimesteps = ds.timestep_dim.dimlen + 1

    @inbounds ds.currtimesteps[1] = currtimesteps
    @inbounds ds.currtimesteps_1[2] = currtimesteps
    @inbounds ds.currtimesteps_1_1[3] = currtimesteps
    @inbounds ds.currtimesteps_1_1_1[4] = currtimesteps
    @inbounds ds.timestep[1] = timestep

    NetCDF.putvar(ds.v_timestep, ds.timestep, start=ds.currtimesteps)

    if !ds.skip_list.skip_number_of_particles
        NetCDF.putvar(ds.v_np, surf_props.np, start=ds.currtimesteps_1_1, count=ds.n_elements_n_species_1)
    end

    if !ds.skip_list.skip_fluxes
        NetCDF.putvar(ds.v_flux_incident, surf_props.flux_incident, start=ds.currtimesteps_1_1, count=ds.n_elements_n_species_1)
        NetCDF.putvar(ds.v_flux_reflected, surf_props.flux_reflected, start=ds.currtimesteps_1_1, count=ds.n_elements_n_species_1)
    end

    if !ds.skip_list.skip_force
        NetCDF.putvar(ds.v_force, surf_props.force, start=ds.currtimesteps_1_1_1, count=ds.n_v_n_elements_n_species_1)
    end

    if !ds.skip_list.skip_normal_pressure
        NetCDF.putvar(ds.v_normal_pressure, surf_props.normal_pressure, start=ds.currtimesteps_1_1, count=ds.n_elements_n_species_1)
    end

    if !ds.skip_list.skip_shear_pressure
        NetCDF.putvar(ds.v_shear_pressure, surf_props.shear_pressure, start=ds.currtimesteps_1_1_1, count=ds.n_v_n_elements_n_species_1)
    end

    if !ds.skip_list.skip_kinetic_energy_flux
        NetCDF.putvar(ds.v_kinetic_energy_flux, surf_props.kinetic_energy_flux, start=ds.currtimesteps_1_1, count=ds.n_elements_n_species_1)
    end

    if (sync_freq > 0) && (currtimesteps % sync_freq == 0)
        NetCDF.sync(ds.filehandle)
    end
end
