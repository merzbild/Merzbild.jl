using NCDatasets

function create_netcdf_phys_props(nc_filename, phys_props, species_data)
    ds = Dataset(nc_filename, "c")

    defDim(ds, "n_cells", phys_props.n_cells)
    defDim(ds, "n_species", phys_props.n_species)
    defDim(ds, "time", Inf)
    defDim(ds, "n_moments", phys_props.n_moments)

    v_spn = defVar(ds, "species_names", String, ("n_species",))
    v_spn[:] = [species.name for species in species_data]

    v_mompows = defVar(ds, "moment_powers", Int64, ("n_moments",))
    v_mompows[:] = phys_props.moment_powers

    v_ndens = defVar(ds, "ndens", Float64, ("n_cells", "n_species", "time"))
    v_vx = defVar(ds, "vx", Float64, ("n_cells", "n_species", "time"))
    v_vy = defVar(ds, "vy", Float64, ("n_cells", "n_species", "time"))
    v_vz = defVar(ds, "vz", Float64, ("n_cells", "n_species", "time"))
    v_T = defVar(ds, "T", Float64, ("n_cells", "n_species", "time"))
    v_time = defVar(ds, "timestep", Float64, ("time",))
    v_moments = defVar(ds, "moments", Float64, ("n_cells", "n_species", "n_moments", "time"))
    return ds
end

function write_netcdf_phys_props(ds, phys_props, timestep)
    currtimesteps = ds.dim["time"] + 1

    ds["timestep"][currtimesteps] = timestep
    for species in 1:phys_props.n_species
        for cell in 1:phys_props.n_cells
            ds["ndens"][cell, species, currtimesteps] = phys_props.n[cell, species]
            ds["vx"][cell, species, currtimesteps] = phys_props.v[cell, species, 1]
            ds["vy"][cell, species, currtimesteps] = phys_props.v[cell, species, 2]
            ds["vz"][cell, species, currtimesteps] = phys_props.v[cell, species, 3]
            ds["T"][cell, species, currtimesteps] = phys_props.T[cell, species, 1]
            # for mom_id in 1:phys_props.n_moments
            if (phys_props.n_moments > 0)
                ds["moments"][cell, species, :, currtimesteps] = phys_props.moments[cell, species, :]
            end
            # end
        end
    end
    # println(currtimesteps)
end