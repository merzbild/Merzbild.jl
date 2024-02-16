using NCDatasets

function create_netcdf_phys_props(nc_filename, phys_props, species_data)
    ds = Dataset(nc_filename, "c")

    defDim(ds, "vel_components", 3)
    defDim(ds, "n_cells", phys_props.n_cells)
    defDim(ds, "n_species", phys_props.n_species)
    defDim(ds, "time", Inf)
    defDim(ds, "n_moments", phys_props.n_moments)

    v_spn = defVar(ds, "species_names", String, ("n_species",))
    v_spn[:] = [species.name for species in species_data]

    v_mompows = defVar(ds, "moment_powers", Int64, ("n_moments",))
    v_mompows[:] = phys_props.moment_powers

    v_lpa = defVar(ds, "length_particle_array", Float64, ("n_species", "time"))
    v_np = defVar(ds, "np", Float64, ("n_cells", "n_species", "time"))
    v_ndens = defVar(ds, "ndens", Float64, ("n_cells", "n_species", "time"))
    v_v = defVar(ds, "v", Float64, ("vel_components", "n_cells", "n_species", "time"))
    v_T = defVar(ds, "T", Float64, ("n_cells", "n_species", "time"))
    v_time = defVar(ds, "timestep", Float64, ("time",))
    v_moments = defVar(ds, "moments", Float64, ("n_moments", "n_cells", "n_species", "time"))
    return ds
end

function write_netcdf_phys_props(ds, phys_props, timestep)
    currtimesteps = ds.dim["time"] + 1

    ds["timestep"][currtimesteps] = timestep
    for species in 1:phys_props.n_species
        for cell in 1:phys_props.n_cells
            ds["length_particle_array"][species, currtimesteps] = phys_props.lpa[species]
            ds["np"][cell, species, currtimesteps] = phys_props.np[cell, species]
            ds["ndens"][cell, species, currtimesteps] = phys_props.n[cell, species]
            ds["v"][:, cell, species, currtimesteps] = phys_props.v[:,cell, species]
            ds["T"][cell, species, currtimesteps] = phys_props.T[cell, species, 1]

            if (phys_props.n_moments > 0)
                ds["moments"][:, cell, species, currtimesteps] = phys_props.moments[:,cell, species]
            end
            # end
        end
    end
    # println(currtimesteps)
end