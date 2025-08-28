@testset "I/O for surf props" begin

    function set_surf_data!(sp, timestep)
        sp.np[1,1] = 1.0 * timestep
        sp.np[2,1] = 2.0 * timestep

        sp.flux_incident[1,1] = 2.0 * timestep
        sp.flux_incident[2,1] = 3.0 * timestep

        sp.flux_reflected[1,1] = 4.0 * timestep
        sp.flux_reflected[2,1] = 5.0 * timestep

        sp.force[1,1,1] = 10.0 * timestep
        sp.force[2,1,1] = 11.0 * timestep
        sp.force[3,1,1] = 12.0 * timestep

        sp.force[1,2,1] = 20.0 * timestep
        sp.force[2,2,1] = 21.0 * timestep
        sp.force[3,2,1] = 22.0 * timestep

        sp.normal_pressure[1,1] = 22.0 * timestep
        sp.normal_pressure[2,1] = 25.0 * timestep

        sp.shear_pressure[1,1,1] = 40.0 * timestep
        sp.shear_pressure[2,1,1] = 41.0 * timestep
        sp.shear_pressure[3,1,1] = 42.0 * timestep

        sp.shear_pressure[1,2,1] = 30.0 * timestep
        sp.shear_pressure[2,2,1] = 31.0 * timestep
        sp.shear_pressure[3,2,1] = 32.0 * timestep

        sp.kinetic_energy_flux[1,1] = 27.0 * timestep
        sp.kinetic_energy_flux[2,1] = 29.0 * timestep
    end

    function compare_readin_to_actual!(sp_from_file, sp, n_timesteps)
        for i in 0:n_timesteps
            set_surf_data!(sp, i)

            @test maximum(abs.(sp.np - sp_from_file["np"][:,:,i+1])) < 2*eps()
            @test maximum(abs.(sp.flux_incident - sp_from_file["flux_incident"][:,:,i+1])) < 2*eps()
            @test maximum(abs.(sp.flux_reflected - sp_from_file["flux_reflected"][:,:,i+1])) < 2*eps()
            @test maximum(abs.(sp.force - sp_from_file["force"][:,:,:,i+1])) < 2*eps()
            @test maximum(abs.(sp.normal_pressure - sp_from_file["normal_pressure"][:,:,i+1])) < 2*eps()
            @test maximum(abs.(sp.shear_pressure - sp_from_file["shear_pressure"][:,:,:,i+1])) < 2*eps()
            @test maximum(abs.(sp.kinetic_energy_flux - sp_from_file["kinetic_energy_flux"][:,:,i+1])) < 2*eps()
        end
    end

    # various basic I/O routines

    iosl = IOSkipListSurf()
    @test iosl.skip_number_of_particles == false
    @test iosl.skip_fluxes == false
    @test iosl.skip_force == false
    @test iosl.skip_normal_pressure == false
    @test iosl.skip_shear_pressure == false
    @test iosl.skip_kinetic_energy_flux == false


    names_skip_list = ["np", "force", "shear_pressure"]
    iosl2 = IOSkipListSurf(names_skip_list)
    @test iosl2.skip_number_of_particles == true
    @test iosl2.skip_fluxes == false
    @test iosl2.skip_force == true
    @test iosl2.skip_normal_pressure == false
    @test iosl2.skip_shear_pressure == true
    @test iosl2.skip_kinetic_energy_flux == false


    names_skip_list_3 = ["fluxes", "normal_pressure", "kinetic_energy_flux"]
    iosl3 = IOSkipListSurf(names_skip_list_3)
    @test iosl3.skip_number_of_particles == false
    @test iosl3.skip_fluxes == true
    @test iosl3.skip_force == false
    @test iosl3.skip_normal_pressure == true
    @test iosl3.skip_shear_pressure == false
    @test iosl3.skip_kinetic_energy_flux == true


    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")

    grid = Grid1DUniform(4.0, 8)
    pia = ParticleIndexerArray(grid.n_cells, 1)
    surf_props = SurfProps(pia, grid)

    sol_path = joinpath(@__DIR__, "data", "tmp_no_skiplist.nc")
    ds = NCDataHolderSurf(sol_path, [], species_data, surf_props)

    set_surf_data!(surf_props, 0)
    write_netcdf(ds, surf_props, 0)

    set_surf_data!(surf_props, 1)
    write_netcdf(ds, surf_props, 1, sync_freq=1)

    set_surf_data!(surf_props, 2)
    write_netcdf(ds, surf_props, 2)
    close_netcdf(ds)

    props_read =  NCDataset(sol_path, "r")

    varnames = []
    for (varname,var) in props_read
        push!(varnames, varname)
    end

    for name in ["np", "flux_incident", "flux_reflected", "force", "normal_pressure", "shear_pressure", "kinetic_energy_flux"]
        @test name in varnames
    end

    @test size(props_read["np"]) == (2, 1, 3)
    @test size(props_read["flux_incident"]) == (2, 1, 3)
    @test size(props_read["flux_reflected"]) == (2, 1, 3)
    @test size(props_read["force"]) == (3, 2, 1, 3)
    @test size(props_read["normal_pressure"]) == (2, 1, 3)
    @test size(props_read["shear_pressure"]) == (3, 2, 1, 3)
    @test size(props_read["kinetic_energy_flux"]) == (2, 1, 3)

    compare_readin_to_actual!(props_read, surf_props, 2)
    close(props_read)
    rm(sol_path)

    
    sol_path = joinpath(@__DIR__, "data", "tmp_skiplist.nc")
    ds = NCDataHolderSurf(sol_path, names_skip_list, species_data, surf_props)

    set_surf_data!(surf_props, 0)
    write_netcdf(ds, surf_props, 0)
    close_netcdf(ds)

    props_read =  NCDataset(sol_path, "r")

    varnames = []
    for (varname,var) in props_read
        push!(varnames, varname)
    end
    
    for name in names_skip_list
        @test !(name in varnames)
    end

    for name in ["flux_incident", "flux_reflected", "normal_pressure", "kinetic_energy_flux"]
        @test name in varnames
    end

    close(props_read)
    rm(sol_path)
end