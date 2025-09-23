@testset "I/O for flux_props" begin

    function set_flux_data!(fp, timestep)
        fp.kinetic_energy_flux[:,1,1] = [1.0, 2.0, 3.0] * timestep
        fp.kinetic_energy_flux[:,2,1] = [-4.0, -5.0, -6.0]  * timestep

        fp.diagonal_momentum_flux[:,1,1] = [21.0, 12.0, -3.0] * timestep
        fp.diagonal_momentum_flux[:,2,1] = [34.0, -1.0, 8.0]  * timestep

        fp.off_diagonal_momentum_flux[:,1,1] = [-1.0, 0.0, 2.0] * timestep
        fp.off_diagonal_momentum_flux[:,2,1] = [-5.0, 4.0, 8.0]  * timestep
    end

    function compare_readin_to_actual!(fp_from_file, fp, n_timesteps)
        for i in 0:n_timesteps
            set_surf_data!(fp, i)

            @test maximum(abs.(sp.kinetic_energy_flux - fp_from_file["kinetic_energy_flux"][:,:,i+1])) < 2*eps()
            @test maximum(abs.(sp.diagonal_momentum_flux - fp_from_file["diagonal_momentum_flux"][:,:,i+1])) < 2*eps()
            @test maximum(abs.(sp.off_diagonal_momentum_flux - fp_from_file["off_diagonal_momentum_flux"][:,:,i+1])) < 2*eps()
        end
    end

    # various basic I/O routines

    iosl = IOSkipListFlux()
    @test iosl.skip_kinetic_energy_flux == false
    @test iosl.skip_diagonal_momentum_flux == false
    @test iosl.skip_off_diagonal_momentum_flux == false


    names_skip_list = ["kinetic_energy_flux", "diagonal_momentum_flux", "off_diagonal_momentum_flux"]
    iosl2 = IOSkipListFlux(names_skip_list)
    @test iosl2.skip_kinetic_energy_flux == true
    @test iosl2.skip_diagonal_momentum_flux == true
    @test iosl2.skip_off_diagonal_momentum_flux == true

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")

    grid = Grid1DUniform(4.0, 8)
    pia = ParticleIndexerArray(grid.n_cells, 1)
    flux_props = FluxProps(pia)

    sol_path = joinpath(@__DIR__, "data", "tmp_no_skiplist.nc")
    ds = NCDataHolderFlux(sol_path, [], species_data, flux_props)

    set_flux_data!(flux_props, 0)
    write_netcdf(ds, flux_props, 0)

    set_flux_data!(flux_props, 1)
    write_netcdf(ds, flux_props, 1, sync_freq=1)

    set_flux_data!(flux_props, 2)
    write_netcdf(ds, flux_props, 2)
    close_netcdf(ds)

    props_read =  NCDataset(sol_path, "r")

    varnames = []
    for (varname,var) in props_read
        push!(varnames, varname)
    end

    for name in ["kinetic_energy_flux", "diagonal_momentum_flux", "off_diagonal_momentum_flux"]
        @test name in varnames
    end

    @test size(props_read["kinetic_energy_flux"]) == (3, 2, 1, 3)
    @test size(props_read["diagonal_momentum_flux"]) == (3, 2, 1, 3)
    @test size(props_read["off_diagonal_momentum_flux"]) == (3, 2, 1, 3)

    compare_readin_to_actual!(props_read, surf_props, 2)
    close(props_read)
    rm(sol_path)

    
    sol_path = joinpath(@__DIR__, "data", "tmp_skiplist.nc")
    ds = NCDataHolderFlux(sol_path, names_skip_list, species_data, flux_props)

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

    for name in []
        @test name in varnames
    end

    close(props_read)
    rm(sol_path)
end