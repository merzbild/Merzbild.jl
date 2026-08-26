@testset "I/O for field_props" begin

    function set_field_data!(fp, timestep)
        for i in 1:fp.n_nodes
            fp.charge_density[i] = 1e-8 * i * (timestep + 1)
            fp.potential[i] = -2.0 * i + 0.5 * timestep
            fp.electric_field[i] = 3.0 * i * timestep - 1.0
        end
    end

    function compare_readin_to_actual!(fp_from_file, fp, n_timesteps)
        for i in 0:n_timesteps
            set_field_data!(fp, i)

            @test maximum(abs.(fp.charge_density - fp_from_file["charge_density"][:,i+1])) < 2*eps()
            @test maximum(abs.(fp.potential - fp_from_file["potential"][:,i+1])) < 2*eps()
            @test maximum(abs.(fp.electric_field - fp_from_file["electric_field"][:,i+1])) < 2*eps()
        end
    end

    iosl = IOSkipListField()
    @test iosl.skip_charge_density == false
    @test iosl.skip_potential == false
    @test iosl.skip_electric_field == false

    iosl2 = IOSkipListField(["rho", "phi", "E"])
    @test iosl2.skip_charge_density == true
    @test iosl2.skip_potential == true
    @test iosl2.skip_electric_field == true

    names_skip_list = ["charge_density", "potential", "electric_field"]
    iosl3 = IOSkipListField(names_skip_list)
    @test iosl3.skip_charge_density == true
    @test iosl3.skip_potential == true
    @test iosl3.skip_electric_field == true

    iosl4 = IOSkipListField(["potential"])
    @test iosl4.skip_charge_density == false
    @test iosl4.skip_potential == true
    @test iosl4.skip_electric_field == false

    grid = Grid1DUniform(4.0, 8)
    field_props = ElectrostaticFieldProps(grid)
    @test field_props.n_nodes == grid.n_cells + 1

    for mode in [NC_64BIT_OFFSET, NC_NETCDF4]
        sol_path = joinpath(@__DIR__, "data", "tmp_field_no_skiplist.nc")
        ds = NCDataHolderField(sol_path, field_props;
                               global_attributes=Dict{Any,Any}("test_attribute" => "test_value"), mode=mode)

        set_field_data!(field_props, 0)
        write_netcdf(ds, field_props, 0)

        set_field_data!(field_props, 1)
        write_netcdf(ds, field_props, 1, sync_freq=1)

        set_field_data!(field_props, 2)
        write_netcdf(ds, field_props, 2)
        close_netcdf(ds)

        props_read = NCDataset(sol_path, "r")

        varnames = []
        for (varname, var) in props_read
            push!(varnames, varname)
        end

        for name in ["timestep", "charge_density", "potential", "electric_field"]
            @test name in varnames
        end

        @test props_read.attrib["test_attribute"] == "test_value"
        @test props_read.dim["n_nodes"] == grid.n_cells + 1

        @test length(props_read["timestep"]) == 3
        @test size(props_read["charge_density"]) == (grid.n_cells + 1, 3)
        @test size(props_read["potential"]) == (grid.n_cells + 1, 3)
        @test size(props_read["electric_field"]) == (grid.n_cells + 1, 3)

        @test maximum(abs.(props_read["timestep"][:] - [0.0, 1.0, 2.0])) < 2*eps()

        compare_readin_to_actual!(props_read, field_props, 2)
        close(props_read)
        rm(sol_path)
    end

    sol_path = joinpath(@__DIR__, "data", "tmp_field_skiplist.nc")
    ds = NCDataHolderField(sol_path, ["rho", "E"], field_props)

    set_field_data!(field_props, 0)
    write_netcdf(ds, field_props, 0)
    close_netcdf(ds)

    props_read = NCDataset(sol_path, "r")

    varnames = []
    for (varname, var) in props_read
        push!(varnames, varname)
    end

    for name in ["charge_density", "electric_field"]
        @test !(name in varnames)
    end

    for name in ["timestep", "potential"]
        @test name in varnames
    end

    @test maximum(abs.(field_props.potential - props_read["potential"][:,1])) < 2*eps()

    close(props_read)
    rm(sol_path)
end
