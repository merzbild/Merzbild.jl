@testset "I/O" begin
    # various basic I/O routines

    iosl = IOSkipList()
    @test iosl.skip_length_particle_array == false
    @test iosl.skip_number_of_particles == false
    @test iosl.skip_number_density == false
    @test iosl.skip_velocity == false
    @test iosl.skip_temperature == false

    names_skip_list = ["length_particle_array", "ndens", "v"]
    iosl = IOSkipList(names_skip_list)
    @test iosl.skip_length_particle_array == true
    @test iosl.skip_number_of_particles == false
    @test iosl.skip_number_density == true
    @test iosl.skip_velocity == true
    @test iosl.skip_temperature == false

    names_skip_list = ["T", "np"]
    iosl = IOSkipList(names_skip_list)
    @test iosl.skip_length_particle_array == false
    @test iosl.skip_number_of_particles == true
    @test iosl.skip_number_density == false
    @test iosl.skip_velocity == false
    @test iosl.skip_temperature == true

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")
    phys_props::PhysProps = PhysProps(4, 2)

    @test phys_props.ndens_not_Np == false

    # netcdf4 format
    sol_path = joinpath(@__DIR__, "data", "tmp_skiplist_nc4.nc")
    ds = NCDataHolder(sol_path, names_skip_list, species_data, phys_props; mode=NC_NETCDF4)

    @test ds.ndens_not_Np == false

    phys_props.T .= -100.0
    phys_props.n .= -100.0
    phys_props.np .= -100.0

    phys_props.lpa .= 30.0
    phys_props.v .= 30.0

    write_netcdf(ds, phys_props, 0)
    write_netcdf(ds, phys_props, 1)
    write_netcdf(ds, phys_props, 2)
    close_netcdf(ds)

    props_read =  NCDataset(sol_path, "r")

    varnames = []
    for (varname,var) in props_read
        push!(varnames, varname)
    end
    
    for name in names_skip_list
        @test !(name in varnames)
    end

    for name in ["v", "ndens", "length_particle_array"]
        @test name in varnames
    end

    @test size(props_read["v"]) == (3, 4, 2, 3)
    @test size(props_read["ndens"]) == (4, 2, 3)
    @test size(props_read["length_particle_array"]) == (2, 3)

    close(props_read)
    rm(sol_path)

    # netcdf classic 64-bit format
    sol_path = joinpath(@__DIR__, "data", "tmp_skiplist_nc_classic.nc")
    ds = NCDataHolder(sol_path, names_skip_list, species_data, phys_props; mode=NC_64BIT_OFFSET)

    @test ds.ndens_not_Np == false

    phys_props.T .= -100.0
    phys_props.n .= -100.0
    phys_props.np .= -100.0

    phys_props.lpa .= 30.0
    phys_props.v .= 30.0

    write_netcdf(ds, phys_props, 0)
    write_netcdf(ds, phys_props, 1)
    write_netcdf(ds, phys_props, 2)
    close_netcdf(ds)

    props_read =  NCDataset(sol_path, "r")

    varnames = []
    for (varname,var) in props_read
        push!(varnames, varname)
    end
    
    for name in names_skip_list
        @test !(name in varnames)
    end

    for name in ["v", "ndens", "length_particle_array"]
        @test name in varnames
    end

    @test size(props_read["v"]) == (3, 4, 2, 3)
    @test size(props_read["ndens"]) == (4, 2, 3)
    @test size(props_read["length_particle_array"]) == (2, 3)

    close(props_read)
    rm(sol_path)

    # now we test that if we start writing number density we can't write a phys_props
    # with # of physical particles

    phys_props_ndens::PhysProps = PhysProps(4, 2; ndens_not_Np=true)

    @test phys_props_ndens.ndens_not_Np == true

    sol_path = joinpath(@__DIR__, "data", "tmp_ndens.nc")
    ds = NCDataHolder(sol_path, [], species_data, phys_props_ndens)

    @test ds.ndens_not_Np == true

    @test_throws ErrorException write_netcdf(ds, phys_props, 0)
    close_netcdf(ds)
    rm(sol_path)

    # test moment I/O
    sol_path = joinpath(@__DIR__, "data", "tmp_moments.nc")

    ds = NCDataHolderMoments(sol_path, species_data,4, 2, [4,6,8])

    moment_vals = zeros((3,4,2))

    for k in 1:2
        for j in 1:4
            for i in 1:3
                moment_vals[i,j,k] = k + j*3 + i*20.0
            end
        end
    end

    moment_vals_2 = moment_vals * 2.0

    write_netcdf(ds, moment_vals, 1)
    write_netcdf(ds, moment_vals_2, 2; sync_freq=1)

    close_netcdf(ds)

    props_read =  NCDataset(sol_path, "r")
    @test size(props_read["timestep"]) == (2,)
    @test size(props_read["moment_powers"]) == (3,)
    @test size(props_read["moments"]) == (3,4,2,2)

    @test props_read["timestep"][:] == [1,2]
    @test props_read["moment_powers"][:] == [4,6,8]
    @test maximum(abs.(props_read["moments"][:,:,:,1] .- moment_vals)) < 2*eps()
    @test maximum(abs.(props_read["moments"][:,:,:,2] .- moment_vals_2)) < 2*eps()

    close(props_read)

    rm(sol_path)
end