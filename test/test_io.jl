@testset "I/O" begin
    # various basic I/O routines

    iosl = IOSkipList()
    @test iosl.skip_length_particle_array == false
    @test iosl.skip_moments == false
    @test iosl.skip_number_of_particles == false
    @test iosl.skip_number_density == false
    @test iosl.skip_velocity == false
    @test iosl.skip_temperature == false


    names_skip_list = ["length_particle_array", "ndens", "v"]
    iosl = IOSkipList(names_skip_list)
    @test iosl.skip_length_particle_array == true
    @test iosl.skip_moments == false
    @test iosl.skip_number_of_particles == false
    @test iosl.skip_number_density == true
    @test iosl.skip_velocity == true
    @test iosl.skip_temperature == false

    names_skip_list = ["T", "np"]
    iosl = IOSkipList(names_skip_list)
    @test iosl.skip_length_particle_array == false
    @test iosl.skip_moments == false
    @test iosl.skip_number_of_particles == true
    @test iosl.skip_number_density == false
    @test iosl.skip_velocity == false
    @test iosl.skip_temperature == true

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")
    phys_props::PhysProps = PhysProps(4, 2, [], Tref=1)

    @test phys_props.ndens_not_Np == false

    sol_path = joinpath(@__DIR__, "data", "tmp_skiplist.nc")
    ds = NCDataHolder(sol_path, names_skip_list, species_data, phys_props)

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

    phys_props_ndens::PhysProps = PhysProps(4, 2, [], Tref=1; ndens_not_Np=true)

    @test phys_props_ndens.ndens_not_Np == true

    sol_path = joinpath(@__DIR__, "data", "tmp_ndens.nc")
    ds = NCDataHolder(sol_path, [], species_data, phys_props_ndens)

    @test ds.ndens_not_Np == true

    @test_throws ErrorException write_netcdf(ds, phys_props, 0)
    close_netcdf(ds)
    rm(sol_path)
end