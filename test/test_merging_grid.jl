@testset "merging_grid" begin

    species_list::Vector{Species} = load_species_list("data/particles.toml", "Ar")

    seed = 1234
    Random.seed!(seed)
    rng::Xoshiro = Xoshiro(seed)

    Nx = 4
    Ny = 3
    Nz = 2
    mg = create_merging_grid(Nx, Ny, Nz, 1.0)

    phys_props::PhysProps = create_props(1, 1, [], Tref=1)

    phys_props.T[1,1] = species_list[1].mass / (2 * k_B)

    @test mg.Ntotal == 32  # Nx * Ny * Nz + 8
    @test mg.NyNz == 6  # Ny * Nz

    Merzbild.compute_velocity_extent!(1, 1, mg, phys_props, species_list)

    counter = 0
    for i in 1:Nx
        for j in 1:Ny
            for k in 1:Nz
                counter += 1
                v = [mg.extent_v_lower[1] + 1e-5 + (i-1) * mg.Δv[1],
                     mg.extent_v_lower[2] + 1e-5 + (j-1) * mg.Δv[2],
                     mg.extent_v_lower[3] + 1e-5 + (k-1) * mg.Δv[3]]
                index = Merzbild.compute_grid_index(mg, v)
                @test index == counter
            end
        end
    end
end