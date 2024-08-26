@testset "species data" begin
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_data(particles_data_path, "Ar")
    @test length(species_data) == 1
    @test species_data[1].name == "Ar"
    @test abs(species_data[1].mass - 6.63e-26)/6.63e-26 < 1e-8
    @test species_data[1].charge == 0
    @test abs(species_data[1].charge_div_mass) == 0

    species_data = load_species_data(particles_data_path, ["Ar", "e-"])
    @test length(species_data) == 2
    @test species_data[1].name == "Ar"
    @test abs(species_data[1].mass - 6.63e-26)/6.63e-26 < 1e-10
    @test species_data[1].charge == 0
    @test abs(species_data[1].charge_div_mass) == 0

    @test species_data[2].name == "e-"
    @test abs(species_data[2].mass - 9.109383632e-31)/9.109383632e-31 < 1e-10
    @test species_data[2].charge == -1
    @test abs(species_data[2].charge_div_mass - (-1) * Merzbild.q_e / species_data[2].mass) < eps()
end