@testset "electron_neutral_data" begin
    # test loading of data from XML file and linear interpolation thereof

    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_list = load_species_list(particles_data_path, ["He", "e-"])

    interaction_data_path = joinpath(@__DIR__, "..", "data", "vhs.toml")
    interaction_data = load_interaction_data(interaction_data_path, species_list, 1e-10, 1.0, 273.0)

    e_n_data_path = joinpath(@__DIR__, "..", "data", "test_neutral_electron_data.xml")
    n_e_interactions = load_electron_neutral_interactions(species_list, e_n_data_path,
                                                          Dict("He" => "LinearDB"),
                                                          Dict("He" => ScatteringIsotropic),
                                                          Dict("He" => ElectronEnergySplitEqual))

    @test length(n_e_interactions.elastic[1].data.E) == 4
    @test length(n_e_interactions.elastic[1].data.sigma) == 4
    @test length(n_e_interactions.ionization[1].data.E) == 3
    @test length(n_e_interactions.ionization[1].data.sigma) == 3

    @test n_e_interactions.elastic[1].data.E[1] == 0.0
    @test n_e_interactions.elastic[1].data.E[end] == 1e3

    @test n_e_interactions.ionization[1].data.sigma[1] == 0.0
    @test n_e_interactions.ionization[1].data.E[end] == 1e2

    n_e_cs = create_computed_crosssections(n_e_interactions)
    @test length(n_e_cs[1].prob_vec) == 2
    @test length(n_e_cs[1].cdf_prob_vec) == 2

    for E_coll in [1.0, 20.0]
        g = sqrt(2 * E_coll * Merzbild.eV_J / interaction_data[1,2].m_r)
        Merzbild.compute_cross_sections!(n_e_cs, interaction_data[1,2], g, n_e_interactions, 1)

        # linear interpolation of linear data is exact
        @test abs(n_e_cs[1].cs_elastic * 1e20 - (4 + E_coll / 100.0)) / 1e20 < eps()
        @test n_e_cs[1].cs_total == n_e_cs[1].cs_elastic
        @test n_e_cs[1].cs_ionization == 0.0
        @test n_e_cs[1].prob_vec[1] == 1.0
        @test n_e_cs[1].prob_vec[2] == 0.0
        @test n_e_cs[1].cdf_prob_vec[1] == 0.0
        @test n_e_cs[1].cdf_prob_vec[2] == 1.0
    end

    for E_coll in [50.0, 75.0]
        # now we have a non-zero ionization cross-section
        g = sqrt(2 * E_coll * Merzbild.eV_J / interaction_data[1,2].m_r)
        Merzbild.compute_cross_sections!(n_e_cs, interaction_data[1,2], g, n_e_interactions, 1)
        @test abs(n_e_cs[1].cs_elastic * 1e20 - (4 + E_coll / 100.0)) / 1e20 < eps()

        @test abs(n_e_cs[1].cs_total - (n_e_cs[1].cs_elastic + n_e_cs[1].cs_ionization)) < eps()
        @test abs(n_e_cs[1].cs_ionization * 1e25 - (1 * (100.0 - E_coll) / 70.0 + 10 * (E_coll - 30) / 70)) / 1e25 < eps()
        @test n_e_cs[1].prob_vec[1] < 1.0
        @test n_e_cs[1].prob_vec[2] > 0.0
        @test abs(sum(n_e_cs[1].prob_vec) - 1.0) < eps()

        @test n_e_cs[1].cdf_prob_vec[1] == 0.0
        @test abs(n_e_cs[1].cdf_prob_vec[2] - n_e_cs[1].cs_elastic / n_e_cs[1].cs_total) < eps()
    end

    for E_coll in [203.5]
        # now we're out of bounds for the ionization cross-section and the value should be the last one in the data set
        # (continued by a constant)
        g = sqrt(2 * E_coll * Merzbild.eV_J / interaction_data[1,2].m_r)
        Merzbild.compute_cross_sections!(n_e_cs, interaction_data[1,2], g, n_e_interactions, 1)
        @test abs(n_e_cs[1].cs_elastic * 1e20 - (4 + E_coll / 100.0)) / 1e20 < eps()


        @test abs(n_e_cs[1].cs_total - (n_e_cs[1].cs_elastic + n_e_cs[1].cs_ionization)) < eps()
        @test abs(n_e_cs[1].cs_ionization * 1e25 - n_e_interactions.ionization[1].data.sigma[end] * 1e25) / 1e25 < eps()

        @test n_e_cs[1].prob_vec[1] < 1.0
        @test n_e_cs[1].prob_vec[2] > 0.0
        @test abs(sum(n_e_cs[1].prob_vec) - 1.0) < eps()
    end

    # 200.0
end