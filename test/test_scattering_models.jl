@testset "elastic scattering models" begin
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    vhs_data_path = joinpath(@__DIR__, "..", "data", "vhs.toml")
    vss_data_path = joinpath(@__DIR__, "..", "data", "vss.toml")

    species_data = load_species_data(particles_data_path, ["Ar", "He"])

    # parsing of the model names in the interaction data files
    @test Merzbild.parse_scattering_model("VHS") == ScatteringVHS
    @test Merzbild.parse_scattering_model("vhs") == ScatteringVHS
    @test Merzbild.parse_scattering_model("VSS") == ScatteringVSS
    @test Merzbild.parse_scattering_model("vss") == ScatteringVSS
    @test_throws ArgumentError Merzbild.parse_scattering_model("Maxwell")

    vhs_data = load_interaction_data(vhs_data_path, species_data)
    vss_data = load_interaction_data(vss_data_path, species_data)

    for i in 1:2
        for k in 1:2
            # the model is stored on the interaction, the array stays concretely typed
            @test vhs_data[i,k].model == ScatteringVHS
            @test vss_data[i,k].model == ScatteringVSS

            # models with isotropic scattering do not use the VSS exponent
            @test vhs_data[i,k].vss_alpha == 1.0
            @test vhs_data[i,k].vss_inv_alpha == 1.0

            @test vss_data[i,k].vss_alpha > 1.0
            @test abs(vss_data[i,k].vss_inv_alpha - 1.0 / vss_data[i,k].vss_alpha) < eps()

            # the VSS model uses the same total cross-section as the VHS model
            @test abs(vss_data[i,k].vhs_o - vhs_data[i,k].vhs_o) < eps()
            @test abs(vss_data[i,k].vhs_exp - vhs_data[i,k].vhs_exp) < eps()
            @test abs(vss_data[i,k].vhs_factor - vhs_data[i,k].vhs_factor) < eps()

            # but a different reference viscosity
            vss_factor = Merzbild.compute_vss_mu_ref_factor(vss_data[i,k].vss_alpha)
            @test vss_factor < 1.0
            @test abs(vss_data[i,k].vhs_muref - vhs_data[i,k].vhs_muref * vss_factor) < eps()

            # the asymmetry of the interaction array is preserved
            @test vss_data[i,k].model == vss_data[k,i].model
            @test abs(vss_data[i,k].vss_alpha - vss_data[k,i].vss_alpha) < eps()
            @test abs(vss_data[i,k].μ1 - vss_data[k,i].μ2) < eps()
        end
    end

    @test Merzbild.compute_vss_mu_ref_factor(1.0) == 1.0

    # the interaction constructors
    m_Ar = species_data[1].mass
    m_He = species_data[2].mass

    @test Interaction(m_Ar, m_He, 3.25e-10, 0.735, 273.0) == Interaction(VHS(), m_Ar, m_He, 3.25e-10, 0.735, 273.0)
    @test Interaction(VHS(), m_Ar, m_He, 3.25e-10, 0.735, 273.0) == vhs_data[1,2]
    @test Interaction(VSS(), m_Ar, m_He, 3.25e-10, 0.735, 273.0, 1.33) == vss_data[1,2]

    # a VSS interaction with alpha = 1 has the same parameters as the VHS one
    vss_alpha1 = Interaction(VSS(), m_Ar, m_He, 3.25e-10, 0.735, 273.0, 1.0)
    @test abs(vss_alpha1.vhs_muref - vhs_data[1,2].vhs_muref) < eps()
    @test abs(vss_alpha1.vhs_factor - vhs_data[1,2].vhs_factor) < eps()

    # a hard sphere gas is the VHS model with omega = 0.5: a constant cross-section
    hard_sphere = Interaction(VHS(), m_Ar, m_He, 3.25e-10, 0.5, 273.0)
    @test hard_sphere.vhs_exp == 0.0
    @test abs(hard_sphere.vhs_factor - π * hard_sphere.vhs_d^2) < eps()

    # cross-sections
    for g in [1.0, 100.0, 1234.5, 1e5]
        @test Merzbild.sigma(VHS(), vhs_data[1,2], g) == Merzbild.sigma_vhs(vhs_data[1,2], g)
        @test Merzbild.sigma(VSS(), vss_data[1,2], g) == Merzbild.sigma_vhs(vhs_data[1,2], g)

        # a hard sphere gas has a cross-section independent of the relative velocity
        @test Merzbild.sigma(VHS(), hard_sphere, g) == π * hard_sphere.vhs_d^2
    end

    # missing VSS data in the interaction file
    bad_vss_path = joinpath(@__DIR__, "data", "tmp_bad_vss.toml")
    open(bad_vss_path, "w") do io
        write(io, "[\"Ar,Ar\"]\nmodel = \"VSS\"\nvhs_d = 4.11e-10\nvhs_o = 0.81\nvhs_Tref = 273.0\n")
    end
    @test_throws KeyError load_interaction_data(bad_vss_path, load_species_data(particles_data_path, ["Ar"]))
    rm(bad_vss_path)

    # unknown model name in the interaction file
    bad_model_path = joinpath(@__DIR__, "data", "tmp_bad_model.toml")
    open(bad_model_path, "w") do io
        write(io, "[\"Ar,Ar\"]\nmodel = \"Maxwell\"\nvhs_d = 4.11e-10\nvhs_o = 0.81\nvhs_Tref = 273.0\n")
    end
    @test_throws ArgumentError load_interaction_data(bad_model_path, load_species_data(particles_data_path, ["Ar"]))
    rm(bad_model_path)
end

@testset "collisions with per-pair scattering models" begin
    particles_data_path = joinpath(@__DIR__, "..", "data", "particles.toml")
    species_data = load_species_data(particles_data_path, ["Ar", "He"])
    n_species = 2

    n_particles = [400, 4000]
    T0_list = [3000.0, 360.0]
    Fnum = 5e12

    n_Ar = Fnum * n_particles[1]
    n_He = Fnum * n_particles[2]
    T_eq = (n_Ar * T0_list[1] + n_He * T0_list[2]) / (n_Ar + n_He)

    function run_collisions(interaction_file, model, n_t)
        rng = StableRNG(1234)

        interaction_data = load_interaction_data(joinpath(@__DIR__, "..", "data", interaction_file), species_data)

        particles = [ParticleVector(n_particles[1]), ParticleVector(n_particles[2])]
        pia = ParticleIndexerArray([0, 0])
        for s in 1:n_species
            sample_particles_equal_weight!(rng, particles[s], pia, 1, s, n_particles[s], species_data[s].mass,
                                           T0_list[s], Fnum, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
        end

        phys_props = PhysProps(1, n_species)
        collision_factors = create_collision_factors_array(n_species)
        collision_data = CollisionData()
        estimate_sigma_g_w_max!(collision_factors, interaction_data, species_data, T0_list, Fnum)

        Δt = 2.5e-3
        V = 1.0

        for _ in 1:n_t
            for s2 in 1:n_species
                for s1 in s2:n_species
                    if (s1 == s2)
                        if model === nothing
                            ntc_equal_weight!(rng, collision_factors[s1,s1,1], collision_data, interaction_data,
                                              particles[s1], pia, 1, s1, Δt, V)
                        else
                            ntc_equal_weight!(rng, model, collision_factors[s1,s1,1], collision_data, interaction_data,
                                              particles[s1], pia, 1, s1, Δt, V)
                        end
                    else
                        if model === nothing
                            ntc_equal_weight!(rng, collision_factors[s1,s2,1], collision_data, interaction_data,
                                              particles[s1], particles[s2], pia, 1, s1, s2, Δt, V)
                        else
                            ntc_equal_weight!(rng, model, collision_factors[s1,s2,1], collision_data, interaction_data,
                                              particles[s1], particles[s2], pia, 1, s1, s2, Δt, V)
                        end
                    end
                end
            end
        end

        compute_props!(particles, pia, species_data, phys_props)
        return phys_props
    end

    # bulk velocity and total energy of the mixture, both conserved by elastic collisions
    function mixture_momentum(phys_props, component)
        return sum(phys_props.n[1,s] * species_data[s].mass * phys_props.v[component,1,s] for s in 1:n_species)
    end

    function mixture_energy(phys_props)
        return sum(phys_props.n[1,s] * (1.5 * k_B * phys_props.T[1,s] +
                                        0.5 * species_data[s].mass * sum(phys_props.v[:,1,s].^2)) for s in 1:n_species)
    end

    # the model stored in the Interaction instance is picked up by the drivers, and passing
    # the model tag explicitly gives exactly the same result
    for (interaction_file, model) in [("vhs.toml", VHS()), ("vss.toml", VSS())]
        props_from_data = run_collisions(interaction_file, nothing, 20)
        props_from_tag = run_collisions(interaction_file, model, 20)

        @test props_from_data.T == props_from_tag.T
        @test props_from_data.v == props_from_tag.v
        @test props_from_data.np == props_from_tag.np
    end

    # all models conserve mass, momentum and energy and relax the two-species mixture
    # towards a common temperature
    props_initial = run_collisions("vhs.toml", nothing, 0)
    mass_total = sum(props_initial.n[1,s] * species_data[s].mass for s in 1:n_species)
    ΔT_initial = abs(props_initial.T[1,1] - props_initial.T[1,2])

    for interaction_file in ["vhs.toml", "vss.toml"]
        phys_props = run_collisions(interaction_file, nothing, 800)

        @test phys_props.np[1,1] == n_particles[1]
        @test phys_props.np[1,2] == n_particles[2]
        @test abs(phys_props.n[1,1] - n_Ar) / n_Ar < eps()
        @test abs(phys_props.n[1,2] - n_He) / n_He < eps()

        for i in 1:3
            Δv = (mixture_momentum(phys_props, i) - mixture_momentum(props_initial, i)) / mass_total
            @test abs(Δv) < 1e-6
        end

        E_initial = mixture_energy(props_initial)
        @test abs(mixture_energy(phys_props) - E_initial) / E_initial < 1e-10

        # the species temperatures have moved towards each other
        @test abs(phys_props.T[1,1] - phys_props.T[1,2]) < 0.2 * ΔT_initial
        @test min(phys_props.T[1,1], phys_props.T[1,2]) > 0.8 * T_eq
        @test max(phys_props.T[1,1], phys_props.T[1,2]) < 1.2 * T_eq
    end
end
