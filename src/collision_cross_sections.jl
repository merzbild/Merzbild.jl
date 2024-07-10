using XML

struct DataMissingException <: Exception
    msg::String
end

struct TabulatedCSData
    n_vals::Int32
    E::Vector{Float64}
    sigma::Vector{Float64}
    ΔE::Float64  # how much energy is lost / threshold energy, eV
end

@enum ScatteringLaw ScatteringIsotropic=1 ScatteringOkhrimovskyy=2
@enum ElectronEnergySplit ElectronEnergySplitEqual=1 ElectronEnergySplitZeroE=2 # Tanh law
struct Ionization
    data::TabulatedCSData
    scattering::ScatteringLaw
    split::ElectronEnergySplit
end
struct ElasticScattering
    data::TabulatedCSData
    scattering::ScatteringLaw
end
struct ExcitationSink
    n_reactions::Int32
    data::Vector{TabulatedCSData}
    scattering::ScatteringLaw
end

struct ElectronNeutralInteractions
   n_neutrals::Int32
   neutral_indexer::Vector{Int32}  # neutral_indexer[species_index] = index of the species in the sub-list of neutral species
   elastic::Vector{ElasticScattering}
   ionization::Vector{Ionization}
   excitation_sink::Vector{ExcitationSink}
end

mutable struct ComputedCrossSections
    n_excitations::Int32
    cs_total::Float64
    cs_elastic::Float64
    cs_ionization::Float64
    cs_excitation::Vector{Float64}
    prob_vec::Vector{Float64}
    cdf_prob_vec::Vector{Float64}
end

function sigma_vhs(interaction, g)
    return interaction.vhs_factor * g^(1.0 - 2 * interaction.vhs_o)
end

function binary_search(x, val)
    low = 1
    high = length(x)
    mid = 0

    if val > x[high]
        return 0
    elseif val < x[low]
        return -1
    end
 
    while low <= high
 
        mid = (high + low) ÷ 2

        if (x[mid] < val) && (x[mid+1] > val)
            return mid
        end
 
        # If x is greater, ignore left half
        if x[mid] < val
            low = mid + 1
        # If x is smaller, ignore right half
        elseif x[mid] > val
            high = mid - 1
 
        # means x is present at mid
        else
            return mid
        end
    end
    # If we reach here, then the element was not present
    return low
end

function linear_interpolation(x, y, val, pos, lower_limit, upper_limit)
    # we found pos such that x[pos] < val < x[pos+1]
    # we want to interpolate y(x): such that y[x[pos]] = y[pos], y[x[pos+1]] = y[pos+1]
    # we interpolate y as y(x) = alpha * y2 + (1.0 - alpha) * y1
    # alpha = (val-x[pos]) / (x[pos+1] - x[pos])
    # if pos is out of bounds (-1 if val < x[1] or 0 if val > x[end]) then return lower_limit / upper_limit

    # TODO: store 1.0 / (x[pos+1] - x[pos]) ?

    if pos == -1
        return lower_limit
    elseif pos == 0
        return upper_limit
    end

    alpha = (val - x[pos]) / (x[pos+1] - x[pos])
    return alpha * y[pos+1] + (1.0 - alpha) * y[pos]
end

function find_species_in_db(xml_data, species_name)
    flag = false
    for i in 1:length(xml_data)
        if tag(xml_data[i]) == "Groups"
            for j in 1:length(xml_data[i])
                if attributes(xml_data[i][j])["id"] == species_name
                    flag = true
                    return (i, j)
                end
            end
        end
    end
    if !flag
        return (-1, -1)
    end
end

function load_ionization_data(xml_data)
    ndata = 0
    ΔE = 0.0
    E_data::Vector{Float64} = []
    cs_data::Vector{Float64} = []

    for i in 1:length(xml_data)
        if tag(xml_data[i]) == "Processes"
            for j in 1:length(xml_data[i])
                if attributes(xml_data[i][j])["type"] == "Ionization"  # found reaction
                    for k in 1:length(xml_data[i][j])
                        if tag(xml_data[i][j][k]) == "Parameters"
                            for l in 1:length(xml_data[i][j][k])
                                if tag(xml_data[i][j][k][l]) == "E"
                                    ΔE = parse(Float64, value(xml_data[i][j][k][l][1]))
                                    ndata += 1
                                end
                            end
                        elseif tag(xml_data[i][j][k]) == "DataX"
                            E_data = parse.(Float64, split(value(xml_data[i][j][k][1])))
                            ndata += 1
                        elseif tag(xml_data[i][j][k]) == "DataY"
                            cs_data = parse.(Float64, split(value(xml_data[i][j][k][1])))
                            ndata += 1
                        end
                    end
                end

                if ndata == 3
                    break
                end
            end
        end
    end
    if ndata != 3
        throw(DataMissingException(""))
    else
        return TabulatedCSData(length(E_data), E_data, cs_data, ΔE)
    end
end

function load_elastic_data(xml_data)
    ndata = 0
    ΔE = 0.0
    E_data::Vector{Float64} = []
    cs_data::Vector{Float64} = []

    for i in 1:length(xml_data)
        if tag(xml_data[i]) == "Processes"
            for j in 1:length(xml_data[i])
                if attributes(xml_data[i][j])["type"] == "Elastic"  # found elastic collision
                    for k in 1:length(xml_data[i][j])
                        if tag(xml_data[i][j][k]) == "DataX"
                            E_data = parse.(Float64, split(value(xml_data[i][j][k][1])))
                            ndata += 1
                        elseif tag(xml_data[i][j][k]) == "DataY"
                            cs_data = parse.(Float64, split(value(xml_data[i][j][k][1])))
                            ndata += 1
                        end
                    end
                end

                if ndata == 2
                    break
                end
            end
        end
    end
    if ndata != 2
        throw(DataMissingException(""))
    else
        # TODO: deal with log
        return TabulatedCSData(length(E_data), E_data, cs_data, ΔE)
    end
end


function load_electron_neutral_interactions(species_data, filename, databases, scattering_laws, energy_splits)
    neutral_indexer::Vector{Int32} = []
    elastic_cs_vector::Vector{ElasticScattering} = []
    ionization_cs_vector::Vector{Ionization} = []
    excitation_sink_cs_vector::Vector{ExcitationSink} = []

    xml_data = read(filename, Node)
    n_db = length(xml_data[end]) - 1 # number of databases

    neutral_subindex = 0
    for species in species_data
        if species.charge == 0
            neutral_subindex += 1
            push!(neutral_indexer, neutral_subindex)
            found = false
            for i in 2:2+n_db-1
                if attributes(xml_data[end][i])["id"] == databases[species.name]
                    (ii, jj) = find_species_in_db(xml_data[end][i], species.name)

                    if ii == -1 || jj == -1
                        throw(DataMissingException("No data found for $(species.name) in DB $(databases[species.name])"))
                    end

                    try
                        push!(ionization_cs_vector, Ionization(load_ionization_data(xml_data[end][i][ii][jj]),
                                                                                    scattering_laws[species.name],
                                                                                    energy_splits[species.name]))
                    catch e
                        throw(DataMissingException("No ionization data found for $(species.name) in DB $(databases[species.name])"))
                    end

                    try
                        push!(elastic_cs_vector, ElasticScattering(load_elastic_data(xml_data[end][i][ii][jj]),
                                                                   scattering_laws[species.name]))
                    catch e
                        throw(DataMissingException("No ionization data found for $(species.name) in DB $(databases[species.name])"))
                    end
                    # TODO: Excitation Sink data loading!

                    push!(excitation_sink_cs_vector, ExcitationSink(0, Vector{TabulatedCSData}[], scattering_laws[species.name]))
                    found = true
                end
            end
            if !found
                throw(DataMissingException("No data found for $(species.name) in DB $(databases[species.name])"))
            end
        else
            push!(neutral_indexer, -1)
        end
    end

    return ElectronNeutralInteractions(neutral_subindex, neutral_indexer, elastic_cs_vector, ionization_cs_vector, excitation_sink_cs_vector)
end

function create_computed_crosssections(electron_neutral_interactions)
    res::Vector{ComputedCrossSections} = []

    for i in 1:electron_neutral_interactions.n_neutrals
        push!(res, ComputedCrossSections(electron_neutral_interactions.excitation_sink[i].n_reactions,
                                         0.0, 0.0, 0.0, zeros(electron_neutral_interactions.excitation_sink[i].n_reactions),
                                         zeros(electron_neutral_interactions.excitation_sink[i].n_reactions+2),
                                         zeros(electron_neutral_interactions.excitation_sink[i].n_reactions+2)))
    end

    return res
end

function compute_tabulated_cs_constant_continuation(tabulated_cs_data, E_coll)
    # linear interpolation + constant values (starting and ending array values) for out-of-bounds energies

    e_index = binary_search(tabulated_cs_data.E, E_coll)
    return linear_interpolation(tabulated_cs_data.E,
                                tabulated_cs_data.sigma, E_coll, e_index,
                                tabulated_cs_data.sigma[1],
                                tabulated_cs_data.sigma[end])
end

function compute_tabulated_cs_zero_continuation(tabulated_cs_data, E_coll)
    # linear interpolation + constant values (starting and ending array values) for out-of-bounds energies

    e_index = binary_search(tabulated_cs_data.E, E_coll)
    return linear_interpolation(tabulated_cs_data.E,
                                tabulated_cs_data.sigma, E_coll, e_index,
                                0.0,
                                0.0)
end

function compute_cross_sections!(computed_cs, interaction, g, electron_neutral_interactions, neutral_species_index)
    # the neutral_species_index is the absolute one (i.e. index in the list of all species)
    E_coll_electron_eV = 0.5 * g^2 * e_mass_div_electron_volt  # convert to eV
    E_coll = 0.5 * g^2 * interaction.m_r * eV_J_inv

    i_neutral = electron_neutral_interactions.neutral_indexer[neutral_species_index]

    computed_cs[i_neutral].n_excitations = electron_neutral_interactions.excitation_sink[i_neutral].n_reactions

    computed_cs[i_neutral].cs_elastic = compute_tabulated_cs_constant_continuation(electron_neutral_interactions.elastic[i_neutral].data, E_coll)
    computed_cs[i_neutral].cs_ionization = compute_tabulated_cs_constant_continuation(electron_neutral_interactions.ionization[i_neutral].data, E_coll_electron_eV)

    for i in 1:computed_cs[i_neutral].n_excitations
        computed_cs[i_neutral].cs_excitation[i] = compute_tabulated_cs_constant_continuation(electron_neutral_interactions.excitation_sink[i_neutral].data, E_coll_electron_eV)
    end

    computed_cs[i_neutral].cs_total = computed_cs[i_neutral].cs_elastic + computed_cs[i_neutral].cs_ionization + sum(computed_cs[i_neutral].cs_excitation)

    computed_cs[i_neutral].prob_vec[1] = computed_cs[i_neutral].cs_elastic / computed_cs[i_neutral].cs_total
    computed_cs[i_neutral].prob_vec[2] = computed_cs[i_neutral].cs_ionization / computed_cs[i_neutral].cs_total
    for i in 1:computed_cs[i_neutral].n_excitations
        computed_cs[i_neutral].prob_vec[2+i] = computed_cs[i_neutral].cs_excitation[i] / computed_cs[i_neutral].cs_total
    end

    computed_cs[i_neutral].cdf_prob_vec[1] = 0.0
    computed_cs[i_neutral].cdf_prob_vec[2] = computed_cs[i_neutral].prob_vec[1]

    for i in 1:computed_cs[i_neutral].n_excitations
        computed_cs[i_neutral].cdf_prob_vec[2+i] = computed_cs[i_neutral].prob_vec[1+i] + computed_cs[i_neutral].cdf_prob_vec[1+i]
    end

    return E_coll_electron_eV
end

function compute_cross_sections_only!(computed_cs, interaction, g, electron_neutral_interactions, neutral_species_index)
    # the neutral_species_index is the absolute one (i.e. index in the list of all species)
    E_coll_electron_eV = 0.5 * g^2 * e_mass_div_electron_volt  # convert to eV
    E_coll = 0.5 * g^2 * interaction.m_r * eV_J_inv

    i_neutral = electron_neutral_interactions.neutral_indexer[neutral_species_index]

    computed_cs[i_neutral].n_excitations = electron_neutral_interactions.excitation_sink[i_neutral].n_reactions

    computed_cs[i_neutral].cs_elastic = compute_tabulated_cs_constant_continuation(electron_neutral_interactions.elastic[i_neutral].data, E_coll)
    computed_cs[i_neutral].cs_ionization = compute_tabulated_cs_constant_continuation(electron_neutral_interactions.ionization[i_neutral].data, E_coll_electron_eV)

    for i in 1:computed_cs[i_neutral].n_excitations
        computed_cs[i_neutral].cs_excitation[i] = compute_tabulated_cs_constant_continuation(electron_neutral_interactions.excitation_sink[i_neutral].data, E_coll_electron_eV)
    end
end

function get_cs_total(computed_cs, neutral_species_index, electron_neutral_interactions)
    return computed_cs[electron_neutral_interactions.neutral_indexer[neutral_species_index]].cs_total
end

function get_cs_elastic(computed_cs, neutral_species_index, electron_neutral_interactions)
    return computed_cs[electron_neutral_interactions.neutral_indexer[neutral_species_index]].cs_elastic
end

function get_cs_ionization(computed_cs, neutral_species_index, electron_neutral_interactions)
    return computed_cs[electron_neutral_interactions.neutral_indexer[neutral_species_index]].cs_ionization
end

function get_ionization_threshold(electron_neutral_interactions, neutral_species_index)
    # the neutral_species_index is the absolute one (i.e. index in the list of all species)
    return electron_neutral_interactions.ionization[electron_neutral_interactions.neutral_indexer[neutral_species_index]].data.ΔE
end

function get_electron_energy_split(electron_neutral_interactions, neutral_species_index)
    # the neutral_species_index is the absolute one (i.e. index in the list of all species)
    return electron_neutral_interactions.ionization[electron_neutral_interactions.neutral_indexer[neutral_species_index]].split
end