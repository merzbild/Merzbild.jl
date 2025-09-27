using XML

@muladd begin

"""
    DataMissingException

Exception for the case of missing tabulated cross-section data

# Fields
* `msg`: error message
"""
struct DataMissingException <: Exception
    msg::String
end

"""
    TabulatedCSData

Structure for storing tabulated cross-section data as a function of relative collision energy

# Fields
* `n_vals`: number of values stored
* `E`: array of relative collision energies
* `sigma`: array of cross-section values at these energies
* `ΔE`: how much energy is lost (gained) in the collision in case it is inelastic
"""
struct TabulatedCSData
    n_vals::Int32
    E::Vector{Float64}
    sigma::Vector{Float64}
    ΔE::Float64  # how much energy is lost / threshold energy, eV
end

"""
    ScatteringLaw ScatteringIsotropic=1 ScatteringOkhrimovskyy=2

Enum for various scattering laws in electron-neutral interactions.
`ScatteringIsotropic` corresponds to isotropic scattering,
`ScatteringOkhrimovskyy` to the scattering model of [A. Okhrimovskyy et al., 2002](https://doi.org/10.1103/PhysRevE.65.037402).
"""
@enum ScatteringLaw ScatteringIsotropic=1 ScatteringOkhrimovskyy=2

"""
    ElectronEnergySplit ElectronEnergySplitEqual=1 ElectronEnergySplitZeroE=2

Enum for various splittings of electron energy in electron-impact ionization reactions.
`ElectronEnergySplitEqual` corresponds to energy being shared equally amongst the electrons,
`ElectronEnergySplitZeroE` corresponds to the one-takes-all sharing model.
"""
@enum ElectronEnergySplit ElectronEnergySplitEqual=1 ElectronEnergySplitZeroE=2 # Tanh law

"""
    CSExtend CSExtendConstant=1 CSExtendZero=2
    
Enum defining how to extend cross-sections in case energy is outside of tabulated range.
# Possible values:
* `CSContinueConstant`: continue with closest value in array
* `CSContinueZero`: continue with zero
"""
@enum CSExtend CSExtendConstant=1 CSExtendZero=2

"""
    Ionization

Structure to hold data on an electron-impact ionization cross-section.

# Fields
* `data`: a `TabulatedCSData` instance holding the cross-section values
* `scattering`: the `ScatteringLaw` model for the scattering of the electrons
* `split`: the `ElectronEnergySplit` model for energy splitting across the primary and secondary electrons
"""
struct Ionization
    data::TabulatedCSData
    scattering::ScatteringLaw
    split::ElectronEnergySplit
end

"""
    ElasticScattering

Structure to hold data on an elastic scattering cross-section.

# Fields
* `data`: a `TabulatedCSData` instance holding the cross-section values
* `scattering`: the `ScatteringLaw` model for the scattering of the particles
"""
struct ElasticScattering
    data::TabulatedCSData
    scattering::ScatteringLaw
end

"""
    Ionization

Structure to hold data on electron-impact electronic excitation cross-sections for a specific species.

# Fields
* `n_reactions`: number of electron-impact electronic excitation reactions
* `data`: an array of `TabulatedCSData` instances (of length `n_reactions`) holding the cross-section values for each reaction
* `scattering`: the `ScatteringLaw` model for the scattering of the particles
"""
struct ExcitationSink
    n_reactions::Int32
    data::Vector{TabulatedCSData}
    scattering::ScatteringLaw
end

"""
    ElectronNeutralInteractions

Structure to hold data on electron-neutral interactions. The `neutral_indexer` field is used to obtain
the index of the species inside the structure, given an index of the neutral species in the full list of species
in the simulation. For example, if we have a following list of species in the simulation: `[e-, He, Ar+, He+, Ar]`,
`n_neutrals=2`, the `ElectronNeutralInteractions` instance stores data for interactions of electrons with `[He, Ar]`,
and `neutral_indexer[2] = 1`, `neutral_indexer[5] = 2`.

# Fields
* `n_neutrals`: number of neutral species for which the data has been loaded
* `neutral_indexer`: array that maps indices of neutral species
    in the full list of species to local indices in the `ElectronNeutralInteractions` structure
* `elastic`: an array of `ElasticScattering` instances (of length `n_neutrals`) holding the cross-section data on elastic scattering
    for each neutral species
* `ionization`: an array of `Ionization` instances (of length `n_neutrals`) holding the cross-section data on electron-impact ionization
    for each neutral species
* `excitation_sink`: an array of `ExcitationSink` instances (of length `n_neutrals`)
    holding the cross-section data on electron-impact electronic excitation for each neutral species
"""
struct ElectronNeutralInteractions
   n_neutrals::Int32
   neutral_indexer::Vector{Int32}  # neutral_indexer[species_index] = index of the species in the sub-list of neutral species
   elastic::Vector{ElasticScattering}
   ionization::Vector{Ionization}
   excitation_sink::Vector{ExcitationSink}
end

"""
    ComputedCrossSections

Structure to hold data on computed cross-sections of electron-neutral interactions for a specific neutral species.

# Fields
* `n_excitations`: number of electron-impact excitation reactions
* `cs_total`: the computed total cross-section (sum of cross-sections of all processes)
* `cs_elastic`: the computed elastic scattering cross-section
* `cs_ionization`: the computed electron-impact ionization cross-section
* `cs_excitation`: the computed electron-impact electronic excitation cross-section
* `prob_vec`: a vector of probabilities of the processes (of length `2+n_excitations`). `prob_vec[1]` is the probability of elastic scattering,
    `prob_vec[2]` is the probability of electron-impact ionization, `prob_vec[3:2+n_excitations]` are the probabilities
    of th 
* `cdf_prob_vec`: a vector of cumulative probabilities of the processes (of length `3+n_excitations`), used for sampling a specific process:
    `cfd_prob_vec[1] = 0.0`, `cfd_prob_vec[n] = cfd_prob_vec[n-1] + prob_vec[n-1], n>1`
"""
mutable struct ComputedCrossSections
    n_excitations::Int32
    cs_total::Float64
    cs_elastic::Float64
    cs_ionization::Float64
    cs_excitation::Vector{Float64}
    prob_vec::Vector{Float64}
    cdf_prob_vec::Vector{Float64}
end

"""
    sigma_vhs(interaction, g)

Computes the VHS cross-section.

# Positional arguments
* `interaction`: the `Interaction` instance
* `g`: the relative velocity of the collision

# Returns
The value of the computed cross-section.
"""
@inline function sigma_vhs(interaction, g)
    return interaction.vhs_factor * g^(1.0 - 2 * interaction.vhs_o)
end

"""
    binary_search(x, val)

Binary search for value `val` in a sorted array `x`. Finds the position `mid` such that
`x[mid] < val < x[mid+1]`.

# Positional arguments
* `x`: the sorted array to be searched
* `val`: the value to search for

# Returns
* `-1` if `val < x[1]`
* `0` if `val > x[end]`
* Otherwise, returns the index `mid` satisfying `x[mid] < val < x[mid+1]`
"""
function binary_search(x, val)
    low = 1
    high = length(x)
    mid = 0

    @inbounds if val > x[high]
        return 0
    elseif val < x[low]
        return -1
    end
 
    @inbounds while low <= high
 
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

"""
    linear_interpolation(x, y, val, pos, lower_limit, upper_limit)

Perform linear interpolation of a tabulated function `y(x)`, given a parameter value `val`,
a sorted array of values of `x` and corresponding values of `y`, as well as the index of the closest values
of `x` to `val`.
In case `val` is outside of the limits of `x`, return placeholder values.
That is, given that `x[pos] < val < x[pos+1]`, we want to interpolate `y(val)` given that
`y[x[pos]] = y[pos], y[x[pos+1]] = y[pos+1]`.

# Positional arguments
* `x`: the vector of values of the parameter of the function to be interpolated
* `y`: the vector of the corresponding function values
* `val`: the value of the parameter at which the function is to be interpolated
* `pos`: the index of the element in `x` such that `x[pos] < val < x[pos+1]`, 
    or `-1` if `val < x[1]`, `0` if `val > x[end]`
* `lower_limit`: the value to return if `val < x[1]`
* `upper_limit`: the value to return if `val > x[end]`

# Returns
* The linearly interpolated value of `y(val)` if `x[1] <= val <= x[end]`
* `lower_limit`, if `val < x[1]`
* `upper_limit`, if `val > x[end]`
"""
@inline function linear_interpolation(x, y, val, pos, lower_limit, upper_limit)
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

    @inbounds alpha = (val - x[pos]) / (x[pos+1] - x[pos])
    @inbounds return alpha * y[pos+1] + (1.0 - alpha) * y[pos]
end


"""
Find a chemical species in an LXCat-format XML; returns the first instance of the species
found and the id of the "Groups" element in which the species was found.
The following conditions need to be satisfied for a tuple `(i,j)` for it to be returned:
`tag(xml_data[i]) == "Groups"` and `attributes(xml_data[i][j])["id"] == species_name`.

# Positional arguments
* `xml_data`: the LXCat-format XML data to be searched
* `species_name`: the name of the species to search for

# Returns
Tuple `(i,j)` denoting the index of the group and group element in which the species is located
"""
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

"""
    load_ionization_data(xml_data)

Load electron-impact ionization data from the part of an LXCat-format XML file for
a specific species.
The following should hold for the data to be loaded:
`tag(xml_data[i]) == "Processes"`, `attributes(xml_data[i][j])["type"] == "Ionization"`.
This will load the first set of data found for the ionization process.

# Positional arguments
* `xml_data`: the LXCat-format XML data to be searched

# Returns
`TabulatedCSData` structure containing the electron-impact ionization cross-section data.

# Throws
`DataMissingException` if data not found or not all required data present.
"""
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

"""
    load_elastic_data(xml_data)

Load electron-neutral elastic scattering data from the part of an LXCat-format XML file for
a specific species.
The following should hold for the data to be loaded:
`tag(xml_data[i]) == "Processes"`, `attributes(xml_data[i][j])["type"] == "Ionization"`.
This will load the first set of data found for the ionization process.

# Positional arguments
* `xml_data`: the LXCat-format XML data to be searched

# Returns
`TabulatedCSData` structure containing the electron-neutral elastic scattering cross-section data.

# Throws
`DataMissingException` if data not found or not all required data present.
"""
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


"""
    load_electron_neutral_interactions(species_data, filename, databases, scattering_laws, energy_splits)

Load electron-neutral interaction data from an LXCAT format XML file for
a set of given neutral species.

# Positional arguments
* `species_data`: vector of `Species` data for the neutral species
* `filename`: path to XML file
* `databases`: dictionary of `(species.name => database)` pairs, specifying the name
of the cross-section database in the XML file to use for the species 
* `scattering_laws`: vector of `ScatteringLaw` instances to use for each species
* `energy_splits`: vector of `ElectronEnergySplit` instances to use for each species

# Returns
`ElectronNeutralInteractions` structure containing the electron-neutral interaction data.

# Throws
`DataMissingException` if data not found or not all required data present.
"""
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
                        throw(DataMissingException("No scattering data found for $(species.name) in DB $(databases[species.name])"))
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

"""
    create_computed_crosssections(electron_neutral_interactions)

Create a vector of `ComputedCrossSection` instances for the electron-neutral interactions.

Positional arguments
* `electron_neutral_interactions`: the `ElectronNeutralInteractions` instance for which the
    cross-sections will be computed

# Returns
Vector of `ComputedCrossSection` of length `electron_neutral_interactions.n_neutrals`.
"""
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

"""
    compute_tabulated_cs_constant_continuation(tabulated_cs_data, E_coll)

Compute an energy-dependent cross-section from tabulated data using linear interpolation.
If the energy is smaller than the minimum energy used for the tabulation, return first element of the table;
if energy is larger than the maximum energy used for the tabulation, return last element of the table.

# Positional arguments
* `tabulated_cs_data`: a `TabulatedCSData` instance with the tabulated energies and cross-section values
* `E_coll`: the collision energy for which to compute the cross-section

# Returns
The value of the cross-section.
"""
function compute_tabulated_cs_constant_continuation(tabulated_cs_data, E_coll)
    # linear interpolation + constant values (starting and ending array values) for out-of-bounds energies

    e_index = binary_search(tabulated_cs_data.E, E_coll)
    return linear_interpolation(tabulated_cs_data.E,
                                tabulated_cs_data.sigma, E_coll, e_index,
                                tabulated_cs_data.sigma[1],
                                tabulated_cs_data.sigma[end])
end

"""
    compute_tabulated_cs_zero_continuation(tabulated_cs_data, E_coll)

Compute an energy-dependent cross-section from tabulated data using linear interpolation.
If the energy is outside of the range of energies used for the tabulation, returns 0.0.

# Positional arguments
* `tabulated_cs_data`: a `TabulatedCSData` instance with the tabulated energies and cross-section values
* `E_coll`: the collision energy for which to compute the cross-section

# Returns
The value of the cross-section.
"""
function compute_tabulated_cs_zero_continuation(tabulated_cs_data, E_coll)
    # linear interpolation + constant values (starting and ending array values) for out-of-bounds energies

    e_index = binary_search(tabulated_cs_data.E, E_coll)
    return linear_interpolation(tabulated_cs_data.E,
                                tabulated_cs_data.sigma, E_coll, e_index,
                                0.0,
                                0.0)
end

"""
    compute_cross_sections_only!(computed_cs, interaction, g, electron_neutral_interactions, neutral_species_index, extend::CSExtend)

Compute electron-impact ionization and excitation cross-sections, and electron-neutral elastic scattering
crosss-sections, and return collision energy in eV. Here `neutral_species_index` is the index of the neutral
species being considered in the overall array of the `Species` instances that includes all species in the simulation
and was used to construct the `ElectronNeutralInteractions` instance.

# Positional arguments
* `computed_cs`: the vector of `ComputedCrossSection` instances in which the computed values will be stored
* `interaction`: the `Interaction` instance describing the electron-neutral interaction being considered
* `g`: the magnitude of the relative collision velocity
* `electron_neutral_interactions`: the `ElectronNeutralInteractions` instance storing the tabulated cross-section
    data
* `neutral_species_index`: the index of the neutral species being considered

# Returns
The electron-neutral collision energy in eV.
"""
function compute_cross_sections_only!(computed_cs, interaction, g, electron_neutral_interactions, neutral_species_index, extend::CSExtend)
    # the neutral_species_index is the absolute one (i.e. index in the list of all species)
    # this is used in the NNLS merging approach
    E_coll_electron_eV = 0.5 * g^2 * e_mass_div_electron_volt  # convert to eV
    E_coll = 0.5 * g^2 * interaction.m_r * eV_J_inv

    i_neutral = electron_neutral_interactions.neutral_indexer[neutral_species_index]

    computed_cs[i_neutral].n_excitations = electron_neutral_interactions.excitation_sink[i_neutral].n_reactions

    if extend == CSExtendConstant
        computed_cs[i_neutral].cs_elastic = compute_tabulated_cs_constant_continuation(electron_neutral_interactions.elastic[i_neutral].data, E_coll)
        computed_cs[i_neutral].cs_ionization = compute_tabulated_cs_constant_continuation(electron_neutral_interactions.ionization[i_neutral].data, E_coll_electron_eV)

        for i in 1:computed_cs[i_neutral].n_excitations
            computed_cs[i_neutral].cs_excitation[i] = compute_tabulated_cs_constant_continuation(electron_neutral_interactions.excitation_sink[i_neutral].data, E_coll_electron_eV)
        end
    else
        computed_cs[i_neutral].cs_elastic = compute_tabulated_cs_zero_continuation(electron_neutral_interactions.elastic[i_neutral].data, E_coll)
        computed_cs[i_neutral].cs_ionization = compute_tabulated_cs_zero_continuation(electron_neutral_interactions.ionization[i_neutral].data, E_coll_electron_eV)

        for i in 1:computed_cs[i_neutral].n_excitations
            computed_cs[i_neutral].cs_excitation[i] = compute_tabulated_cs_zero_continuation(electron_neutral_interactions.excitation_sink[i_neutral].data, E_coll_electron_eV)
        end
    end

    return E_coll_electron_eV
end


"""
    compute_cross_sections!(computed_cs, interaction, g, electron_neutral_interactions, neutral_species_index; extend::CSExtend=CSExtendConstant)

Compute electron-impact ionization and excitation cross-sections, electron-neutral elastic scattering
cross-sections, total collision cross-section, probabilities of the different processes, and return collision energy in eV.
Here `neutral_species_index` is the index of the neutral
species being considered in the overall array of the `Species` instances that includes all species in the simulation
and was used to construct the `ElectronNeutralInteractions` instance.
Out-of-tabulation-range values are treated according to what `extend` method is used.

# Positional arguments
* `computed_cs`: the vector of `ComputedCrossSection` instances in which the computed values will be stored
* `interaction`: the `Interaction` instance describing the electron-neutral interaction being considered
* `g`: the magnitude of the relative collision velocity
* `electron_neutral_interactions`: the `ElectronNeutralInteractions` instance storing the tabulated cross-section
    data
* `neutral_species_index`: the index of the neutral species being considered

# Keyword arguments
* `extend`: enum of `CSExtend` type that sets how out-of-range energy values are treated when computing cross-sections

# Returns
The electron-neutral collision energy in eV.
"""
function compute_cross_sections!(computed_cs, interaction, g, electron_neutral_interactions, neutral_species_index; extend::CSExtend=CSExtendConstant)
    # the neutral_species_index is the absolute one (i.e. index in the list of all species)
    E_coll_electron_eV = compute_cross_sections_only!(computed_cs, interaction, g,
                                                      electron_neutral_interactions, neutral_species_index, extend)

    i_neutral = electron_neutral_interactions.neutral_indexer[neutral_species_index]

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

"""
    get_cs_total(electron_neutral_interactions, computed_cs, neutral_species_index)

Get the value of a computed total collision cross-section for electron-neutral interactions.
Here `neutral_species_index` is the index of the neutral
species being considered in the overall array of the `Species` instances that includes all species in the simulation
and was used to construct the `ElectronNeutralInteractions` instance.

# Positional arguments
* `electron_neutral_interactions`: the `ElectronNeutralInteractions` instance storing the tabulated cross-section
    data
* `computed_cs`: the vector of `ComputedCrossSection` instances holding the computed cross-section values
* `neutral_species_index`: the index of the neutral species being considered

# Returns
The value of the total collision cross-section.
"""
@inbounds function get_cs_total(electron_neutral_interactions, computed_cs, neutral_species_index)
    return computed_cs[electron_neutral_interactions.neutral_indexer[neutral_species_index]].cs_total
end

"""
    get_cs_elastic(electron_neutral_interactions, computed_cs, neutral_species_index)

Get the value of a computed electron-neutral elastic scattering cross-section.
Here `neutral_species_index` is the index of the neutral
species being considered in the overall array of the `Species` instances that includes all species in the simulation
and was used to construct the `ElectronNeutralInteractions` instance.

# Positional arguments
* `electron_neutral_interactions`: the `ElectronNeutralInteractions` instance storing the tabulated cross-section
    data
* `computed_cs`: the vector of `ComputedCrossSection` instances holding the computed cross-section values
* `neutral_species_index`: the index of the neutral species being considered

# Returns
The value of the elastic scattering cross-section.
"""
@inbounds function get_cs_elastic(electron_neutral_interactions, computed_cs, neutral_species_index)
    return computed_cs[electron_neutral_interactions.neutral_indexer[neutral_species_index]].cs_elastic
end

"""
    get_cs_ionization(electron_neutral_interactions, computed_cs, neutral_species_index)

Get the value of a computed electron-impact ionization cross-section.
Here `neutral_species_index` is the index of the neutral
species being considered in the overall array of the `Species` instances that includes all species in the simulation
and was used to construct the `ElectronNeutralInteractions` instance.

# Positional arguments
* `electron_neutral_interactions`: the `ElectronNeutralInteractions` instance storing the tabulated cross-section
    data
* `computed_cs`: the vector of `ComputedCrossSection` instances holding the computed cross-section values
* `neutral_species_index`: the index of the neutral species being considered

# Returns
The value of the electron-impact ionization scattering cross-section.
"""
@inbounds function get_cs_ionization(electron_neutral_interactions, computed_cs, neutral_species_index)
    return computed_cs[electron_neutral_interactions.neutral_indexer[neutral_species_index]].cs_ionization
end

"""
    get_ionization_threshold(electron_neutral_interactions, neutral_species_index)

Get the ionization threshold energy.
Here `neutral_species_index` is the index of the neutral
species being considered in the overall array of the `Species` instances that includes all species in the simulation
and was used to construct the `ElectronNeutralInteractions` instance.

# Positional arguments
* `electron_neutral_interactions` the `ElectronNeutralInteractions` instance storing the tabulated cross-section
    data
* `neutral_species_index`: the index of the neutral species being considered

# Returns
The ionization threshold energy of a specific neutral species.
"""
@inbounds function get_ionization_threshold(electron_neutral_interactions, neutral_species_index)
    # the neutral_species_index is the absolute one (i.e. index in the list of all species)
    return electron_neutral_interactions.ionization[electron_neutral_interactions.neutral_indexer[neutral_species_index]].data.ΔE
end

"""
    get_electron_energy_split(electron_neutral_interactions, neutral_species_index)

Get the way electron energy is split during ionization for a specific electron-neutral interaction.
Here `neutral_species_index` is the index of the neutral
species being considered in the overall array of the `Species` instances that includes all species in the simulation
and was used to construct the `ElectronNeutralInteractions` instance.

# Positional arguments
* `electron_neutral_interactions` the `ElectronNeutralInteractions` instance storing the tabulated cross-section
    data
* `neutral_species_index`: the index of the neutral species being considered

# Returns
The `ElectronEnergySplit` value for the specific electron-neutral interaction.
"""
@inbounds function get_electron_energy_split(electron_neutral_interactions, neutral_species_index)
    # the neutral_species_index is the absolute one (i.e. index in the list of all species)
    return electron_neutral_interactions.ionization[electron_neutral_interactions.neutral_indexer[neutral_species_index]].split
end

end