struct TabulatedCSData
    E::Vector{Float64}
    sigma::Vector{Float64}
    ΔE::Float64  # how much energy is lost / threshold energy
    species_loss::Int8  # which species loses energy
end

@enum ScatteringLaw Isotropic=1 Okhrimovskyy=2
@enum ElectronEnergySplit Equal=1 ZeroE=2 # Tanh law
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
    data::Vector{TabulatedCSData}
    scattering::ScatteringLaw
end

function sigma_vhs(interaction, g)
    return interaction.vhs_factor * g^(1.0 - 2 * interaction.vhs_o)
end

