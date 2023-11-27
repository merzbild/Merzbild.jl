function sigma_vhs(interaction, g)
    return interaction.vhs_factor * g^(2.0 - 2 * interaction.vhs_o)
end