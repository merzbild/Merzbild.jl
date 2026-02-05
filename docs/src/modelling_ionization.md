# Modelling ionization reactions

## Overview
Currently, it is possible to model electron-impact ionization of neutrals, as well as electron-neutral scattering.
Only isotropic (VHS) scattering is currently implemented.

In order to do so, a structure for storing tabulated cross-section data, as well as a structure for storing computed
cross-section values and corresponding probabilities, are required. Cross-section data is assumed to be given as a function
of the collision energy in electron-volt, and a binary search is used to find the closest energy values to the given collision energy,
with subsequent linear interpolation of the cross-sections.

## Loading cross-section data
Currently, it is possible to load XML data in the [LXCat](https://nl.lxcat.net/home/) format. As the XML format for LXCat data
is not finalized, some data exported from LXCat need to be edited manually, as they might miss species identifiers, etc.
Examples of synthetic data used for testing can be found in `data/test_neutral_electron_data.xml`.

In general, the XML file has to have the following structure
```xml
<lxcat version="" created="" message="example structure">
    <Database name="Linear dependence of cross-section on energy, He" id="LinearDB">
        <Groups>
            <Group id="He">
                <Processes>
                    <Process class="Scattering Cross Sections" type="Elastic">
                        <Species>
                            <Reactant>e</Reactant><Reactant>He</Reactant>
                            <Product>E</Product><Product>He</Product>
                        </Species>
                        
                        <Reaction>
                        E + He -&gt; E + He
                        </Reaction>
                        
                        <Parameters>
                            <mM>
                            1.360000e-4
                            </mM>
                            <Parameter>
                            complete set
                            </Parameter>
                        </Parameters>
                        <DataX type="Energy" units="eV" size="4">
                        0.000000e+0 1.0e2 5.0e2 1.0e3
                        </DataX>
                        <DataY type="Cross section" units="m2" size="4">
                        4.0e-20 5.0e-20 9.0e-20 1.4e-19
                        </DataY>
                    </Process>

                    <Process class="Scattering Cross Sections" type="Ionization">
                    <Species>
                        <Reactant>e</Reactant><Reactant>He</Reactant>
                        <Product>E</Product><Product>E</Product><Product>He^+</Product>
                    </Species>
                    
                    <Reaction>
                    E + He -&gt; E + E + He^+
                    </Reaction>
                    
                    <Parameters>
                        <E units="eV">
                        2.458740e+1
                        </E>
                        <Parameter>
                        complete set
                        </Parameter>
                    </Parameters>
                    <DataX type="Energy" units="eV" size="3">
                    2.458739e+1 3.0e+1 1.0e+2
                    </DataX>
                    
                    <DataY type="Cross section" units="m2" size="3">
                    0.000000e+0 1.0e-25 1.0e-24 
                    </DataY>
                    </Process>
                </Processes>
            </Group>
        </Groups>
    </Database>
</lxcat>
```

The `Database` refers to a specific source for the cross-section data. Within a `Database`, the collision partner species are
grouped by `Group`, within a specific `Group` a specific `Process` has an associated type, and potentially a threshold energy
(i.e. for ionization or electronic excitation reactions) given in `Parameters`. The energies used for the tabulation are in
a `DataX` field, whereas the cross-section values are in a `DataY` field.

The tabulated values for the different processes, along with the threshold energy values are stored in
a [`ElectronNeutralInteractions`](@ref) structure, which holds all the electron-neutral interactions for a system.
To instantiate such an instance, the [`load_electron_neutral_interactions`](@ref) is called. The type of scattering
used for each species is set by the `scattering_laws` parameter, and how the energy split across primary and secondary electrons
in ionization reactions is governed by the `energy_splits` parameter.

A second type, [`ComputedCrossSections`](@ref) stores information about the computed cross-sections and probabilities for a specific collision
pair; each one is species-specific, therefore a vector for all the neutral species can be instantiated by calling
[`create_computed_crosssections`](@ref) and passing an [`ElectronNeutralInteractions`](@ref) instance.

## Performing ionizing collisions
Electron-neutral collisions with elastic scattering and ionization reactions are modelled by the [`ntc_n_e!`] and
[`ntc_n_e_es!`] functions. The latter implements the Event Splitting
of [Oblapenko et al. (2022)](https://doi.org/10.1016/j.jcp.2022.111390) and reduces the level of stochastic noise,
but is suitable only for variable-weight simulations with particle merging.
The function implements the collision mechanics as described in [Nanbu (2000)](https://doi.org/10.1109/27.887765);
apart from taking a `ElectronNeutralInteractions` and a `Vector` of `ComputedCrossSections` (with
`n_neutral_species` elements). Additionally, an `extend` keyword argument of enum type [`Merzbild.CSExtend`](@ref)
can be set that governs how cross-sections will be computed for energies that are outside of the range of tabulated values.
In this case the cross-sections are either assumed to be equal to the nearest in-range value (`extend=CSExtendConstant`, default
behaviour), or are set to 0 (`extend=CSExtendZero`).

**Important**: due to the large velocity discrepancies between neutrals and electrons, as well as potentially non-trivial cross-section
behaviour, use of a standard estimate for `collision_factors` (that is based on thermal velocities and VHS cross-section values)
can lead to a significant over- or under-estimation of the number of collisions to perform.
To avoid this, it is suggested to use [`estimate_sigma_g_w_max_ntc_n_e!`](@ref), a function that 
samples multiple collision pairs (without performing the collisions)
and updates the collision factors based on the actual computed cross-section values.
It samples ``N_{coll} + min_{coll}`` pairs, where ``N_{coll}`` is the number of pairs computed using the standard
formulas for the NTC routine, and ``min_{coll}}`` is the minimum number of pairs to sample (default values of 15).
The whole process (compute ``N_{coll}`` and sample, updating collision factors) is repeated `n_loops` times
(default value of 6).

An example 0-D simulation with a constant applied electric field can be found under
`simulations/0D/ionization/0D_ionization_1neutralspecies_es.jl`. By default it uses IST Lisbon data from LXCat which should
be downloaded separately.