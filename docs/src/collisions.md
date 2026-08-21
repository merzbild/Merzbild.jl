# Collision algorithms

Merzbild.jl implements several algorithms for modelling binary collisions:
either directly or via the Fokker-Planck approach. All of them operate on the particles of a single grid cell,
using the [`ParticleIndexerArray`](@ref) to access the particles of the relevant species.
The variable-weight NTC variants split particles when the collision partners have unequal
weights; the equal-weight variant skips these checks for performance.

## Overview

| Algorithm | Weights | Species | Function |
| --- | --- | --- | --- |
| NTC (No-Time-Counter) DSMC | equal-weight | single species | [`ntc_equal_weight!`](@ref) |
| NTC DSMC | equal-weight | two species | [`ntc_equal_weight!`](@ref) |
| NTC DSMC | variable-weight | single species | [`ntc!`](@ref) |
| NTC DSMC | variable-weight | two species | [`ntc!`](@ref) |
| SWPM (Stochastic Weighted Particle Method) | variable-weight | single species | [`swpm!`](@ref) |
| Linear Fokker–Planck | variable-weight | single species | [`fp_linear!`](@ref) |
| NTC electron–neutral (elastic scattering + electron-impact ionization) | variable-weight | electron/neutral/ion | [`ntc_n_e!`](@ref) |
| NTC electron–neutral with event splitting | variable-weight | electron/neutral/ion | [`ntc_n_e_es!`](@ref) |

## Elastic scattering models

The DSMC and SWPM routines listed above support several models of elastic scattering, which define
both the total collision cross-section ``\sigma(g)`` (and thus the collision probability) and the
angular scattering law (and thus the post-collisional velocities):

| Model | Type tag | Enum value | Total cross-section | Scattering |
| --- | --- | --- | --- | --- |
| Variable Hard Sphere | [`VHS`](@ref) | `ScatteringVHS` | ``C g^{1 - 2\omega}`` | isotropic |
| Variable Soft Sphere | [`VSS`](@ref) | `ScatteringVSS` | ``C g^{1 - 2\omega}`` | ``\cos\chi = 2 R^{1/\alpha} - 1`` |

The model is fixed for a species pair and is stored in the pair's [`Interaction`](@ref) instance.
It is set by the optional `model` key of the pair's entry in the interaction data TOML file,
defaulting to `"VHS"` if the key is absent:

```toml
["Ar,Ar"]
model = "VSS"
vhs_d = 4.11e-10
vhs_o = 0.81
vhs_Tref = 273.0
vss_alpha = 1.40
```

The VSS model requires the additional `vss_alpha` key; a bundled interaction data file using it is
`vss.toml` (see [`MERZBILD_DATA_PATH`](@ref)). Since the VSS model modifies the angular scattering law
but not the total cross-section, the reference viscosity of a VSS interaction is corrected by the
factor ``(\alpha + 1)(\alpha + 2)/(6\alpha)``; ``\alpha = 1`` recovers isotropic scattering, and
the VSS model then coincides with the VHS model.

A hard sphere gas is a special case of the VHS model, obtained by setting ``\omega = 1/2`` in the
interaction data file: the cross-section is then the constant ``\pi d^2``.

The collision routines resolve the model of the species pair once per call, and are then compiled
for that single model, so that no branching on the model, no dynamic dispatch and no allocations
are left in the collision loop itself. The model can also be passed explicitly as the second
positional argument (e.g. `ntc!(rng, VSS(), collision_factors, ...)`), which overrides the model
stored in the `Interaction` instance.

## DSMC collisions

The No-Time-Counter (NTC) algorithm is the default DSMC collision routine. Four variants are
available:

* [`ntc_equal_weight!`](@ref) — single-species collisions assuming all particles share the same
  computational weight; no weight checks or particle splitting are performed. This is the classic
  NTC scheme of [Bird (1994)](https://doi.org/10.1093/oso/9780198561958.001.0001).
* [`ntc!`](@ref) — single-species collisions for variable-weight particles. When the collision
  partners' weights differ by more than `dw_tol`, particle splitting is used to conserve mass,
  momentum, and energy, following the variable-weight NTC approach of
  [Schmidt and Rutland (2000)](https://doi.org/10.1006/jcph.2000.6568).
* [`ntc!`](@ref) (two-species method) — collisions between particles of two different species,
  again supporting variable weights and particle splitting
  ([Schmidt and Rutland (2000)](https://doi.org/10.1006/jcph.2000.6568)).
* [`ntc_equal_weight!`](@ref) (two-species method) — collisions between particles of two different
  species with equal computational weights.

## SWPM

[`swpm!`](@ref) implements the Stochastic Weighted Particle Method of
[Rjasanow and Wagner (2005)](https://doi.org/10.1007/3-540-27689-0) for single-species,
variable-weight collisions. During a collision of particles with weights ``w_i``, ``w_j``, the
weights are depleted by ``\min(w_i, w_j) / (1+G)``, where ``G \geq 0`` is a user-defined
parameter controlling the weight-transfer function.

## Fokker–Planck

[`fp_linear!`](@ref) models single-species elastic collisions using a linear Fokker–Planck
approximation instead of resolving individual binary collisions, following the model of
[Gorji, Torrilhon, and Jenny (2011)](https://doi.org/10.1017/jfm.2011.188).

## Electron–neutral collisions

For plasma simulations, dedicated NTC routines model electron–neutral elastic scattering together
with electron-impact ionization, using tabulated cross-sections stored in an
`ElectronNeutralInteractions` instance:

* [`ntc_n_e!`](@ref) — electron–neutral elastic scattering and electron-impact ionization without
  event splitting.
* [`ntc_n_e_es!`](@ref) — the same processes modelled with the event-splitting method of
  [Oblapenko et al. (2022)](https://doi.org/10.1016/j.jcp.2022.111390), which reduces statistical
  noise by splitting collision outcomes rather than sampling a single outcome per pair.

Because the electron–neutral cross-sections vary strongly with collision energy, the value of
``(\sigma g w)_{max}`` needed by the NTC algorithm can be estimated stochastically with
[`estimate_sigma_g_w_max_ntc_n_e!`](@ref) before performing the collisions.
