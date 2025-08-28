# Overview of capabilities

This page provides an overview of the main capabilities of the code.

## Grids, fixed and variable weight DSMC
|                        | **0D**                                        | **1D** |
|:----------------------:|:-----------------------------------------:|:----:|
| Fixed-weight DSMC      | ✅                                        | ✅ |
| Variable-weight DSMC   | ✅ |  ✅ |

Currently, simulations are either 0D (spatially homogeneous) or on a uniform 1D grid.
For 1D simulations, specular and diffuse reflection models are available (along with a mixture of the two
via an accommodation coefficient).

## Collisions
The No-Time-Counter (NTC) approach of [Bird (1994)](https://doi.org/10.1093/oso/9780198561958.001.0001) is implemented for fixed-weight DSMC.
The variable-weight NTC approach of [Schmidt and Rutland (2000)](https://doi.org/10.1006/jcph.2000.6568) is implemented for variable-weight DSMC.
Event splitting ([Oblapenko et al. (2022)](https://doi.org/10.1016/j.jcp.2022.111390)) is implemented for neutral-electron interactions.

## Fokker-Planck collisions
As an alternative to DSMC, one can use the stochastic Fokker-Planck algorithm to simulate the particle collisions.
Currently, the linear Fokker-Planck model of [Gorhi, Torrilhon, and Jenny (2011)](https://doi.org/10.1017/jfm.2011.188) is implemented for fixed-weight particles.

|                        | **Linear**                                        | **Cubic** |
|:----------------------:|:-----------------------------------------:|:----:|
| Fixed-weight Fokker-Planck      | ✅                                        | ❌ |
| Variable-weight Fokker-Planck   | ❌ | ❌ |

## Particle merging algorithms
The following particle merging algorithms are available for variable-weight DSMC
simulations:

  1. A grid-based merging algorithm as described in [Oblapenko et al. (2020)](https://doi.org/10.1016/j.jcp.2020.109302) (see also [Vranic et al. (2015)](https://doi.org/10.1016/j.cpc.2015.01.020))
  2. The octree merging algorithm of [Martin and Cambier (2016)](https://doi.org/10.1016/j.jcp.2016.01.020)
  3. A Non-Negative Least Squares (NNLS)-based merging approach described in [Oblapenko (2024)](https://doi.org/10.48550/arXiv.2412.12354)

## Cross-sections
The Variable-Hard Sphere (VHS) model is implemented for collisions of neutral particles.
For neutral-electron collisions, [LXCat](http://www.lxcat.net) data in XML format needs to be provided
for the elastic scattering and electron-impact ionization cross-sections. Currently, only
isotropic scattering is implemented for neutral-electron collisions.

## Inelastic collisions
In flows with neutrals, ions, and electrons, electron-impact ionization is supported.
Variable weight DSMC simulations also support the Event Splitting (ES) collision algorithm
of [Oblapenko et al. (2022)](https://doi.org/10.1016/j.jcp.2022.111390).

## External fields
Acceleration of charged particles by a constant electric field is supported.

## I/O
The code assumes the TOML format for the particle and VHS interaction data. XML is used
for the LXCat data.
Output of the computed macroscopic properties is in NetCDF4 format.

## Parallel computations
Multithreading is supported via domain decomposition,
see [the section on setting up multithreaded simulations](@ref "Multithreaded simulations") for more details.
Multithreaded variable-weight DSMC computations have not been tested yet.