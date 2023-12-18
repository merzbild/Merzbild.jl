# NNLSQ

## Grid in Cartesian coordinates

Two ways to generate noise/remove bias: coordinate-wise, or for each grid point (more expensive)
Construct velocity moment matrix and solve

## Grid in spherical coordinates

Equal volume/equal radius?

# TODO

- [x] NTC collisions

- [x] I/O of particle properties (basic)

- [x] Particle generation (equal weights)

- [x] Properties computation

    - [x] basic ones

    - [x] moments

- [x] I/O of output

- [x] Maxwellian test case (equal weights)

- [ ] BKW test case (equal weights)

- [ ] Particle generation (variable weight)

- [ ] Grid merging

- [ ] Octree merging

- [ ] mixing rule VHS creator

- [ ] The science begins

# TODO: tests

- [ ] energy / momentum conservation in scattering

- [ ] correct indexing

- [ ] 1D - no merging, particles don't switch cells during variable weight collisions!!!