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

- [x] Two species relaxation (equal weights)

- [x] BKW test case (equal weights)

- [x] Particle generation (variable weight)

- [x] Grid merging

- [ ] Octree merging

- [ ] mixing rule VHS creator

- [ ] The science begins


# TODO: features

- [ ] Add time to output (since we can change dt on the fly)

- [ ] Compute sigma_g_vhs directly (to avoid additional multiplication)

- [ ] Avoid duplicate computation of particle indices in grid-based merging? Avoid second loop? Use Welford's algorithm + initial estimate of mean as middle of box? Or just shift data?

# TODO: tests

- [ ] energy / momentum conservation in scattering

- [x] correct indexing

- [ ] BKW var weight reference solution

- [ ] 1D - no merging, particles don't switch cells during variable weight collisions!!!
