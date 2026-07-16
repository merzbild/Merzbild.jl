# Merzbild.jl benchmarks

Various benchmarks and comparisons to other open-source codes are provided here for reference test cases. An overview of benchmarks between different versions of Merzbild is provided in the Summary section.

## Summary

### Couette flow, serial, small grid

|                         | **M1 Pro, 32 GB RAM** |  **Intel Core i9-13900K, 128 GB RAM** |  **AMD EPYC 9374F, 378 GB RAM** |
|:-----------------------:|:---------------------:|:-------------------------------------:|:-------------------------------:|
| v0.7.0 (Julia 1.11)     |       33.5s           |    28.4s                              |             33.0s               |
| v0.8.0 (Julia 1.12)     |       30.0s           |    23.4s                              |             31.1s               |    


### Couette flow, serial, large grid

|                         | **M1 Pro, 32 GB RAM** |  **Intel Core i9-13900K, 128 GB RAM** |  **AMD EPYC 9374F, 378 GB RAM** |
|:-----------------------:|:---------------------:|:-------------------------------------:|:-------------------------------:|
| v0.7.0 (Julia 1.11)     |       1164s           |    795s                               |             1467s               |
| v0.8.0 (Julia 1.12)     |       550s            |    418s                               |             902s                |    

### Couette flow threading speed-up
Note that this is 1) without any dynamical load-balancing (speed-up is less than expected due to density changes in the domain leading to different amount of work spent on collisions)
2) using a serial particle exchange procedure. Speed-up is compared to the "Couette flow, serial, large grid" case.

|                     | **M1 Pro, 32 GB RAM, 2/4/8 threads** | **Intel Core i9-13900K, 128 GB RAM, 2/4/8 threads** | **AMD EPYC 9374F, 378 GB RAM 2/4/8/16/32 threads** |
|:-------------------:|:------------------------------------:|:---------------------------------------------------:|:--------------------------------------------------:|
| v0.7.0 (Julia 1.11) | 2.3 / 3.8 / 5.8                      |    1.9 / 3.2 / 5.1                                  | 1.9 / 3.9 / 5.2 / 5.7 / 6.5                        |
| v0.8.0 (Julia 1.12) |       2.0 / 3.3 / 4.5                |    1.6 / 2.7 / 3.7                                  | 1.7 / 3.3 / 3.6 / 4.1 / 3.6                        |    


### Couette flow, variable-weight particles, octree merging

|                         | **M1 Pro, 32 GB RAM** |  **Intel Core i9-13900K, 128 GB RAM** |  **AMD EPYC 9374F, 378 GB RAM** |
|:-----------------------:|:---------------------:|:-------------------------------------:|:-------------------------------:|
| v0.7.10 (Julia 1.12)    |       35.2s           |    27.3s                              |             38.6s               |
| v0.8.0 (Julia 1.12)     |       28.1s           |    23.7s                              |             31.8s               |    

### 0D ionization, variable-weight particles, octree merging

|                         | **M1 Pro, 32 GB RAM** |  **Intel Core i9-13900K, 128 GB RAM** |  **AMD EPYC 9374F, 378 GB RAM** |
|:-----------------------:|:---------------------:|:-------------------------------------:|:-------------------------------:|
| v0.7.10 (Julia 1.12)    |       21.1s           |    15.1s                              |             23.0s               |
| v0.8.0 (Julia 1.12)     |       6.62s           |    4.72s                              |             6.83s               |    

### 0D ionization, variable-weight particles, NNLS merging

|                         | **M1 Pro, 32 GB RAM** |  **Intel Core i9-13900K, 128 GB RAM** |  **AMD EPYC 9374F, 378 GB RAM** |
|:-----------------------:|:---------------------:|:-------------------------------------:|:-------------------------------:|
| v0.7.10 (Julia 1.12)    |       22.9s           |    17.4s                              |             26.5s               |
| v0.8.0 (Julia 1.12)     |       4.95s           |    3.98s                              |             5.98s               |    

## Couette flow, serial, small grid

Comparison with SPARTA are provided for a single-species (argon) Couette flow test case with 50000 particles and 50 cells (averaging over 36k timesteps after t>14000). The computation is serial. Timing in Merzbild.jl providedd by [TimerOutputs.jl](https://github.com/KristofferC/TimerOutputs.jl), timing in SPARTA provided by the inbuilt timers. No surface quantities are being computed.
The input can be found in `simulations/1D/couette_benchmarking.jl`.

Merzbild.jl version 0.8.0, run with  `--check-bounds=no -O3`.

SPARTA version 20Jan2025, compiled with `-O3`.

### Intel Core i9-13900K, 128 GB RAM

Ubuntu 22.04.5, Julia version 1.12.5, SPARTA compiled with gcc version 11.4.0.

#### Merzbild.jl
```
─────────────────────────────────────────────────────────────────────────────
                                    Time                    Allocations      
                           ───────────────────────   ────────────────────────
     Tot / % measured:          23.4s /  94.9%            188MiB /   6.8%    

Section            ncalls     time    %tot     avg     alloc    %tot      avg
─────────────────────────────────────────────────────────────────────────────
sort                50.0k    7.18s   32.4%   144μs     0.00B    0.0%    0.00B
collide             2.50M    5.96s   26.9%  2.38μs     0.00B    0.0%    0.00B
convect             50.0k    4.28s   19.3%  85.6μs     0.00B    0.0%    0.00B
props compute       36.0k    3.72s   16.8%   103μs     0.00B    0.0%    0.00B
restore ordering    5.00k    949ms    4.3%   190μs     0.00B    0.0%    0.00B
I/O                     1   90.0ms    0.4%  90.0ms   10.5MiB   82.1%  10.5MiB
avg physprops       36.0k   5.96ms    0.0%   166ns     0.00B    0.0%    0.00B
sampling                1   2.44ms    0.0%  2.44ms   2.29MiB   17.9%  2.29MiB
─────────────────────────────────────────────────────────────────────────────
```

#### SPARTA
```
Loop time of 30.4417 on 1 procs for 50000 steps with 50000 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 8.0618     | 8.0618     | 8.0618     |   0.0 | 26.48
Coll    | 10.862     | 10.862     | 10.862     |   0.0 | 35.68
Sort    | 2.793      | 2.793      | 2.793      |   0.0 |  9.18
Comm    | 0.0034328  | 0.0034328  | 0.0034328  |   0.0 |  0.01
Modify  | 8.719      | 8.719      | 8.719      |   0.0 | 28.64
Output  | 0.00057459 | 0.00057459 | 0.00057459 |   0.0 |  0.00
Other   |            | 0.002392   |            |       |  0.01
```

### M1 Pro (Macbook Pro), 32 GB RAM

MacOS 15.4.1, Julia version 1.12.5, SPARTA compiled with Apple clang version 17.0.0.

#### Merzbild.jl
```
─────────────────────────────────────────────────────────────────────────────
                                    Time                    Allocations      
                           ───────────────────────   ────────────────────────
     Tot / % measured:          30.0s /  95.9%            177MiB /   7.3%    

Section            ncalls     time    %tot     avg     alloc    %tot      avg
─────────────────────────────────────────────────────────────────────────────
sort                50.0k    11.1s   38.6%   222μs     0.00B    0.0%    0.00B
collide             2.50M    6.93s   24.1%  2.77μs     0.00B    0.0%    0.00B
convect             50.0k    4.85s   16.9%  96.9μs     0.00B    0.0%    0.00B
props compute       36.0k    4.52s   15.7%   125μs     0.00B    0.0%    0.00B
restore ordering    5.00k    1.25s    4.3%   250μs     0.00B    0.0%    0.00B
I/O                     1    109ms    0.4%   109ms   10.6MiB   82.2%  10.6MiB
avg physprops       36.0k   4.58ms    0.0%   127ns     0.00B    0.0%    0.00B
sampling                1   3.65ms    0.0%  3.65ms   2.29MiB   17.8%  2.29MiB
─────────────────────────────────────────────────────────────────────────────
```

#### SPARTA
```
Loop time of 47.4825 on 1 procs for 50000 steps with 50000 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 12.001     | 12.001     | 12.001     |   0.0 | 25.28
Coll    | 13.299     | 13.299     | 13.299     |   0.0 | 28.01
Sort    | 2.8246     | 2.8246     | 2.8246     |   0.0 |  5.95
Comm    | 0.0020463  | 0.0020463  | 0.0020463  |   0.0 |  0.00
Modify  | 19.352     | 19.352     | 19.352     |   0.0 | 40.76
Output  | 0.0019748  | 0.0019748  | 0.0019748  |   0.0 |  0.00
Other   |            | 0.0009062  |            |       |  0.00
```

## Couette flow, serial, large grid
The physical parameters for this test case are the same as for the previous one, but a larger (2000 cells) grid is used, with 250 particles per cell at `t=0`.
So the number of grid cells is 40x higher than for the small grid test case, and the number of particles is 10x higher.

In addition, surface properties are also computed and averaged. The numerical setup corresponds to the `simulations/1D/couette_with_surface_quantities.jl` file
with the following parameters parameters for the `run` command:

`run(1234, 300.0, 500.0, 5e-4, 5e22, 2000, 250, 2.59e-9, 1000, 50000, 14000; do_benchmark=true)`

Setting `do_benchmark` to `true` turns off computation of the degree of particle index fragmentation.

### Intel Core i9-13900K, 128 GB RAM

Ubuntu 22.04.5, Julia version 1.12.5, SPARTA compiled with gcc version 11.4.0.

#### Merzbild.jl
```
────────────────────────────────────────────────────────────────────────────────────────
                                               Time                    Allocations      
                                      ───────────────────────   ────────────────────────
          Tot / % measured:                 418s /  99.8%            202MiB /  11.7%    

Section                       ncalls     time    %tot     avg     alloc    %tot      avg
────────────────────────────────────────────────────────────────────────────────────────
main loop                          1     418s  100.0%    418s    631KiB    2.6%   631KiB
  sort                         50.0k     122s   29.3%  2.44ms     0.00B    0.0%    0.00B
  collide                       100M     106s   25.5%  1.06μs     0.00B    0.0%    0.00B
  convect + surface compute    36.0k    69.8s   16.7%  1.94ms     0.00B    0.0%    0.00B
  props compute                36.0k    64.8s   15.5%  1.80ms     0.00B    0.0%    0.00B
  convect                      14.0k    28.5s    6.8%  2.04ms     0.00B    0.0%    0.00B
  restore ordering             5.00k    22.8s    5.5%  4.57ms     0.00B    0.0%    0.00B
  avg physprops                36.0k    220ms    0.1%  6.12μs     0.00B    0.0%    0.00B
  avg surfprops                36.0k   6.39ms    0.0%   177ns     0.00B    0.0%    0.00B
  I/O                             13   1.73ms    0.0%   133μs   3.05KiB    0.0%     240B
sampling                           1   54.4ms    0.0%  54.4ms   22.9MiB   97.4%  22.9MiB
I/O final                          2    147μs    0.0%  73.6μs      544B    0.0%     272B
────────────────────────────────────────────────────────────────────────────────────────
```

#### SPARTA
```
Loop time of 1155.34 on 1 procs for 50000 steps with 500000 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 598.55     | 598.55     | 598.55     |   0.0 | 51.81
Coll    | 262.52     | 262.52     | 262.52     |   0.0 | 22.72
Sort    | 87.103     | 87.103     | 87.103     |   0.0 |  7.54
Comm    | 0.018935   | 0.018935   | 0.018935   |   0.0 |  0.00
Modify  | 207.14     | 207.14     | 207.14     |   0.0 | 17.93
Output  | 0.002845   | 0.002845   | 0.002845   |   0.0 |  0.00
Other   |            | 0.01218    |            |       |  0.00
```

### M1 Pro (Macbook Pro), 32 GB RAM

MacOS 15.4.1, Julia version 1.12.5, SPARTA compiled with Apple clang version 17.0.0.

#### Merzbild.jl
```
────────────────────────────────────────────────────────────────────────────────────────
                                               Time                    Allocations      
                                      ───────────────────────   ────────────────────────
          Tot / % measured:                 550s /  99.9%            202MiB /  11.6%    

Section                       ncalls     time    %tot     avg     alloc    %tot      avg
────────────────────────────────────────────────────────────────────────────────────────
main loop                          1     549s  100.0%    549s    636KiB    2.6%   636KiB
  sort                         50.0k     164s   29.9%  3.28ms     0.00B    0.0%    0.00B
  collide                       100M     115s   20.9%  1.15μs     0.00B    0.0%    0.00B
  convect + surface compute    36.0k     113s   20.6%  3.14ms     0.00B    0.0%    0.00B
  props compute                36.0k    89.9s   16.4%  2.50ms     0.00B    0.0%    0.00B
  convect                      14.0k    37.8s    6.9%  2.70ms     0.00B    0.0%    0.00B
  restore ordering             5.00k    26.9s    4.9%  5.39ms     0.00B    0.0%    0.00B
  avg physprops                36.0k    180ms    0.0%  4.99μs     0.00B    0.0%    0.00B
  I/O                             13   3.38ms    0.0%   260μs   3.05KiB    0.0%     240B
  avg surfprops                36.0k   3.05ms    0.0%  84.8ns     0.00B    0.0%    0.00B
sampling                           1   24.8ms    0.0%  24.8ms   22.9MiB   97.4%  22.9MiB
I/O final                          2    272μs    0.0%   136μs      544B    0.0%     272B
────────────────────────────────────────────────────────────────────────────────────────
```

#### SPARTA
```
Loop time of 1418.29 on 1 procs for 50000 steps with 500000 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 885.53     | 885.53     | 885.53     |   0.0 | 62.44
Coll    | 274.42     | 274.42     | 274.42     |   0.0 | 19.35
Sort    | 45.409     | 45.409     | 45.409     |   0.0 |  3.20
Comm    | 0.0055656  | 0.0055656  | 0.0055656  |   0.0 |  0.00
Modify  | 212.91     | 212.91     | 212.91     |   0.0 | 15.01
Output  | 0.0032787  | 0.0032787  | 0.0032787  |   0.0 |  0.00
Other   |            | 0.006827   |            |       |  0.00
```

### AMD EPYC 9374F, 378 GB RAM
Ubuntu 24.04.3, Julia version 1.12.5.

#### Merzbild.jl
```
TODO
```

## Couette flow, variable-weight particles, octree merging
This simulation uses 200 cells, 250 initial particles per cell, octree N:2 merging down to 100 particles when particle count exceeds 150 particles in a cell.
The simulation file is `simulations/1D/couette_varweight_octree.jl`, run with `run(1234, 300.0, 500.0, 5e-4, 5e22, 200, 250, 150, 100, 2.59e-9, 1000, n_t, 14000; debug=false)`.

### Intel Core i9-13900K, 128 GB RAM
Ubuntu 22.04.5, Julia version 1.12.5.
```
──────────────────────────────────────────────────────────────────────────────────────
                                             Time                    Allocations      
                                    ───────────────────────   ────────────────────────
         Tot / % measured:               23.7s /  95.4%            144MiB /   2.1%    

Section                     ncalls     time    %tot     avg     alloc    %tot      avg
──────────────────────────────────────────────────────────────────────────────────────
collide                      10.0M    10.1s   44.8%  1.01μs     0.00B    0.0%    0.00B
sort                         50.0k    3.86s   17.0%  77.2μs     0.00B    0.0%    0.00B
merge                         569k    2.80s   12.4%  4.92μs     0.00B    0.0%    0.00B
props compute                36.0k    2.11s    9.3%  58.7μs     0.00B    0.0%    0.00B
convect + surface compute    36.0k    1.87s    8.2%  51.9μs     0.00B    0.0%    0.00B
restore ordering             5.00k    729ms    3.2%   146μs     0.00B    0.0%    0.00B
convect                      14.0k    724ms    3.2%  51.7μs     0.00B    0.0%    0.00B
squash                       50.0k    415ms    1.8%  8.30μs     0.00B    0.0%    0.00B
sampling                         1   2.63ms    0.0%  2.63ms   3.05MiB   99.6%  3.05MiB
I/O                             52   2.18ms    0.0%  42.0μs   12.2KiB    0.4%     241B
merge (t=0)                    200   1.29ms    0.0%  6.46μs     0.00B    0.0%    0.00B
squash (t=0)                     1   21.1μs    0.0%  21.1μs     0.00B    0.0%    0.00B
──────────────────────────────────────────────────────────────────────────────────────
```

### M1 Pro (Macbook Pro), 32 GB RAM
```
──────────────────────────────────────────────────────────────────────────────────────
                                             Time                    Allocations      
                                    ───────────────────────   ────────────────────────
         Tot / % measured:               28.1s /  95.8%            145MiB /   2.1%    

Section                     ncalls     time    %tot     avg     alloc    %tot      avg
──────────────────────────────────────────────────────────────────────────────────────
collide                      10.0M    11.7s   43.3%  1.17μs     0.00B    0.0%    0.00B
sort                         50.0k    5.33s   19.8%   107μs     0.00B    0.0%    0.00B
merge                         569k    3.28s   12.2%  5.77μs     0.00B    0.0%    0.00B
props compute                36.0k    2.22s    8.2%  61.6μs     0.00B    0.0%    0.00B
convect + surface compute    36.0k    2.19s    8.1%  60.9μs     0.00B    0.0%    0.00B
restore ordering             5.00k    915ms    3.4%   183μs     0.00B    0.0%    0.00B
convect                      14.0k    713ms    2.6%  50.9μs     0.00B    0.0%    0.00B
squash                       50.0k    628ms    2.3%  12.6μs     0.00B    0.0%    0.00B
I/O                             52   3.00ms    0.0%  57.7μs   12.2KiB    0.4%     241B
sampling                         1   2.65ms    0.0%  2.65ms   3.05MiB   99.6%  3.05MiB
merge (t=0)                    200   1.60ms    0.0%  7.98μs     0.00B    0.0%    0.00B
squash (t=0)                     1   26.4μs    0.0%  26.4μs     0.00B    0.0%    0.00B
──────────────────────────────────────────────────────────────────────────────────────
```

### AMD EPYC 9374F, 378 GB RAM

TODO

## 0D ionization, variable-weight particles, octree merging
This simulation can be found under `simulations/0D/0D_ionization_1neutralspecies_es.jl`. Lisbon IST data was used for the cross-sections,
the run is `run(1234, 400.0, 500000, 500, 350, cs_n_e_filepath; merging_bin_split=OctreeBinMidSplit, adds=0, do_es=true)`.
Octree N:2 merging is used with a threshold of 500 particles and a target of 350 particles, along with event splitting.

### Intel Core i9-13900K, 128 GB RAM
```
──────────────────────────────────────────────────────────────────────────
                                 Time                    Allocations      
                        ───────────────────────   ────────────────────────
   Tot / % measured:         4.72s /  82.0%            244MiB /  41.2%    

Section         ncalls     time    %tot     avg     alloc    %tot      avg
──────────────────────────────────────────────────────────────────────────
props             500k    1.63s   42.1%  3.26μs     0.00B    0.0%    0.00B
coll n-e ES       500k    725ms   18.7%  1.45μs     0.00B    0.0%    0.00B
acc e             500k    689ms   17.8%  1.38μs     0.00B    0.0%    0.00B
I/O               500k    668ms   17.3%  1.34μs   99.2MiB   98.5%     208B
merge e            697   98.9ms    2.6%   142μs     0.00B    0.0%    0.00B
merge n          6.10k   34.2ms    0.9%  5.61μs     0.00B    0.0%    0.00B
coll n-n          500k   18.8ms    0.5%  37.7ns     0.00B    0.0%    0.00B
merge i          1.80k   4.09ms    0.1%  2.27μs     0.00B    0.0%    0.00B
merge e (t=0)        1   1.42ms    0.0%  1.42ms   1.48MiB    1.5%  1.48MiB
merge i (t=0)        1    148μs    0.0%   148μs     0.00B    0.0%    0.00B
merge n (t=0)        1    123μs    0.0%   123μs     0.00B    0.0%    0.00B
──────────────────────────────────────────────────────────────────────────
```

### M1 Pro (Macbook Pro), 32 GB RAM

```
──────────────────────────────────────────────────────────────────────────
                                 Time                    Allocations      
                        ───────────────────────   ────────────────────────
   Tot / % measured:         6.62s /  86.8%            245MiB /  41.1%    

Section         ncalls     time    %tot     avg     alloc    %tot      avg
──────────────────────────────────────────────────────────────────────────
props             500k    2.48s   43.0%  4.95μs     0.00B    0.0%    0.00B
acc e             500k    1.09s   19.0%  2.19μs     0.00B    0.0%    0.00B
I/O               500k    1.01s   17.6%  2.02μs   99.2MiB   98.5%     208B
coll n-e ES       500k    973ms   16.9%  1.95μs     0.00B    0.0%    0.00B
merge e            697    132ms    2.3%   190μs     0.00B    0.0%    0.00B
merge n          6.10k   40.4ms    0.7%  6.62μs     0.00B    0.0%    0.00B
coll n-n          500k   17.5ms    0.3%  35.0ns     0.00B    0.0%    0.00B
merge i          1.80k   5.45ms    0.1%  3.02μs     0.00B    0.0%    0.00B
merge e (t=0)        1   1.70ms    0.0%  1.70ms   1.50MiB    1.5%  1.50MiB
merge i (t=0)        1    239μs    0.0%   239μs     0.00B    0.0%    0.00B
merge n (t=0)        1    203μs    0.0%   203μs     0.00B    0.0%    0.00B
──────────────────────────────────────────────────────────────────────────
```

### AMD EPYC 9374F, 378 GB RAM

TODO

## 0D ionization, variable-weight particles, NNLS merging
This simulation can be found under `simulations/0D/0D_ionization_1neutralspecies_nnls_es.jl`. Lisbon IST data was used for the cross-sections
This uses event splitting and non-rate-preserving NNLS merging conserving all mixed-order moments up to total order 6.
Event splitting is used.
`run(1234, 400.0, 500000, 6, 95, 88, cs_n_e_filepath; adds=0, rate_preserving=:off, do_es=true)`

### Intel Core i9-13900K, 128 GB RAM
```
────────────────────────────────────────────────────────────────────────
                               Time                    Allocations      
                      ───────────────────────   ────────────────────────
  Tot / % measured:        3.98s /  76.8%            802MiB /  80.4%    

Section       ncalls     time    %tot     avg     alloc    %tot      avg
────────────────────────────────────────────────────────────────────────
NNLSmerge e    3.04k    987ms   32.3%   325μs    391MiB   60.6%   132KiB
merge e t=0        1    940ms   30.8%   940ms    126MiB   19.6%   126MiB
I/O             500k    686ms   22.5%  1.37μs    114MiB   17.7%     240B
props           500k    261ms    8.6%   523ns     0.00B    0.0%    0.00B
coll n-e ES     500k   60.5ms    2.0%   121ns     0.00B    0.0%    0.00B
NNLSinit           2   54.7ms    1.8%  27.4ms   13.2MiB    2.0%  6.59MiB
acc e           500k   46.3ms    1.5%  92.6ns     0.00B    0.0%    0.00B
coll n-n        500k   16.5ms    0.5%  32.9ns     0.00B    0.0%    0.00B
merge n          390   2.53ms    0.1%  6.49μs     0.00B    0.0%    0.00B
merge i          184    489μs    0.0%  2.66μs     0.00B    0.0%    0.00B
merge i t=0        1    167μs    0.0%   167μs     0.00B    0.0%    0.00B
merge n t=0        1    127μs    0.0%   127μs     0.00B    0.0%    0.00B
NNLSinit RP        1   16.6μs    0.0%  16.6μs    142KiB    0.0%   142KiB
────────────────────────────────────────────────────────────────────────
```

### M1 Pro (Macbook Pro), 32 GB RAM

```
────────────────────────────────────────────────────────────────────────
                               Time                    Allocations      
                      ───────────────────────   ────────────────────────
  Tot / % measured:        4.95s /  79.6%            809MiB /  80.5%    

Section       ncalls     time    %tot     avg     alloc    %tot      avg
────────────────────────────────────────────────────────────────────────
NNLSmerge e    3.03k    1.46s   37.1%   483μs    396MiB   60.8%   134KiB
I/O             500k    1.08s   27.5%  2.16μs    114MiB   17.6%     240B
merge e t=0        1    809ms   20.5%   809ms    127MiB   19.5%   127MiB
props           500k    384ms    9.8%   769ns     0.00B    0.0%    0.00B
coll n-e ES     500k   76.0ms    1.9%   152ns     0.00B    0.0%    0.00B
NNLSinit           2   57.4ms    1.5%  28.7ms   13.3MiB    2.0%  6.64MiB
acc e           500k   51.5ms    1.3%   103ns     0.00B    0.0%    0.00B
coll n-n        500k   12.1ms    0.3%  24.3ns     0.00B    0.0%    0.00B
merge n          393   3.12ms    0.1%  7.95μs     0.00B    0.0%    0.00B
merge i          183    714μs    0.0%  3.90μs     0.00B    0.0%    0.00B
merge i t=0        1    246μs    0.0%   246μs     0.00B    0.0%    0.00B
merge n t=0        1    237μs    0.0%   237μs     0.00B    0.0%    0.00B
NNLSinit RP        1   12.8μs    0.0%  12.8μs    142KiB    0.0%   142KiB
────────────────────────────────────────────────────────────────────────
```

### AMD EPYC 9374F, 378 GB RAM

TODO