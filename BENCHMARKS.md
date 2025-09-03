# Merzbild.jl benchmarks

Various benchmarks and comparisons to other open-source codes are provided here for reference test cases.

## Couette flow, serial, small grid

Comparison with SPARTA are provided for a single-species (argon) Couette flow test case with 50000 particles and 50 cells (averaging over 36k timesteps after t>14000). The computation is serial. Timing in Merzbild.jl providedd by [TimerOutputs.jl](https://github.com/KristofferC/TimerOutputs.jl), timing in SPARTA provided by the inbuilt timers. No surface quantities are being computed.
The input can be found in `simulations/1D/couette_benchmarking.jl`.

Merzbild.jl version 0.7.0, run with  `--check-bounds=no -O3`.

SPARTA version 20Jan2025, compiled with `-O3`.

### Intel Core i9-13900K, 128 GB RAM

Ubuntu 22.04.5, Julia version 1.11.2, SPARTA compiled with gcc version 11.4.0.

#### Merzbild.jl
```
──────────────────────────────────────────────────────────────────────────
                                 Time                    Allocations      
                        ───────────────────────   ────────────────────────
   Tot / % measured:         28.4s /  97.3%            126MiB /   9.8%    

Section         ncalls     time    %tot     avg     alloc    %tot      avg
──────────────────────────────────────────────────────────────────────────
sort             50.0k    10.5s   37.9%   209μs     0.00B    0.0%    0.00B
convect          50.0k    6.29s   22.8%   126μs     0.00B    0.0%    0.00B
collide          2.50M    5.90s   21.3%  2.36μs     0.00B    0.0%    0.00B
props compute    36.0k    4.88s   17.7%   136μs     0.00B    0.0%    0.00B
I/O                  1   85.4ms    0.3%  85.4ms   9.32MiB   75.3%  9.32MiB
avg physprops    36.0k   5.54ms    0.0%   154ns     0.00B    0.0%    0.00B
sampling             1   2.52ms    0.0%  2.52ms   3.05MiB   24.7%  3.05MiB
──────────────────────────────────────────────────────────────────────────
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

MacOS 15.4.1, Julia version 1.11.2, SPARTA compiled with Apple clang version 17.0.0.

#### Merzbild.jl
```
──────────────────────────────────────────────────────────────────────────
                                 Time                    Allocations      
                        ───────────────────────   ────────────────────────
   Tot / % measured:         33.5s /  97.7%            125MiB /   9.9%    

Section         ncalls     time    %tot     avg     alloc    %tot      avg
──────────────────────────────────────────────────────────────────────────
sort             50.0k    11.6s   35.5%   232μs     0.00B    0.0%    0.00B
collide          2.50M    8.01s   24.4%  3.20μs     0.00B    0.0%    0.00B
convect          50.0k    7.44s   22.7%   149μs     0.00B    0.0%    0.00B
props compute    36.0k    5.60s   17.1%   156μs     0.00B    0.0%    0.00B
I/O                  1   91.8ms    0.3%  91.8ms   9.30MiB   75.3%  9.30MiB
avg physprops    36.0k   5.33ms    0.0%   148ns     0.00B    0.0%    0.00B
sampling             1   2.78ms    0.0%  2.78ms   3.05MiB   24.7%  3.05MiB
──────────────────────────────────────────────────────────────────────────
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

Ubuntu 22.04.5, Julia version 1.11.2, SPARTA compiled with gcc version 11.4.0.

#### Merzbild.jl
```
──────────────────────────────────────────────────────────────────────────────────────
                                             Time                    Allocations      
                                    ───────────────────────   ────────────────────────
         Tot / % measured:                795s /  99.5%            192MiB /  17.1%    

Section                     ncalls     time    %tot     avg     alloc    %tot      avg
──────────────────────────────────────────────────────────────────────────────────────
sort                         50.0k     269s   34.0%  5.37ms     0.00B    0.0%    0.00B
convect + surface compute    36.0k     209s   26.4%  5.80ms   2.20MiB    6.7%    64.0B
props compute                36.0k     131s   16.6%  3.63ms     0.00B    0.0%    0.00B
collide                       100M     114s   14.4%  1.14μs     0.00B    0.0%    0.00B
convect                      14.0k    68.1s    8.6%  4.86ms     0.00B    0.0%    0.00B
avg physprops                36.0k    196ms    0.0%  5.43μs     0.00B    0.0%    0.00B
sampling                         1   35.1ms    0.0%  35.1ms   30.5MiB   93.3%  30.5MiB
avg surfprops                36.0k   12.6ms    0.0%   351ns     0.00B    0.0%    0.00B
I/O                             15   2.26ms    0.0%   151μs   3.58KiB    0.0%     244B
──────────────────────────────────────────────────────────────────────────────────────
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

MacOS 15.4.1, Julia version 1.11.2, SPARTA compiled with Apple clang version 17.0.0.

#### Merzbild.jl
```
──────────────────────────────────────────────────────────────────────────────────────
                                             Time                    Allocations      
                                    ───────────────────────   ────────────────────────
         Tot / % measured:               1164s /  99.5%            190MiB /  17.2%    

Section                     ncalls     time    %tot     avg     alloc    %tot      avg
──────────────────────────────────────────────────────────────────────────────────────
convect + surface compute    36.0k     359s   31.0%  10.0ms   2.20MiB    6.7%    64.0B
sort                         50.0k     337s   29.0%  6.73ms     0.00B    0.0%    0.00B
props compute                36.0k     181s   15.6%  5.02ms     0.00B    0.0%    0.00B
collide                       100M     171s   14.8%  1.71μs     0.00B    0.0%    0.00B
convect                      14.0k     110s    9.5%  7.87ms     0.00B    0.0%    0.00B
avg physprops                36.0k    198ms    0.0%  5.51μs     0.00B    0.0%    0.00B
sampling                         1   26.1ms    0.0%  26.1ms   30.5MiB   93.3%  30.5MiB
avg surfprops                36.0k   14.8ms    0.0%   410ns     0.00B    0.0%    0.00B
I/O                             15   1.77ms    0.0%   118μs   3.58KiB    0.0%     244B
──────────────────────────────────────────────────────────────────────────────────────
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
Ubuntu 24.04.3, Julia version 1.11.6.

#### Merzbild.jl
```
──────────────────────────────────────────────────────────────────────────────────────
                                             Time                    Allocations      
                                    ───────────────────────   ────────────────────────
         Tot / % measured:               1467s /  99.5%            192MiB /  17.0%    

Section                     ncalls     time    %tot     avg     alloc    %tot      avg
──────────────────────────────────────────────────────────────────────────────────────
convect + surface compute    36.0k     584s   40.0%  16.2ms   2.20MiB    6.7%    64.0B
sort                         50.0k     369s   25.3%  7.39ms     0.00B    0.0%    0.00B
collide                       100M     173s   11.8%  1.73μs     0.00B    0.0%    0.00B
convect                      14.0k     170s   11.6%  12.1ms     0.00B    0.0%    0.00B
props compute                36.0k     164s   11.2%  4.54ms     0.00B    0.0%    0.00B
avg physprops                36.0k    306ms    0.0%  8.51μs     0.00B    0.0%    0.00B
sampling                         1   31.4ms    0.0%  31.4ms   30.5MiB   93.3%  30.5MiB
avg surfprops                36.0k   21.8ms    0.0%   605ns     0.00B    0.0%    0.00B
I/O                             15   7.44ms    0.0%   496μs   3.58KiB    0.0%     244B
──────────────────────────────────────────────────────────────────────────────────────
```

## Couette flow, multi-threaded, large grid
The numerical and physical parameters are the same as for the serial large grid case (2000 cells, 250 particles per cell at `t=0`).
The simulation file is `simulations/1D/couette_multithreaded.jl`.

### Intel Core i9-13900K, 128 GB RAM

Ubuntu 22.04.5, Julia version 1.11.2.
Shown is the speed-up compared to a serial execution on the same computer (see above). `DLB` denotes dynamic load balancing (currently not used).

|                               | **2 cores** |  **4 cores** |  **8 cores** |
|:-----------------------------:|:-----------:|:------------:|:------------:|
| `n_chunks=n_threads`, no DLB  |   1.85      |    3.24      |     5.06     |         



### M1 Pro (Macbook Pro), 32 GB RAM
MacOS 15.4.1, Julia version 1.11.2.
Shown is the speed-up compared to a serial execution on the same computer (see above). `DLB` denotes dynamic load balancing (currently not used).

|                               | **2 cores** |  **4 cores** |  **8 cores** |
|:-----------------------------:|:-----------:|:------------:|:------------:|
| `n_chunks=n_threads`, no DLB  |    2.27     |    3.75      |     5.76     |      



### AMD EPYC 9374F, 378 GB RAM
Ubuntu 24.04.3, Julia version 1.11.6.
Shown is the speed-up compared to a serial execution on the same computer (see above). `DLB` denotes dynamic load balancing (currently not used).

|                               | **2 cores** |  **4 cores** |  **8 cores** |  **16 cores** | **32 cores ** |
|:-----------------------------:|:-----------:|:------------:|:------------:|:-------------:|:-------------:|
| `n_chunks=n_threads`, no DLB  |   1.88      |    3.68      |     5.24     |      5.70     |     6.57      |
