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
   Tot / % measured:         28.1s /  97.3%            126MiB /   9.8%    

Section         ncalls     time    %tot     avg     alloc    %tot      avg
──────────────────────────────────────────────────────────────────────────
sort             50.0k    10.4s   38.0%   208μs     0.00B    0.0%    0.00B
convect          50.0k    5.92s   21.6%   118μs     0.00B    0.0%    0.00B
collide          2.50M    5.89s   21.5%  2.36μs     0.00B    0.0%    0.00B
props compute    36.0k    5.07s   18.5%   141μs     0.00B    0.0%    0.00B
I/O                  1   84.8ms    0.3%  84.8ms   9.32MiB   75.3%  9.32MiB
avg physprops    36.0k   5.56ms    0.0%   154ns     0.00B    0.0%    0.00B
sampling             1   2.44ms    0.0%  2.44ms   3.05MiB   24.7%  3.05MiB
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
   Tot / % measured:         33.4s /  97.6%            125MiB /   9.9%    

Section         ncalls     time    %tot     avg     alloc    %tot      avg
──────────────────────────────────────────────────────────────────────────
sort             50.0k    11.6s   35.7%   233μs     0.00B    0.0%    0.00B
collide          2.50M    8.02s   24.6%  3.21μs     0.00B    0.0%    0.00B
convect          50.0k    7.29s   22.4%   146μs     0.00B    0.0%    0.00B
props compute    36.0k    5.58s   17.1%   155μs     0.00B    0.0%    0.00B
I/O                  1   95.2ms    0.3%  95.2ms   9.30MiB   75.3%  9.30MiB
avg physprops    36.0k   5.19ms    0.0%   144ns     0.00B    0.0%    0.00B
sampling             1   2.66ms    0.0%  2.66ms   3.05MiB   24.7%  3.05MiB
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
         Tot / % measured:                816s /  99.4%            334MiB /  11.8%    

Section                     ncalls     time    %tot     avg     alloc    %tot      avg
──────────────────────────────────────────────────────────────────────────────────────
sort                         50.0k     271s   33.5%  5.43ms     0.00B    0.0%    0.00B
convect + surface compute    36.0k     217s   26.8%  6.03ms   8.79MiB   22.4%     256B
props compute                36.0k     134s   16.5%  3.72ms     0.00B    0.0%    0.00B
collide                       100M     115s   14.2%  1.15μs     0.00B    0.0%    0.00B
convect                      14.0k    73.1s    9.0%  5.22ms     0.00B    0.0%    0.00B
avg physprops                36.0k    207ms    0.0%  5.75μs     0.00B    0.0%    0.00B
sampling                         1   27.8ms    0.0%  27.8ms   30.5MiB   77.6%  30.5MiB
avg surfprops                36.0k   12.4ms    0.0%   345ns     0.00B    0.0%    0.00B
I/O                             52   8.21ms    0.0%   158μs   12.2KiB    0.0%     241B
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
         Tot / % measured:               1169s /  99.6%            197MiB /  20.0%    

Section                     ncalls     time    %tot     avg     alloc    %tot      avg
──────────────────────────────────────────────────────────────────────────────────────
convect + surface compute    36.0k     362s   31.1%  10.1ms   8.79MiB   22.4%     256B
sort                         50.0k     336s   28.9%  6.72ms     0.00B    0.0%    0.00B
props compute                36.0k     183s   15.7%  5.08ms     0.00B    0.0%    0.00B
collide                       100M     174s   15.0%  1.74μs     0.00B    0.0%    0.00B
convect                      14.0k     108s    9.3%  7.74ms     0.00B    0.0%    0.00B
avg physprops                36.0k    199ms    0.0%  5.52μs     0.00B    0.0%    0.00B
sampling                         1   34.2ms    0.0%  34.2ms   30.5MiB   77.6%  30.5MiB
avg surfprops                36.0k   11.0ms    0.0%   306ns     0.00B    0.0%    0.00B
I/O                             15   2.14ms    0.0%   143μs   3.58KiB    0.0%     244B
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

## Couette flow, multi-threaded, large grid
The numerical and physical parameters are the same as for the serial large grid case (2000 cells, 250 particles per cell at `t=0`).
The simulation file is `simulations/1D/couette_multithreaded.jl`.

### Intel Core i9-13900K, 128 GB RAM

Ubuntu 22.04.5, Julia version 1.11.2.
Shown is the speed-up compared to a serial execution on the same computer (see above). `DLB` denotes dynamic load balancing (currently not used).

|                               | **2 cores** |  **4 cores** |  **8 cores** |
|:-----------------------------:|:-----------:|:------------:|:------------:|
| `n_chunks=n_threads`, no DLB  |     1.89    |   3.31       |    5.16      |         



### M1 Pro (Macbook Pro), 32 GB RAM
MacOS 15.4.1, Julia version 1.11.2.
Shown is the speed-up compared to a serial execution on the same computer (see above). `DLB` denotes dynamic load balancing (currently not used).

|                               | **2 cores** |  **4 cores** |  **8 cores** |
|:-----------------------------:|:-----------:|:------------:|:------------:|
| `n_chunks=n_threads`, no DLB  |     2.26    |   3.78       |    5.7      |      