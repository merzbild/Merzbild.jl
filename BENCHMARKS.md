# Merzbild.jl benchmarks

Various benchmarks and comparisons to other open-source codes are provided here for reference test cases.

## Couette flow, serial

Comparison with SPARTA are provided for a single-species (argon) Couette flow test case with 50000 particles and 50 cells (averaging over 36k timesteps after t>14000). The computation is serial. Timing in Merzbild.jl providedd by [TimerOutputs.jl](https://github.com/KristofferC/TimerOutputs.jl), timing in SPARTA provided by the inbuilt timers. No surface quantities are being computed.
The input can be found in `simulations/1D/couette_benchmarking.jl`.

Merzbild.jl version 0.6.6, run with  `--check-bounds=no -O3`.

SPARTA version 4Sep2024, compiled with `-O3`.

### Intel Core i9-13900K, 128 GB RAM

Ubuntu 22.04.5, Julia version 1.11.2, SPARTA compiled with gcc version 11.4.0.

#### Merzbild.jl
```
──────────────────────────────────────────────────────────────────────────
                                 Time                    Allocations      
                        ───────────────────────   ────────────────────────
   Tot / % measured:         28.9s /  97.4%            127MiB /   9.8%    

Section         ncalls     time    %tot     avg     alloc    %tot      avg
──────────────────────────────────────────────────────────────────────────
sort             50.0k    10.4s   36.9%   208μs     0.00B    0.0%    0.00B
convect          50.0k    6.66s   23.7%   133μs     0.00B    0.0%    0.00B
collide          2.50M    5.91s   21.0%  2.37μs     0.00B    0.0%    0.00B
props compute    36.0k    5.08s   18.0%   141μs     0.00B    0.0%    0.00B
I/O                  1   85.7ms    0.3%  85.7ms   9.29MiB   75.3%  9.29MiB
sampling             1   13.4ms    0.0%  13.4ms   3.05MiB   24.7%  3.05MiB
avg physprops    36.0k   6.17ms    0.0%   172ns     0.00B    0.0%    0.00B
──────────────────────────────────────────────────────────────────────────
```

#### SPARTA

```
Loop time of 31.1462 on 1 procs for 50000 steps with 50000 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 8.1848     | 8.1848     | 8.1848     |   0.0 | 26.28
Coll    | 11.155     | 11.155     | 11.155     |   0.0 | 35.81
Sort    | 2.8435     | 2.8435     | 2.8435     |   0.0 |  9.13
Comm    | 0.0032952  | 0.0032952  | 0.0032952  |   0.0 |  0.01
Modify  | 8.9578     | 8.9578     | 8.9578     |   0.0 | 28.76
Output  | 0.00046062 | 0.00046062 | 0.00046062 |   0.0 |  0.00
Other   |            | 0.001451   |            |       |  0.00
```

### M1 Pro (Macbook Pro), 32 GB RAM

MacOS 15.4.1, Julia version 1.11.2, SPARTA compiled with Apple clang version 14.0.0.

#### Merzbild.jl
```
──────────────────────────────────────────────────────────────────────────
                                 Time                    Allocations      
                        ───────────────────────   ────────────────────────
   Tot / % measured:         33.6s /  97.8%            125MiB /   9.9%    

Section         ncalls     time    %tot     avg     alloc    %tot      avg
──────────────────────────────────────────────────────────────────────────
sort             50.0k    11.7s   35.5%   233μs     0.00B    0.0%    0.00B
collide          2.50M    8.00s   24.3%  3.20μs     0.00B    0.0%    0.00B
convect          50.0k    7.51s   22.8%   150μs     0.00B    0.0%    0.00B
props compute    36.0k    5.62s   17.1%   156μs     0.00B    0.0%    0.00B
I/O                  1   91.7ms    0.3%  91.7ms   9.30MiB   75.3%  9.30MiB
avg physprops    36.0k   4.99ms    0.0%   139ns     0.00B    0.0%    0.00B
sampling             1   2.74ms    0.0%  2.74ms   3.05MiB   24.7%  3.05MiB
──────────────────────────────────────────────────────────────────────────
```

#### SPARTA

```
Loop time of 51.1924 on 1 procs for 50000 steps with 50000 particles

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Move    | 13.622     | 13.622     | 13.622     |   0.0 | 26.61
Coll    | 13.132     | 13.132     | 13.132     |   0.0 | 25.65
Sort    | 2.8252     | 2.8252     | 2.8252     |   0.0 |  5.52
Comm    | 0.002234   | 0.002234   | 0.002234   |   0.0 |  0.00
Modify  | 21.607     | 21.607     | 21.607     |   0.0 | 42.21
Output  | 0.0025806  | 0.0025806  | 0.0025806  |   0.0 |  0.01
Other   |            | 0.0009735  |            |       |  0.00
```