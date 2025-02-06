# Merzbild.jl benchmarks

Various benchmarks and comparisons to other open-source codes are provided here for reference test cases.

## Couette flow, serial

Comparison with SPARTA are provided for a single-species (argon) Couette flow test case with 50000 particles and 50 cells (averaging over 36k timesteps after t>14000). The computation is serial. Timing in Merzbild.jl providedd by [TimerOutputs.jl](https://github.com/KristofferC/TimerOutputs.jl), timing in SPARTA provided by the inbuilt timers.

Merzbild.jl version 0.6.2, run with  `--check-bounds=no -O3`.

SPARTA version 4Sep2024, compiled with `-O3`.

### Intel Core i9-13900K, 128 GB RAM

Ubuntu 22.04.5, Julia version 1.11.2, gcc version 11.4.0.

#### Merzbild.jl
```
 ──────────────────────────────────────────────────────────────────────────
                                  Time                    Allocations      
                         ───────────────────────   ────────────────────────
    Tot / % measured:         28.9s /  95.7%            698MiB /   0.4%    

 Section         ncalls     time    %tot     avg     alloc    %tot      avg
 ──────────────────────────────────────────────────────────────────────────
 sort             50.0k    10.3s   37.2%   206μs     0.00B    0.0%    0.00B
 convect          50.0k    6.48s   23.4%   130μs     0.00B    0.0%    0.00B
 collide          2.50M    5.93s   21.5%  2.37μs     0.00B    0.0%    0.00B
 props compute    36.0k    4.95s   17.9%   138μs     0.00B    0.0%    0.00B
 sampling             1   2.33ms    0.0%  2.33ms   3.05MiB  100.0%  3.05MiB
 I/O                  1    427μs    0.0%   427μs      240B    0.0%     240B
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

MacOS 12.7.5, Julia version 1.11.2, Apple clang version 14.0.0.

#### Merzbild.jl
```
──────────────────────────────────────────────────────────────────────────
                                 Time                    Allocations      
                        ───────────────────────   ────────────────────────
   Tot / % measured:         33.8s /  96.3%            693MiB /   0.4%    

Section         ncalls     time    %tot     avg     alloc    %tot      avg
──────────────────────────────────────────────────────────────────────────
sort             50.0k    11.6s   35.7%   233μs     0.00B    0.0%    0.00B
collide          2.50M    8.57s   26.3%  3.43μs     0.00B    0.0%    0.00B
convect          50.0k    6.86s   21.0%   137μs     0.00B    0.0%    0.00B
props compute    36.0k    5.51s   16.9%   153μs     0.00B    0.0%    0.00B
sampling             1   2.64ms    0.0%  2.64ms   3.05MiB  100.0%  3.05MiB
I/O                  1    139μs    0.0%   139μs      240B    0.0%     240B
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