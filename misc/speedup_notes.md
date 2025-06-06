12.6.24
Octree merging profiling: loop3 is the slowest (bin shifting)
Checking condition for best refinement is also slow

```
 Section          ncalls     time    %tot     avg     alloc    %tot      avg
 ───────────────────────────────────────────────────────────────────────────
 merge                50    1.26s   82.1%  25.3ms    554KiB    2.3%  11.1KiB
   condition       53.6k    628ms   40.8%  11.7μs     0.00B    0.0%    0.00B
   split           53.5k    588ms   38.2%  11.0μs   6.61KiB    0.0%    0.13B
     loop 3        53.5k    507ms   32.9%  9.48μs     0.00B    0.0%    0.00B
     loop 1        53.5k   20.0ms    1.3%   375ns     0.00B    0.0%    0.00B
     loop 5        53.5k   7.34ms    0.5%   137ns     0.00B    0.0%    0.00B
     ifbounds      53.5k   4.10ms    0.3%  76.6ns     0.00B    0.0%    0.00B
     loop 6        53.5k   3.19ms    0.2%  59.7ns     0.00B    0.0%    0.00B
     loop 2        53.5k   2.62ms    0.2%  49.0ns     0.00B    0.0%    0.00B
     loop 4        53.5k   1.63ms    0.1%  30.5ns     0.00B    0.0%    0.00B
     ifsplit       53.5k   1.63ms    0.1%  30.4ns     0.00B    0.0%    0.00B
     bounds1       53.5k   1.23ms    0.1%  23.0ns     0.00B    0.0%    0.00B
   binprops         260k   16.1ms    1.0%  61.8ns     0.00B    0.0%    0.00B
   newparticles       50   8.03ms    0.5%   161μs     0.00B    0.0%    0.00B
   init               50   1.68ms    0.1%  33.7μs     0.00B    0.0%    0.00B
   resize             50   37.9μs    0.0%   758ns    542KiB    2.2%  10.8KiB
   clear              50   7.71μs    0.0%   154ns     0.00B    0.0%    0.00B
 I/O                 500    134ms    8.7%   268μs   16.0MiB   67.1%  32.7KiB
 props compute       500    118ms    7.7%   237μs     0.00B    0.0%    0.00B
 collide             500   23.0ms    1.5%  45.9μs   7.28MiB   30.6%  14.9KiB
 ───────────────────────────────────────────────────────────────────────────```

 With new bin positioning/without shifting of bins:
```
 Section          ncalls     time    %tot     avg     alloc    %tot      avg
 ───────────────────────────────────────────────────────────────────────────
 merge                51    765ms   73.5%  15.0ms    553KiB    2.3%  10.8KiB
   condition       54.8k    635ms   61.1%  11.6μs     0.00B    0.0%    0.00B
   split           54.7k   81.6ms    7.8%  1.49μs   5.88KiB    0.0%    0.11B
     loop 1        54.7k   18.9ms    1.8%   345ns     0.00B    0.0%    0.00B
     ifbounds      54.7k   7.00ms    0.7%   128ns     0.00B    0.0%    0.00B
     loop 5        54.7k   6.82ms    0.7%   125ns     0.00B    0.0%    0.00B
     loop 6        54.7k   3.31ms    0.3%  60.5ns     0.00B    0.0%    0.00B
     loop 2        54.7k   2.66ms    0.3%  48.6ns     0.00B    0.0%    0.00B
     loop 4        54.7k   2.09ms    0.2%  38.2ns     0.00B    0.0%    0.00B
     ifsplit       54.7k   1.63ms    0.2%  29.8ns     0.00B    0.0%    0.00B
     bounds1       54.7k   1.20ms    0.1%  21.9ns     0.00B    0.0%    0.00B
   binprops         266k   15.9ms    1.5%  59.7ns     0.00B    0.0%    0.00B
   newparticles       51   8.24ms    0.8%   162μs     0.00B    0.0%    0.00B
   init               51   1.73ms    0.2%  33.9μs     0.00B    0.0%    0.00B
   resize             51   23.8μs    0.0%   467ns    542KiB    2.2%  10.6KiB
   clear              51   5.54μs    0.0%   109ns     0.00B    0.0%    0.00B
 I/O                 500    134ms   12.9%   269μs   16.0MiB   67.2%  32.7KiB
 props compute       500    119ms   11.4%   237μs     0.00B    0.0%    0.00B
 collide             500   22.7ms    2.2%  45.4μs   7.27MiB   30.6%  14.9KiB
 ───────────────────────────────────────────────────────────────────────────```

 condition       54.8k    143ms   26.1%  2.61μs     0.00B    0.0%    0.00B

NNLS merging
```
Section               ncalls     time    %tot     avg     alloc    %tot      avg
 ────────────────────────────────────────────────────────────────────────────────
 NNLSmerge: 1st time        1    810ms   93.1%   810ms   25.7MiB   77.1%  25.7MiB
   solve                    1    783ms   90.0%   783ms     16.0B    0.0%    16.0B
   LHS/RHS                  1   15.9ms    1.8%  15.9ms      816B    0.0%     816B
   load                     1   7.08ms    0.8%  7.08ms   8.77MiB   26.3%  8.77MiB
   scaling                  1   3.49ms    0.4%  3.49ms   8.48MiB   25.4%  8.48MiB
   residual                 1    362μs    0.0%   362μs      848B    0.0%     848B
   lhs alloc                1   78.7μs    0.0%  78.7μs   8.47MiB   25.4%  8.47MiB
   LHS+                     1   29.6μs    0.0%  29.6μs     0.00B    0.0%    0.00B
 NNLSmerge                 11   59.9ms    6.9%  5.44ms   7.21MiB   21.6%   672KiB
   solve                   11   52.1ms    6.0%  4.74ms      176B    0.0%    16.0B
   LHS/RHS                 11   4.51ms    0.5%   410μs   8.77KiB    0.0%     816B
   residual                11   1.98ms    0.2%   180μs   9.11KiB    0.0%     848B
   scaling                 11    488μs    0.1%  44.3μs   2.64MiB    7.9%   246KiB
   LHS+                    11    289μs    0.0%  26.3μs     0.00B    0.0%    0.00B
   load                    11    193μs    0.0%  17.5μs   2.05MiB    6.1%   191KiB
   lhs alloc               11    179μs    0.0%  16.3μs   2.50MiB    7.5%   233KiB
 NNLSinit                   1    336μs    0.0%   336μs    435KiB    1.3%   435KiB```

```
 n_moms = 6
 append!(mim, [[8, 0, 0], [0, 8, 0], [0, 0, 8]])
 threshold = 300
 ────────────────────────────────────────────────────────────────────────────────
                                        Time                    Allocations      
                               ───────────────────────   ────────────────────────
       Tot / % measured:            1.06s /  79.8%           53.9MiB /  61.9%    

 Section               ncalls     time    %tot     avg     alloc    %tot      avg
 ────────────────────────────────────────────────────────────────────────────────
 NNLSmerge: 1st time        1    795ms   93.9%   795ms   25.7MiB   77.1%  25.7MiB
 NNLSmerge                 11   50.8ms    6.0%  4.61ms   7.21MiB   21.6%   671KiB
 NNLSinit                   1    933μs    0.1%   933μs    435KiB    1.3%   435KiB
 ────────────────────────────────────────────────────────────────────────────────```

vs
```

threshold = 300
Ntarget = 100
 ────────────────────────────────────────────────────────────────────────────
                                    Time                    Allocations      
                           ───────────────────────   ────────────────────────
     Tot / % measured:          150ms /   1.1%           23.4MiB /   2.4%    

 Section           ncalls     time    %tot     avg     alloc    %tot      avg
 ────────────────────────────────────────────────────────────────────────────
 collide              500   1.14ms   67.2%  2.28μs    170KiB   29.4%     348B
 merge: 1st time        1    357μs   21.0%   357μs    408KiB   70.6%   408KiB
 merge                 12    201μs   11.8%  16.8μs     0.00B    0.0%    0.00B
 ────────────────────────────────────────────────────────────────────────────```

```Ntarget=200
 ────────────────────────────────────────────────────────────────────────────
                                    Time                    Allocations      
                           ───────────────────────   ────────────────────────
     Tot / % measured:          157ms /   1.6%           23.4MiB /   2.6%    

 Section           ncalls     time    %tot     avg     alloc    %tot      avg
 ────────────────────────────────────────────────────────────────────────────
 collide              500   1.54ms   61.4%  3.09μs    211KiB   34.1%     432B
 merge                 26    586μs   23.4%  22.6μs     0.00B    0.0%    0.00B
 merge: 1st time        1    382μs   15.2%   382μs    408KiB   65.9%   408KiB
 ────────────────────────────────────────────────────────────────────────────```



 ────────────────────────────────────────────────────────────────────
                            Time                    Allocations      
                   ───────────────────────   ────────────────────────
 Tot / % measured:      2.22s /  98.0%            256MiB / 100.0%    

 Section   ncalls     time    %tot     avg     alloc    %tot      avg
 ────────────────────────────────────────────────────────────────────
 ntc            1    2.17s  100.0%   2.17s    256MiB  100.0%   256MiB
   I/O      10.0k    2.15s   99.1%   215μs    256MiB  100.0%  26.2KiB
 ────────────────────────────────────────────────────────────────────
5.0e17
 ────────────────────────────────────────────────────────────────────
                            Time                    Allocations      
                   ───────────────────────   ────────────────────────
 Tot / % measured:      2.37s /  98.1%            256MiB /  99.9%    

 Section   ncalls     time    %tot     avg     alloc    %tot      avg
 ────────────────────────────────────────────────────────────────────
 ntc            1    2.32s  100.0%   2.32s    256MiB  100.0%   256MiB
   I/O      10.0k    2.15s   92.6%   215μs    256MiB  100.0%  26.2KiB
 ────────────────────────────────────────────────────────────────────
5.0e16
 ────────────────────────────────────────────────────────────────────
                            Time                    Allocations      
                   ───────────────────────   ────────────────────────
 Tot / % measured:      3.99s /  98.9%            257MiB /  99.7%    

 Section   ncalls     time    %tot     avg     alloc    %tot      avg
 ────────────────────────────────────────────────────────────────────
 ntc            1    3.94s  100.0%   3.94s    256MiB  100.0%   256MiB
   I/O      10.0k    2.16s   54.8%   216μs    256MiB  100.0%  26.2KiB
 ────────────────────────────────────────────────────────────────────
5.0e15
 ────────────────────────────────────────────────────────────────────
                            Time                    Allocations      
                   ───────────────────────   ────────────────────────
 Tot / % measured:      23.7s /  99.8%            263MiB /  97.4%    

 Section   ncalls     time    %tot     avg     alloc    %tot      avg
 ────────────────────────────────────────────────────────────────────
 ntc            1    23.7s  100.0%   23.7s    256MiB  100.0%   256MiB
   I/O      10.0k    2.77s   11.7%   277μs    256MiB  100.0%  26.2KiB
 ────────────────────────────────────────────────────────────────────


 NCDatasets
  ─────────────────────────────────────────────────────────────────────
                             Time                    Allocations      
                    ───────────────────────   ────────────────────────
  Tot / % measured:      1.52s /   0.7%            505MiB /   1.8%    

 Section    ncalls     time    %tot     avg     alloc    %tot      avg
 ─────────────────────────────────────────────────────────────────────
 merge e         1   6.37ms   64.3%  6.37ms   8.43MiB   94.8%  8.43MiB
 I/O            10   2.28ms   23.0%   228μs    456KiB    5.0%  45.6KiB
 coll n-e       10    989μs   10.0%  98.9μs   16.1KiB    0.2%  1.61KiB
 merge i         1    140μs    1.4%   140μs      400B    0.0%     400B
 merge n         1    126μs    1.3%   126μs     0.00B    0.0%    0.00B
 props          10   7.30μs    0.1%   730ns     0.00B    0.0%    0.00B
 coll n-n       10   1.95μs    0.0%   195ns     0.00B    0.0%    0.00B
 acc e          10   1.15μs    0.0%   115ns     0.00B    0.0%    0.00B
 ─────────────────────────────────────────────────────────────────────
  1.524122 seconds (8.23 M allocations: 505.519 MiB, 8.68% gc time, 97.12% compilation time)
20000
40000
60000
80000
100000
 ─────────────────────────────────────────────────────────────────────
                             Time                    Allocations      
                    ───────────────────────   ────────────────────────
  Tot / % measured:      22.5s /  99.7%           4.41GiB /  99.5%    

 Section    ncalls     time    %tot     avg     alloc    %tot      avg
 ─────────────────────────────────────────────────────────────────────
 I/O          100k    22.3s   99.4%   223μs   4.38GiB   99.8%  45.9KiB
 props        100k   93.5ms    0.4%   935ns     0.00B    0.0%    0.00B
 merge e        29   17.5ms    0.1%   604μs   9.90MiB    0.2%   350KiB
 coll n-e     100k   14.1ms    0.1%   141ns    896KiB    0.0%    9.18B
 coll n-n     100k   10.4ms    0.0%   104ns     0.00B    0.0%    0.00B
 acc e        100k   9.25ms    0.0%  92.5ns     0.00B    0.0%    0.00B
 merge n        14    261μs    0.0%  18.6μs     0.00B    0.0%    0.00B
 merge i         4    151μs    0.0%  37.9μs   1.56KiB    0.0%     400B
 ─────────────────────────────────────────────────────────────────────


  ─────────────────────────────────────────────────────────────────────
                             Time                    Allocations      
                    ───────────────────────   ────────────────────────
  Tot / % measured:      1.89s /  96.7%            110MiB /  78.0%    

 Section    ncalls     time    %tot     avg     alloc    %tot      avg
 ─────────────────────────────────────────────────────────────────────
 I/O          100k    1.70s   93.1%  17.0μs   74.8MiB   87.4%     784B
 props        100k   85.9ms    4.7%   859ns     0.00B    0.0%    0.00B
 merge e        29   14.8ms    0.8%   509μs   9.90MiB   11.6%   350KiB
 coll n-e     100k   10.8ms    0.6%   108ns    896KiB    1.0%    9.18B
 coll n-n     100k   7.65ms    0.4%  76.5ns     0.00B    0.0%    0.00B
 acc e        100k   6.45ms    0.4%  64.5ns     0.00B    0.0%    0.00B
 merge n        14    258μs    0.0%  18.4μs     0.00B    0.0%    0.00B
 merge i         4    157μs    0.0%  39.2μs   1.56KiB    0.0%     400B
 ─────────────────────────────────────────────────────────────────────

  ─────────────────────────────────────────────────────────────────────
                             Time                    Allocations      
                    ───────────────────────   ────────────────────────
  Tot / % measured:      1.90s /  95.9%           86.7MiB /  72.2%    

 Section    ncalls     time    %tot     avg     alloc    %tot      avg
 ─────────────────────────────────────────────────────────────────────
 I/O          100k    1.71s   93.4%  17.1μs   51.9MiB   82.8%     544B
 props        100k   87.5ms    4.8%   875ns     0.00B    0.0%    0.00B
 coll n-e     100k   11.0ms    0.6%   110ns    896KiB    1.4%    9.18B
 coll n-n     100k   8.30ms    0.5%  83.0ns     0.00B    0.0%    0.00B
 merge e        29   6.69ms    0.4%   231μs   9.90MiB   15.8%   350KiB
 acc e        100k   6.54ms    0.4%  65.4ns     0.00B    0.0%    0.00B
 merge n        14    239μs    0.0%  17.1μs     0.00B    0.0%    0.00B
 merge i         4    145μs    0.0%  36.3μs   1.56KiB    0.0%     400B
 ─────────────────────────────────────────────────────────────────────

  ─────────────────────────────────────────────────────────────────────
                             Time                    Allocations      
                    ───────────────────────   ────────────────────────
  Tot / % measured:      1.80s /  96.6%           51.6MiB /  53.4%    

 Section    ncalls     time    %tot     avg     alloc    %tot      avg
 ─────────────────────────────────────────────────────────────────────
 I/O          100k    1.62s   93.3%  16.2μs   16.8MiB   60.9%     176B
 props        100k   86.1ms    5.0%   861ns     0.00B    0.0%    0.00B
 coll n-e     100k   10.4ms    0.6%   104ns    896KiB    3.2%    9.18B
 coll n-n     100k   7.07ms    0.4%  70.7ns     0.00B    0.0%    0.00B
 merge e        29   6.61ms    0.4%   228μs   9.90MiB   35.9%   350KiB
 acc e        100k   6.50ms    0.4%  65.0ns     0.00B    0.0%    0.00B
 merge n        14    240μs    0.0%  17.1μs     0.00B    0.0%    0.00B
 merge i         4    160μs    0.0%  40.1μs   1.56KiB    0.0%     400B
 ─────────────────────────────────────────────────────────────────────