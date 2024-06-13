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