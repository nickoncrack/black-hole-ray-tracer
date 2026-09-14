Run the following commands to compile and run:

```
gcc main.c -O3 -ffast-math -march=native -fopenmp -lm
./a.out
```

Ray tracing takes ~10 minutes (on my setup) for a 400x400 image using the Schwarzschild metric.
