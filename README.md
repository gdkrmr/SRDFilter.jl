# Modified Sinc Kernel Smoothing

Details see [Schmid et al. (2022)](https://pubs.acs.org/doi/10.1021/acsmeasuresciau.1c00054). Code translated from their Matlab
implementation and tested against it.

# Usage

```julia
using SRDFilter

deg = 4
m = 4

data = Float64[0, 1, -2, 3, -4, 5, -6, 7, -8, 9, 10, 6, 3, 1, 0];
out = smoothMS(data, deg, m);
out = smoothMS1(data, deg, m);
```
