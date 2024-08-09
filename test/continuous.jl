using Pkg
Pkg.offline(true)
Pkg.activate("srdtest", shared = true)


using GLMakie
using SRDFilter

n = 4
m = 4

# look at kernel
x = -20:0.01:20
y = SRDFilter.a(x, n, m)
lines(x, y)

# look at weights at the edge
x = 0:0.01:20
y = SRDFilter.wfit.(x, n, m)
lines(x, y)


# look at extended measurements
