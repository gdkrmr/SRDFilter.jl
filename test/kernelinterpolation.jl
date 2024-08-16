#### - this doesn't work with KernelInterpolation.jl as is because there is a
####   normalization factor A (eq. 3 and 7) that is a / sum(a).
####
#### - Also, these kernels are indefinite and may just not work as is. Work
####   around could be regularization.
####
#### -

using Pkg
Pkg.offline(true)
Pkg.activate("srdtest", shared = true)


using GLMakie
using SRDFilter
using KernelInterpolation
using CSV
using DataFrames

n = 10
m = 60

# look at kernel
x = -20:0.01:20
y = SRDFilter.a(x, n, m)
lines(x, y)

# look at weights at the edge
x = 0:0.01:20
y = SRDFilter.wfit.(x, n, m)
lines(x, y)


# look at extended measurements

struct ModifiedSincKernel{Dim} <: KernelInterpolation.RadialSymmetricKernel{Dim} end

KernelInterpolation.phi(::ModifiedSincKernel, r) = SRDFilter.a(r, n, m)
KernelInterpolation.order(::ModifiedSincKernel) = 0

df = CSV.read("data/ResultsLai.csv", DataFrame, header = false)

x = 1.0:nrow(df)
y = df[:, 1]

nodes = NodeSet(x)
k = ModifiedSincKernel{1}()
itp = interpolate(nodes, y, k)
ycont = itp.(x)

yfir = smoothMS(y, n, m)

fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, x, y, label = "y")
lines!(ax, x, ycont, label = "ycont")
lines!(ax, x, yfir, label = "yfir")
Legend(fig[1, 2], ax)

ax2 = Axis(fig[2, 1])
lines!(ax2, x, y, label = "y")

ax3 = Axis(fig[3, 1])
lines!(ax3, x, ycont, label = "ycont")

ax4 = Axis(fig[4, 1])
lines!(ax4, x, yfir, label = "yfir")

ax5 = Axis(fig[5, 1])
lines!(ax5, x, y .- ycont)

ax6 = Axis(fig[6, 1])
lines!(ax6, x, SRDFilter.x.(x, m))

ax7 = Axis(fig[7, 1])
lines!(ax7, x, SRDFilter.w.(SRDFilter.x.(x, m), 4))
vlines!(ax7, m)

ax8 = Axis(fig[8, 1])
lines!(ax8, x, SRDFilter.a.(x, n, m))
vlines!(ax8, m)
