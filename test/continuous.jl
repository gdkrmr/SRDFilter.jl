using Pkg
Pkg.offline(true)
Pkg.activate("srdtest", shared = true)


using GLMakie
using SRDFilter
using CSV
using DataFrames

n = 8
m = 40

df = CSV.read("data/ResultsLai.csv", DataFrame, header = false)

x = 1.0:nrow(df)
y = df[:, 1]

# using JuliaInterpreter
# breakpoint(SRDFilter.interpolate)

@time ycont = SRDFilter.interpolate(x, x, y, n, m);
@time yfir = smoothMS(y, n, m);
extrema(ycont .- yfir)

fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, x, y, label = "y")
lines!(ax, x, yfir, label = "yfir")
lines!(ax, x, ycont, label = "ycont")
ax2 = Axis(fig[2, 1])
lines!(ax2, x, yfir .- ycont, label = "yfir - ycont")
Legend(fig[1, 2], ax)
Legend(fig[2, 2], ax2)

fig2 = Figure()
ax1 = Axis(fig2[1, 1])
ax2 = Axis(fig2[2, 1])
ax3 = Axis(fig2[3, 1])
a1 = SRDFilter.kernelMS(n, m, Float64)
a2 = SRDFilter.a(-m:m, n, m)
lines!(ax1, -m:m, a1, label = "kernelMS")
lines!(ax2, -m:m, a2, label = "a")
lines!(ax3, -m:m, a1 .- a2, label = "kernelMS - a")
extrema(a1 .- a2)




w = [ SRDFilter.wfit(i, n, m) for i in 0:m-1 ]
w2 = filter(x -> x > 0, w)
lw = length(w2)
a, b = SRDFilter.linreg(x[1:m], y[1:m], w)
a, b = SRDFilter.linreg(x[1:lw], y[1:lw], w2)
a, b = SRDFilter.lm_left(x, y, n, m)



ext11 = (x[1] - m: x[1] - 1) .* a .+ b  |> collect;
a, b = SRDFilter.linreg(x[end - m + 1:end], y[end - m + 1:end], reverse(w))
ext12 = (x[end] + 1:x[end] + m) .* a .+ b |> collect;

ext2 = SRDFilter.extendData(y, m, w)
ext21 = ext2[1:m]
ext22 = ext2[end - m + 1:end]
[ext11 ext12 ext21 ext22]


fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, -m + 1:length(x) + m, ext2)
lines!(ax, -m + 1:length(x) + m, [ext11..., y..., ext12...])




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
