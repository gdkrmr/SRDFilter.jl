using Pkg
# Pkg.offline(true)
Pkg.activate("srdtest", shared = true)


using LinearAlgebra
using GLMakie
using SRDFilter
using CSV
using DataFrames

n = 8
m = 40

df = CSV.read("data/ResultsLai.csv", DataFrame, header = false)

x = 1.0:nrow(df)
xnew = 1.0:0.01:nrow(df)
y = df[:, 1]

# using JuliaInterpreter
# breakpoint(SRDFilter.interpolate)

@time ycont = SRDFilter.interpolate(x, x, y, n, m);
@time ycontnew = SRDFilter.interpolate(xnew, x, y, n, m);
@time yfir = smoothMS(y, n, m);
extrema(ycont .- yfir)

fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, x, y, label = "y")
lines!(ax, x, yfir, label = "yfir")
lines!(ax, x, ycont, label = "ycont")
lines!(ax, xnew, ycontnew, label = "ycontnew")
ax2 = Axis(fig[2, 1])
lines!(ax2, x, yfir .- ycont, label = "yfir - ycont")
Legend(fig[1, 2], ax)
Legend(fig[2, 2], ax2)


n = 10
m = 10
T = Float32
fig2 = Figure()
ax1 = Axis(fig2[1, 1])
ax2 = Axis(fig2[2, 1])
ax3 = Axis(fig2[3, 1])
a1 = SRDFilter.kernelMS(n, m, Float64)
a2 = SRDFilter.a(T.(collect(-m:m)), n, m)
lines!(ax1, -m:m, a1, label = "kernelMS")
lines!(ax2, -m:m, a2, label = "a")
lines!(ax3, -m:m, a1 .- a2, label = "kernelMS - a")
extrema(a1 .- a2)



n = 10
m = 10
T = Float32
ci = 2
file = "../test/data/smooth/MS_ResultsLai_deg_10_m_10.csv"
data_smooth = CSV.File(file, header = false) |>
    CSV.Tables.matrix |>
    x -> x[:, ci] .|>
    T
t = eltype(data_smooth).(collect(1:length(data_smooth)))
y = SRDFilter.smoothMS(data_smooth, n, m)
ycont = SRDFilter.interpolate(t, t, data_smooth, n, m)
extrema(y .- ycont)
lines(y .- ycont)

varX2 = 0.06217410988057037
offset = -16210.454835929333
slope = 39661.63646306257
fitX = 1.0:1.0:2.0
fitY = [-793446.285239884, -771325.532923279]
fitWeights = [1.0, 0.06217410988057054]
varX2 = 0.06217410988057037
offset = -815567.0375564876
slope = 22120.75231660364


a_left = 39661.63646306236
b_left = -16210.454835929097
a_right = -22120.75231086372
b_right = 1.91152307945378e7

begin
T = Float32
xlm = T[900.0, 899.0]
w = T[1.0, 0.06217411]
ylm = T[-793446.285239884, -771325.532923279]
sum(xlm .* xlm .* w) * sum(w) - sum(xlm .* w) * sum(xlm .* w)
([xlm [1, 1]]' * Diagonal(w) * [xlm [1, 1]]) \ ([xlm [1, 1]]' * Diagonal(w) * ylm)
end;
SRDFilter.linreg(xlm, ylm, w)
SRDFilter.fitWeighted(xlm, ylm, w)


file = "data/smooth/MS_ResultsLai_deg_2_m_1.csv"






ai .- kernel
yi .- extData[end - length(yi) + 1:end]


degs = [2, 4, 6, 8, 10]
ms = collect(1:10)
for n in degs
    for m in ms
        a1 = SRDFilter.kernelMS(n, m, Float64)
        a2 = SRDFilter.a(-m:m, n, m)
        @show n m extrema(a1 .- a2)
    end
end

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
