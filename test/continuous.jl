using Pkg
# Pkg.offline(true)
Pkg.activate("srdtest", shared = true)


using LinearAlgebra
using GLMakie
using SRDFilter
using CSV
using DataFrames
using SRDFilter


function par_to_filename(method, filename, deg, m)
    basename = split(filename, ".")[1]
    return "../test/data/smooth/$(method)_$(basename)_deg_$(deg)_m_$(m).csv"
end

par_to_testresult(method, filename, deg, m, T) =
    CSV.File(par_to_filename(method, filename, deg, m), header = false) |>
    CSV.Tables.matrix .|>
    T

i = 1
ci = 1
T = Float32
m = 1
deg = 2
file = "ResultsLai.csv"
data = CSV.File("../test/data/$(file)", header = false) |> CSV.Tables.matrix |> x -> x[:, ci] .|> T
res = par_to_testresult("MS", file, deg, m, T)[:, ci]
smoothfil = SRDFilter.smoothMS(data, deg, m)
x = collect(T(1):length(data)) .+ 1_000_000
smothcont = SRDFilter.interpolate(x, x, data, deg, m)


@code_warntype SRDFilter.expand_right!(x, data, deg, m)
@code_warntype SRDFilter.lm_left(x, data, deg, m)
w = T[1, 2, 3, 4]
@code_warntype SRDFilter.linreg(x, data, w)

@code_warntype SRDFilter.interpolate(x, x, data, deg, m)
@code_warntype SRDFilter.interpolate(x[1], x, data, deg, m)

using JET
@report_opt SRDFilter.interpolate(x, x, data, deg, m)
@report_opt SRDFilter.linreg(x, data, w)

using Cthulhu
@descend SRDFilter.interpolate(x, x, data, deg, m)

a, b = [1, 2]

function f(x::AbstractVector{T}, y, w) where T
    x2 = [x ones(T, length(x))]
    x3 = (x2' * Diagonal(w) * x2)
    y2 = (x' * Diagonal(w) * y)
     x3 \ y2
end

x = T[1, 2, 3]
y = T[1, 2, 3]
w = T[1, 2, 3]
@report_opt f(x, y, w)


f(::Type{T}) where T = string(T)
f(Float32)

@code_warntype SRDFilter.a(1f0, 2, 2)
@code_warntype SRDFilter.a([1f0], 2, 2)
@code_warntype SRDFilter.a_default(1f0, 4, 2)
@code_warntype SRDFilter.a_better(1f0, 4, 2, 2)
@code_warntype SRDFilter.w(1f0, 4)


fig = Figure()
ax0 = Axis(fig[1, 1])
lines!(ax0, x, data, label = "data")
lines!(ax0, x, smoothfil, label = "smoothfil")
lines!(ax0, x, smothcont, label = "smothcont")
lines!(ax0, x, res, label = "res")
Legend(fig[1, 2], ax0)
ax1 = Axis(fig[2, 1])
lines!(ax1, x, smoothfil .- res, label = "smoothfil - smothcont")
lines!(ax1, x, smothcont .- res, label = "smothcont - res")
Legend(fig[2, 2], ax1)









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

data_rough = "../test/data/ResultsLai.csv" |>
    x -> CSV.File(x, header = false) |>
    CSV.Tables.matrix |>
    x -> x[:, ci] .|>
    T


SRDFilter.interpolate(T(1), collect(T(1):length(data_rough)), data_rough, n, m)

x_left = Float32[-9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0]
y_left = [-277908.95248296217, -247151.44336747847, -216393.93425199474, -185636.42513651104, -154878.9160210273, -124121.4069055436, -93363.89779005988, -62606.388674576156, -31848.879559092442, -1091.3704436087235]
23451.18275928195


Evaluated:
[0.15835884531613056, -0.46217518548785597, -1.008098621556303, -1.5204907403192778, -1.9821019170810614, -2.512756040996683, -2.957128829118962, -0.5906786605713916, -1.2192869459451743, 5.286105202110527, 10.461619519603236, 6.826742464105779, 2.492367430378482, 1.0422038091960275, 0.032646599192894275]
[0.15835885,           0.11657466,          -0.09224721,         0.03165689,         -0.05481473,         -0.05436219,         0.51054827,        -0.59067866,         -1.21928695,         5.2861052,         10.46161952, 6.82674246, 2.49236743, 1.04220381, 0.0326466]




deg = 6 # degree
m = 7   # kernel halfwidth
data = Float64[0, 1, -2, 3, -4, 5, -6, 7, -8, 9, 10, 6, 3, 1, 0];
x = collect(1.0:length(data))
res = [0.15835885, 0.11657466, -0.09224721, 0.03165689, -0.05481473, -0.05436219,
       0.51054827, -0.59067866, -1.21928695, 5.28610520, 10.46161952, 6.82674246,
       2.49236743, 1.04220381, 0.03264660]

SRDFilter.interpolate(x, x, data, deg, m)
SRDFilter.smoothMS(data, deg, m)


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


(xlm, ylm, w) = ([1.0], [1.19787587285916e6], [1.0])

SRDFilter.linreg(xlm, ylm, w)
SRDFilter.fitWeighted(xlm, ylm, w)

[xlm [1]]' * Diagonal(w) * [xlm [1]] |> rank

sum(xlm .* xlm .* w) * sum(w) - sum(xlm .* w) * sum(xlm .* w)




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
