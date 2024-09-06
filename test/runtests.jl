using Test
using SRDFilter
using CSV
using DataFrames
using Base.Filesystem
using Statistics

degs = [2, 4, 6, 8, 10]
ms = collect(1:10)
Ts = [Float64, Float32]



@testset "example test" begin
    deg = 6 # degree
    m = 7   # kernel halfwidth
    data = Float64[0, 1, -2, 3, -4, 5, -6, 7, -8, 9, 10, 6, 3, 1, 0];
    res = [0.15835885, 0.11657466, -0.09224721, 0.03165689, -0.05481473, -0.05436219,
           0.51054827, -0.59067866, -1.21928695, 5.28610520, 10.46161952, 6.82674246,
           2.49236743, 1.04220381, 0.03264660]
    # This one runs the example from the original code
    @testset "example test discrete" begin
        out = smoothMS(data, deg, m);
        @test out ≈ res
    end

    # This one runs the example from the original code
    @testset "example test continuous" begin
        # TODO: make this work for ranges and Int
        x = collect(1.0:length(data))
        out = SRDFilter.interpolate(x, x, data, deg, m);
        @test out ≈ res

    end
end

@testset "a() vs kernelMS()" begin
    for n in degs
        for m in ms
            a1 = SRDFilter.kernelMS(n, m, Float64)
            a2 = SRDFilter.a(-m:m, n, m)
            @test a1 ≈ a2
        end
    end
end


function par_to_filename(method, filename, deg, m)
    basename = split(filename, ".")[1]
    return "data/smooth/$(method)_$(basename)_deg_$(deg)_m_$(m).csv"
end

par_to_testresult(method, filename, deg, m, T) =
    CSV.File(par_to_filename(method, filename, deg, m), header = false) |>
    CSV.Tables.matrix .|>
    T

err_disc = Float64[]
err_cont = Float64[]

@testset "csv files" begin
    data_files = Filesystem.readdir("data")
    data_files = filter(x -> occursin(r".csv", x), data_files)
    for file in data_files
        data = CSV.File("data/$(file)", header = false) |> CSV.Tables.matrix
        for deg in degs
            for m in ms
                for T in Ts
                    testdata = T.(data)
                    testres_ms = par_to_testresult("MS", file, deg, m, T)
                    testres_ms1 = par_to_testresult("MS1", file, deg, m, T)
                    for ci in 1:size(data, 2)
                        x = testdata[:, ci]
                        y = smoothMS(x, deg = deg, m = m)
                        y1 = smoothMS1(x, deg = deg, m = m)
                        @test y ≈ testres_ms[:, ci]
                        @test y1 ≈ testres_ms1[:, ci]
                        @test eltype(testres_ms) == eltype(y)
                        @test eltype(testres_ms1) == eltype(y1)

                        # TODO: make this work for ranges, Int, and other types
                        t = T.(collect(1:length(x)))
                        ycont = SRDFilter.interpolate(t, t, x, deg, m)
                        @show par_to_filename("MS", file, deg, m) deg m T ci extrema(ycont .- testres_ms[:, ci])
                        # TODO: This is still way to inaccurate
                        @test ycont ≈ testres_ms[:, ci] atol = 1e4
                        # TODO: make this work for arbitrary types
                        #@test eltype(ycont) == eltype(y)

                        push!(err_disc, maximum(abs.(y     .- testres_ms[:, ci])))
                        push!(err_cont, maximum(abs.(ycont .- testres_ms[:, ci])))
                    end
                end
            end
        end
    end
end

@show maximum(err_disc) mean(err_disc) maximum(err_cont) mean(err_cont)
