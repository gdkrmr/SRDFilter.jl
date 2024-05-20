using Test
using SRDFilter
using CSV
using DataFrames
using Base.Filesystem

degs = [2, 4, 6, 8, 10]
ms = collect(1:10)
Ts = [Float64, Float32]



@testset "example test" begin
    deg = 6 # degree
    m = 7   # kernel halfwidth

    data = Float64[0, 1, -2, 3, -4, 5, -6, 7, -8, 9, 10, 6, 3, 1, 0];
    out = smoothMS(data, deg, m);
    res = [0.15835885, 0.11657466, -0.09224721, 0.03165689, -0.05481473, -0.05436219,
           0.51054827, -0.59067866, -1.21928695, 5.28610520, 10.46161952, 6.82674246,
           2.49236743, 1.04220381, 0.03264660]
    @test out ≈ res
end

function par_to_filename(method, filename, deg, m)
    basename = split(filename, ".")[1]
    return "data/smooth/$(method)_$(basename)_deg_$(deg)_m_$(m).csv"
end

par_to_testresult(method, filename, deg, m, T) =
    CSV.File(par_to_filename(method, filename, deg, m), header = false) |>
    CSV.Tables.matrix .|>
    T

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
                        y = smoothMS(x, deg, m)
                        y1 = smoothMS1(x, deg, m)
                        @test y ≈ testres_ms[:, ci]
                        @test y1 ≈ testres_ms1[:, ci]
                    end
                end
            end
        end
    end
end
