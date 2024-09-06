# This function smooths the data all the way to the boundaries
# by convolution with a MS kernel. Near-boundary values are handled
# by (weighted linear) extrapolation before convolution.
#
# 'data' should be a row vector (1xN).
#
# 'deg' degree, determines the sharpness of the cutoff in the frequency domain,
# similar to the degree of Savitzky-Golay filters.
#
# 'm' is the halfwidth of the kernel; higher values lead to stronger smoothing.
#
#
# This code was taken and translated from the matlab version of the original
# publication (Schmid, Rath, and Diebold, 2022). The code ist tested numerically
# against that version (see test/ directory).

"""
    smoothMS(data::AbstractVector{T}, deg::Int, m::Int) where T

Smoothes a vector `data` with a modified sinc kernel. The `deg` parameter
specifies the degree of the polinomial and `m` the halfwidth of the kernel.

"""
function smoothMS(data::AbstractVector{T}, deg::Int, m::Int) where T
    kernel = kernelMS(deg, m, T)
    fitWeights = edgeWeights(deg, m, T)
    extData = extendData(data, m, fitWeights)
    smoothedExtData = conv(extData, kernel, "same")
    return smoothedExtData[m + 1:end - m]
end

smoothMS(data; deg=4, m=6) = smoothMS(data, deg, m)

# The same with the shorter MS1 kernel
"""
    smoothMS1(data::AbstractVector{T}, deg::Int, m::Int) where T

Smoothes a vector `data` with a modified sinc kernel. The `deg` parameter
specifies the degree of the polinomial and `m` the halfwidth of the kernel.

"""
function smoothMS1(data::AbstractVector{T}, deg::Int, m::Int) where {T}
    kernel = kernelMS1(deg, m, T)
    fitWeights = edgeWeights1(deg, m, T)
    extData = extendData(data, m, fitWeights)
    smoothedExtData = conv(extData, kernel, "same")
    return smoothedExtData[m+1:end-m]
end
smoothMS1(data; deg=4, m=6) = smoothMS1(data, deg, m)

# TODO: Make these static?
# Correction coeffficients for a flat passband of the MS kernel
function corrCoeffsMS(deg::Int, T::Type)
    if deg == 2
        return Matrix{T}(undef, 0, 0)
    elseif deg == 4
        return Matrix{T}(undef, 0, 0)
    elseif deg == 6
        return T[0.001717576 0.02437382 1.64375;]
    elseif deg == 8
        return T[0.0043993373 0.088211164 2.359375;
                 0.006146815 0.024715371 3.6359375]
    elseif deg == 10
        return T[0.0011840032 0.04219344 2.746875;
                 0.0036718843 0.12780383 2.7703125]
    else
        error("Invalid deg")
    end
end

# Correction coeffficients for a flat passband of the MS1 kernel
function corrCoeffsMS1(deg::Int, T::Type)
    if deg == 2
        return Matrix{T}(undef, 0, 0)
    elseif deg == 4
        return T[0.021944195 0.050284006 0.765625;]
    elseif deg == 6
        return T[0.001897730 0.00847681 1.2625;
                 0.023064667 0.13047926 1.2265625]
    elseif deg == 8
        return T[0.006590300 0.05792946 1.915625;
                 0.002323448 0.01029885 2.2726562;
                 0.021046653 0.16646601 1.98125]
    elseif deg == 10
        return T[9.749618e-4 0.00207429 3.74375;
                 0.008975366 0.09902466 2.707812;
                 0.002419541 0.01006486 3.296875;
                 0.019185117 0.18953617 2.784961]
    else
        error("Invalid deg")
    end
end

# Calculates the MS convolution kernel
function kernelMS(deg::Int, m::Int, T::Type)
    coeffs = corrCoeffsMS(deg, T)
    kappa = Vector{T}(undef, size(coeffs, 1))
    for (i, row) in enumerate(eachrow(coeffs))
        kappa[i] = row[1] + row[2] / (row[3] - m)^3
    end
    nuMinus2 = (rem(T(deg) / 2, 2) == 1) ? T(-1) : T(0)
    kernel = zeros(T, m * 2 + 1)
    kernel[m+1] = windowMS(T(0), 4) # center element
    for i in 1:m
        x = T(i) / (m + 1)
        w = windowMS(x, 4)
        a = sin((T(deg) / 2 + 2) * pi * x) / ((T(deg) / 2 + 2) * pi * x)
        for (j, k) in enumerate(kappa)
            a += k * x * sin((2 * j + nuMinus2) * pi * x)
        end
        a *= w
        kernel[m + 1 - i] = a
        kernel[m + 1 + i] = a
    end
    norm = sum(kernel)
    return kernel ./ norm
end

# Calculates the MS1 convolution kernel
function kernelMS1(deg::Int, m::Int, ::Type{T}) where T
    coeffs = corrCoeffsMS1(deg, T)
    kappa = Vector{T}(undef, size(coeffs, 1))
    for (i, row) in enumerate(eachrow(coeffs))
        kappa[i] = row[1] + row[2] / (row[3] - m)^3
    end
    kernel = zeros(T, m * 2 + 1)
    kernel[m + 1] = windowMS(0, 2) # center element
    for i in 1:m
        x = T(i) / (m + 1)
        w = windowMS(x, 2)
        a = sin((T(deg) / 2 + 1) * pi * x) / ((T(deg) / 2 + 1) * pi * x)
        for (j, k) in enumerate(kappa)
            a += k * x * sin(x * j * pi)
        end
        a *= w
        kernel[m + 1 - i] = a
        kernel[m + 1 + i] = a
    end
    norm = sum(kernel)
    return kernel ./ norm
end

"""

This function behaves the same way as the conv function in
Matlab/Octave, here only the "same" option is implemented.

"""
function conv(x::Vector{T}, y::Vector{T}, shape = "same") where T
    nx = length(x)
    ny = length(y)
    nz = nx + ny - 1
    z = zeros(T, nz)
    for i in 1:nz
        for j in max(1, i - ny + 1):min(i, nx)
            z[i] += x[j] * y[i - j + 1]
        end
    end
    # Return only the "same" part of the result
    istart = ceil(Int, (ny - 1) / 2) + 1
    iend = istart + nx - 1
    return z[istart:iend]
end

# Weighted linear fit of the data.
# All inputs must be row vectors of equal length.
function  fitWeighted(xData::AbstractVector{T}, yData::Vector{T}, weights) where T
    n = length(xData)
    sumWeights = sum(weights)
    sumX  = sum(i -> xData[i] * weights[i], 1:n)
    sumY  = sum(i -> yData[i] * weights[i], 1:n)
    sumX2 = sum(i -> xData[i] * xData[i] * weights[i] , 1:n)
    sumXY = sum(i -> xData[i] * yData[i] * weights[i], 1:n)
    varX2 = sumX2 * sumWeights - sumX * sumX
    if varX2 == 0
        slope = T(0)
    else
        slope = (sumXY * sumWeights - sumX * sumY) / varX2
    end
    offset = (sumY - slope * sumX) / sumWeights;
    return offset, slope
end

function extendData(data::Vector{T}, m, fitWeights::Vector{T}) where T
    datLength = length(data)
    extData = zeros(T, datLength + 2*m)
    fitLength = length(fitWeights)
    fitY = data[1:fitLength]
    fitX = T(1):T(fitLength)
    offset, slope = fitWeighted(fitX, fitY, fitWeights)
    #fitBasis = [ones(1, fitLength) 1:fitLength];
    #[params] = LinearRegression (fitBasis', fitY', fitWeights')
    extData[1:m] = offset .+ (-m + 1:0) * slope
    extData[m + 1:datLength + m] = data
    fitY = reverse(data[(datLength - fitLength + 1):datLength])
    offset, slope = fitWeighted(fitX, fitY, fitWeights)
    #[params] = LinearRegression (fitBasis', fitY', fitWeights')
    #extData(datLength+m+1 : datLength+2*m) = [ones(1, m) 0:-1:-m+1]' * params;
    extData[datLength+m+1 : datLength+2*m] = offset .+ (0:-1:-m+1) * slope
    return extData
end

# Gaussian-like window function for the MS and MS1 kernels.
# The function reaches 0 at x=+1 and x=-1 (where the kernel ends);
# at these points also the 1st derivative is very close to zero.
function windowMS(x::T, alpha) where T
    alphaT = T(alpha)
    return exp(-alphaT * x * x) + exp(-alphaT * (x + 2) * (x + 2)) + exp(-alphaT * (x - 2) * (x - 2)) - (2 * exp(-alphaT) - exp(-9 * alphaT))
end

# Hann-square weights for linear fit at the edges, for MS smoothing.
function edgeWeights(deg, m, ::Type{T}) where T
    beta = T(0.70) + T(0.14) * exp(T(-0.6) * (deg - 4))
    fitLengthD = ((m + 1) * beta)/(T(1.5) + T(0.5) * deg)
    fitLength = floor(Int, fitLengthD)
    w = zeros(T, fitLength + 1)
    for i in 1:fitLength + 1
        cosine = cos(T(0.5) * pi * (i - 1) / fitLengthD)
        w[i] = cosine^2
    end
    return w
end

# Hann-square weights for linear fit at the edges, for MS1 smoothing.
function edgeWeights1(deg, m, ::Type{T}) where T
    beta = T(0.65) + T(0.35) * exp(T(-0.55) * (deg - 4))
    fitLengthD = ((m + 1) * beta) / (1 + T(0.5) * deg)
    fitLength = floor(Int, fitLengthD)
    w = zeros(T, fitLength + 1)
    for i in 1:fitLength + 1
        cosine = cos(T(0.5) * pi * (i - 1) / fitLengthD)
        w[i] = cosine^2
    end
    return w
end
