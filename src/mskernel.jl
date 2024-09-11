# NOTE: The filter in Schmid et al. (2022) is supposed to be FIR and only works
# on data with equal spacing. Here we expand it to be able to interpolate
# anything. This is an experiment and I have no idea how well this works.
#

"""
returns the values from table 1
a, b, c
"""
function abc(j, n, ::Type{T}) where T
    if j == 0 && n == 6
        T(0.001717576), T(0.02437382), T(1.64375)
    elseif j == 0 && n == 8
        T(0.0043993373), T(0.088211164), T(2.359375)
    elseif j == 1 && n == 8
        T(0.006146815), T(0.024715371), T(3.6359375)
    elseif j == 0 && n == 10
        T(0.0011840032), T(0.04219344), T(2.746875)
    elseif j == 1 && n == 10
        T(0.0036718843), T(0.12780383), T(2.7703125)
    else
        # this makes iteration easier because kappa(j, n, m) becomes 0
        T(0), T(0), T(0)
    end
end


"""
eq. 8
"""
function kappa(j, n, m, ::Type{T}) where T
    a, b, c = abc(j, n, T)
    kappa = a + b / (c - m) ^ 3
    return kappa
end

"""
eq. 4
"""
function w(x::T, alpha) where T
    alphaT = T(alpha)
    exp1 = exp(-alphaT * x * x)
    exp2 = exp(-alphaT * (x + 2) ^ 2)
    exp3 = exp(-alphaT * (x - 2) ^ 2)
    exp4 = 2 * exp(-alphaT)
    exp5 = exp(-9 * alphaT)
    return exp1 + exp2 + exp3 - exp4 - exp5
    #                                ^ corrected this from + which is in the octave code.
end

"""
eq. 5
"""
x(i, m) = i / (m + 1)


"""
eq. 3

kernel function for sinc

w: window function
n: degree

used for n <= 4
"""
function a_default(x::T, alpha, n::Int) where T
    xsinc = T(n + 4) / 2 * x
    return w(x, alpha) * Base.sinc(xsinc)
end


"""
eq. 7

like a but with better passband response

only needed for n > 4
"""
function a_better(x::T, alpha, n::Int, m::Int) where T
    nu = isodd(n / 2) ? 1 : 2
    xsinc = T(n + 4) / 2 * x
    # paper says \sum_j but j \in {0, 1}
    correction0 = kappa(0, n, m, T) * x * sinpi(nu * x)
    correction1 = kappa(1, n, m, T) * x * sinpi((2 + nu) * x)
    correction = correction0 + correction1
    # the matlab code uses nu - 2, this is because they iterate for j in 1:2 and
    # not 0:1
    return w(x, alpha) * (Base.sinc(xsinc) + correction)
end


"""
the kernel function, i are supposed to be discrete values

n: degree
m: kernel half-width
"""
function a(i::T, n, m) where T <: Number
    alpha = 4
    if -m <= i <= m
        if n <= 4
            vals = a_default(x(i, m), alpha, n)
        else
            vals = a_better(x(i, m), alpha, n, m)
        end
    else
        vals = T(0.0)
    end
    return vals
end

function a(i::AbstractArray, n, m)
    vals = a.(i, n, m)
    Ainv = sum(vals) # eq. 6
    return vals ./ Ainv
end

########## b and m_ms translate the m used in SG-Filter to an equivalent m_ms
########## for the sinc filter
"""
Schmid et al. (2022), eq. 14

n: degree
"""
function b(m, n)
    a = 6.352 * (m + 0.5) / (n + 1.379) - (0.512 + 0.316 * n) / (m + 0.5)
    return 1 / a
end


"""
Schmid et a. (2022) eq. 1

modified sinc kernel with linear extrapolation

m value for b

kernel half width m should be rounded to nearest integer
"""
m_ms(n, b) = (0.745 + 0.249 * n) / b - 1
#########################################

