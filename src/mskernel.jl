# NOTE: The filter in Schmid et al. (2022) is supposed to be FIR and only works
# on data with equal spacing. Here we expand it to be able to interpolate
# anything. This is an experiment and I have no idea how well this works.
#


"""
returns the values from table 1
a, b, c
"""
function abc(j, n)
    if j == 0 && n == 6
        0.00172, 0.02437, 1.64375
    elseif j == 0 && n == 8
        0.00440, 0.08821, 2.35938
    elseif j == 1 && n == 8
        0.00615, 0.02472, 3.63594
    elseif j == 0 && n == 10
        0.00118, 0.04219, 2.74688
    elseif j == 1 && n == 10
        0.00367, 0.12780, 2.77031
    else
        # this makes iteration easier because kappa(j, n, m) becomes 0
        0.0, 0.0, 0.0
    end
end


"""
eq. 8
"""
function kappa(j, n, m)
    a, b, c = abc(j, n)
    a + b / (c - m) ^ 3
end

"""
eq. 4
"""
function w(x, alpha)
    exp1 = exp(-alpha * x * x)
    exp2 = exp(-alpha * (x + 2) ^ 2)
    exp3 = exp(-alpha * (x - 2) ^ 2)
    exp4 = 2 * exp(-alpha)
    exp5 = exp(-9 * alpha)
    return exp1 + exp2 + exp3 - exp4 - exp5
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
function a_default(x, alpha, n::Int)
        xsinc = (n + 4) / 2 * pi * x
        return w(x, alpha) * Base.sinc(xsinc)
end


"""
eq. 7

like a but with better passband response

only needed for n > 4
"""
function a_better(x, alpha, n::Int, m::Int)
    nu = isodd(n / 2) ? 1 : 2
    xsinc = (n + 4) / 2 * x
    correction_f(j) = kappa(n, j, m) * x * sinpi((2 * j + nu) * x)
    correction = sum(correction_f, 0:1)
    return w(x, alpha) * (Base.sinc(xsinc) + correction)
end


"""
the kernel function, i are supposed to be discrete values

n: degree
m: kernel half-width
"""
function a(i::AbstractArray, n, m)
    alpha = 4
    if n <= 4
        vals = a_default.(x.(i, m), alpha, n)
    else
        vals = a_better.(x.(i, m), alpha, n, m)
    end
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

