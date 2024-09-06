# NOTE: The filter in Schmid et al. (2022) is supposed to be FIR and only works
# on data with equal spacing. Here we expand it to be able to interpolate
# anything. This is an experiment and I have no idea how well this works.
#

# With the high precision numbers from the table:
# # maximum(err_disc) = 8.5
# # mean(err_disc) = 0.25217936249226114
# # maximum(err_cont) = 8.0987135e6
# # mean(err_cont) = 153912.13688329546
#
# with the low precision numbers from the matlab code:
# # maximum(err_disc) = 8.5
# # mean(err_disc) = 0.25217936249226114
# # maximum(err_cont) = 8.0970035e6
# # mean(err_cont) = 153882.95751185901
#
#
# With high precision numbers and the corrected w() function
# # maximum(err_disc) = 8.5
# # mean(err_disc) = 0.25217936249226114
# # maximum(err_cont) = 8.0987135e6
# # mean(err_cont) = 153912.13688329546

"""
returns the values from table 1
a, b, c
"""
function abc(j, n)
    if j == 0 && n == 6
        # 0.00172, 0.02437, 1.64375
        0.001717576, 0.02437382, 1.64375
    elseif j == 0 && n == 8
        # 0.00440, 0.08821, 2.35938
        0.0043993373, 0.088211164, 2.359375
    elseif j == 1 && n == 8
        # 0.00615, 0.02472, 3.63594
        0.006146815, 0.024715371, 3.6359375
    elseif j == 0 && n == 10
        # 0.00118, 0.04219, 2.74688
        0.0011840032, 0.04219344, 2.746875
    elseif j == 1 && n == 10
        # 0.00367, 0.12780, 2.77031
        0.0036718843, 0.12780383, 2.7703125
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
    kappa = a + b / (c - m) ^ 3
    return kappa
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
function a_default(x, alpha, n::Int)
    xsinc = (n + 4) / 2 * x
    return w(x, alpha) * Base.sinc(xsinc)
end


"""
eq. 7

like a but with better passband response

only needed for n > 4
"""
function a_better(x, alpha, n::Int, m::Int)
    nu = isodd(n / 2) ? 1 : 2
    # nu = nu - 2
    xsinc = (n + 4) / 2 * x
    # paper says \sum_j but j \in {0, 1}
    correction0 = kappa(0, n, m) * x * sinpi(nu * x)
    correction1 = kappa(1, n, m) * x * sinpi((2 + nu) * x)
    correction = correction0 + correction1
    # the matlab code uses nu - 2, no idea why, this is the correct formula and
    # numerically closer to the matlab version than nu - 2
    return w(x, alpha) * (Base.sinc(xsinc) + correction)
end


"""
the kernel function, i are supposed to be discrete values

n: degree
m: kernel half-width
"""
function a(i::Number, n, m)
    alpha = 4
    if -m <= i <= m
        if n <= 4
            vals = a_default(x(i, m), alpha, n)
        else
            vals = a_better(x(i, m), alpha, n, m)
        end
    else
        vals = 0.0
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

