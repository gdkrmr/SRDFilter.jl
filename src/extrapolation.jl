
##################### code for extrapolation
"""
eq. 18

beta for the boundary extrapolation
"""
beta(n) = 0.7 + 0.14 * exp(-0.6 * (n - 4))

"""

weights for the linear extrapolation of the data for the MS filter

eq. 17
"""
function wfit(i, n, m)
    num = i * (n + 3)
    denom = 4 * beta(n) * (m + 1)
    x = num / denom
    if x < 0.5
        return cospi(x) ^ 2
    else
        return 0.0
    end
end

"""
weightedlinear regression
"""
function linreg(x, y, w)
    # TODO: look up if this is correct
    @assert length(x) == length(y) == length(w)
    l = length(x)
    wx = w .* x
    wx2 = [ wx zeros(l) ] # do the intercept
    res = (wx2' * wx2) \ (wx2' * y)
    a, b = res
    return a, b
end

"""

extrapolate x, y to the right with steplength dx and weights wfit(i, n, m)

- x and y have to be the same length
- x has to be sorted
"""
function lm_right(x::AbstractArray{T},
                  y::AbstractArray{T},
                  n::Int, m::Int) where {T}
    @assert length(x) == length(y)
    @assert issorted(x)

    ###### extend right side
    i = length(x) - 1
    w = T[1] # this is the point at the border
    while true
        # TODO: do we need to test for i < 1 here?
        if i < 1
            break
        end
        xi = x[i] # x must be sorted
        dist_border = x[end] - xi
        wi = wfit(dist_border, n, m)
        if wi > 0 && i >= x[1]
            #^^^ make sure that we are inside the data range
            #  ^^^^ assume that x is sorted
            #       ^^^^^^ wi will be 0 at a certain distance from the border
            push!(w, wi)
        else
            # we end the loop either when wi has become 0 or when we run out of data
            break
        end
        i -= 1
    end

    xlm = @view x[end:-1:end - i]
    ylm = @view x[end:-1:end - 1]

   # TODO: to the linear regression and return a, b
   a, b = linreg(xlm, ylm, w)
   return a, b
end

"""

extrapolate x, y to the left with steplength dx and weights wfit(i, n, m)

- x and y have to be the same length
- x has to be sorted
"""
function lm_left(x::AbstractArray{T},
                 y::AbstractArray{T},
                 n::Int, m::Int) where {T}
    @assert length(x) == length(y)
    @assert issorted(x)

    ###### extend right side
    i = 2
    w = T[1] # this is the point at the border
    while true
        # TODO: do we need to test for i > length(x) here?
        if i > length(x)
            break
        end
        xi = x[i] # x must be sorted
        dist_border = xi - x[1]
        wi = wfit(dist_border, n, m)
        if wi > 0 && i >= x[1]
            #^^^ make sure that we are inside the data range
            #  ^^^^ assume that x is sorted
            #       ^^^^^^ wi will be 0 at a certain distance from the border
            push!(w, wi)
        else
            # we end the loop either when wi has become 0 or when we run out of data
            break
        end
        i += 1
    end

    xlm = @view x[1:lw]
    ylm = @view x[1:lw]

    a, b = linreg(xlm, ylm, w)

    return a, b
end
