

"""
interpolate for a single number
"""
function interpolate(xnew::Number,
                     x, y::AbstractVector{T},
                     n, m) where {T}
    @assert issorted(x)
    @assert length(x) == length(y)
    @assert x[1] <= xnew <= x[end]

    # do we need to extrapolate to the left? We need to extrapolate if the left
    # end of the kernel sticks out of the range of x. Points for extrapolation
    # are [xnew - m, integers < x[1]...]
    if xnew < x[1] + m
        #### create x values for extrapolaton
        a_left, b_left = lm_left(x, y, n, m)
        x_left = collect(ceil(xnew - m):floor(x[1]))
        if x_left[1] != xnew - m
            pushfirst!(x_left, xnew - m)
        end
        if x_left[end] == x[1]
            pop!(x_left)
        end
        #### extrapolate
        y_left = a_left .* x_left .+ b_left
        # we change the reference of x here so we don't mutate the input
        # arguments
        x = [x_left..., x...]
        y = [y_left..., y...]
    end

    # do we need to extrapolate to the right? We need to extrapolate if the
    # right side of the kernel sticks out of the range of x. Points for
    # extrapolation are [integers > x[end]..., xnew + m]
    if xnew > x[end] - m
        #### create x values for extrapolaton
        a_right, b_right = lm_right(x, y, n, m)
        x_right = collect(ceil(x[end]):floor(xnew + m))
        if x_right[end] != xnew + m
            push!(x_right, xnew + m)
        end
        if x_right[1] == x[end]
            popfirst!(x_right)
        end
        #### extrapolate
        y_right = a_right .* x_right .+ b_right
        # we change the reference of x here so we don't mutate the input
        # arguments
        x = [x..., x_right...]
        y = [y..., y_right...]
    end

    # first index of x inside the kernel
    # TODO: optimize this, x is sorted
    i = findall(xi -> xnew - m <= xi <= xnew + m, x)
    ai = a(xnew .- x[i], n, m)
    ynew = sum(ai .* y[i])

    if isnan(ynew)
        @show xnew ynew ai i x[i] n m
        error("ynew is NaN")
    end

    return ynew
end

"""
interpolate for a vector of new values
"""
function interpolate(xnew::AbstractVector,
                     x,
                     y::AbstractVector{T},
                     n, m) where T
    @assert issorted(x)
    @assert issorted(xnew)
    @assert length(x) == length(y)
    l = length(xnew)

    # TODO: only do this if required
    # a_right, b_right = lm_right(x, y, n, m)
    # a_left, b_left = lm_left(x, y, n, m)

    ynew = Array{T}(undef, l)

    for i in 1:l
        xnewi = xnew[i]
        ynew[i] = interpolate(xnewi, x, y, n, m)
    end

    return ynew
end
