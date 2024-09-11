
function expand_left!(x, y, n, m)
    # we use the fact that x is sorted here
    xi_left = findfirst(xi -> xi > x[1] + m, x)::Int - 1
    x_left0 = @view(x[1:xi_left]) .- x[1] .+ 1
    y_left0 = @view(y[1:xi_left])
    #### create x values for extrapolaton
    a_left, b_left = lm_left(x_left0, y_left0, n, m)
    #                        ^^^^^^^^^^^^^^ move x so that x[1] == 1.
    #                        Because slopes can be very large, having x far
    #                        a way from 0 can lead to numerical
    #                        instabilities. The matlab code starts at 1 and
    #                        slight change here can lead to large
    #                        differences in the result. TODO: optimization,
    #                        we only need a couple of values

    # mirror x values to the left of x[1] for extrapolation. The user is
    # responsible to have enough values.
    #
    x_left = reverse(1 .- (@view(x[2:xi_left]) .- x[1]))
    @assert length(x_left) > 0 "need values in range x[1]..(x[1] + m), no values for extending to the left"

    #### extrapolate
    y_left = a_left .* x_left .+ b_left

    # These are the actual values for x_left, not the ones used fo the regression
    x_left_real = x[1] .+ x_left .- 1

    # [x_left_real..., x...]
    prepend!(x, x_left_real)
    # y = [y_left..., y...]
    prepend!(y, y_left)
    return nothing
end

function expand_right!(x, y, n, m)

    #### create x values for extrapolaton
    xi_right = findlast(xi -> xi < x[end] - m, x)::Int + 2
    ### end at -1 so that predictions can start at 0
    x_right0 = @view(x[xi_right:end]) .- x[end] .- 1
    y_right0 = @view(y[xi_right:end])
    a_right, b_right = lm_right(x_right0, y_right0, n, m)

    x_right = reverse((.-x_right0) .- 1)
    #### extrapolate
    y_right = a_right .* x_right .+ b_right

    x_right_real = x_right .+ x[end] .+ 1
    # we change the reference of x here so we don't mutate the input
    # arguments
    # x = [x..., x_right_real...]
    append!(x, x_right_real)
    # y = [y..., y_right...]
    append!(y, y_right)
    return nothing
end

"""
interpolate for a single number
"""
function interpolate(xnew::Number,
                     x::AbstractVector,
                     y::AbstractVector{T},
                     n, m) where {T}
    @assert issorted(x)
    @assert length(x) == length(y)
    @assert x[1] <= xnew <= x[end]

    # we do not want to modify the input arguments!
    x2 = copy(x)
    y2 = copy(y)

    # do we need to extrapolate to the left? We need to extrapolate if the left
    # end of the kernel sticks out of the range of x. Points for extrapolation
    # are [xnew - m, integers < x[1]...]
    if xnew < x2[1] + m
        expand_left!(x2, y2, n, m)
    end

    # do we need to extrapolate to the right? We need to extrapolate if the
    # right side of the kernel sticks out of the range of x. Points for
    # extrapolation are [integers > x[end]..., xnew + m]
    # Same comments as above apply just to the right.
    if xnew > x2[end] - m
        expand_right!(x2, y2, n, m)
    end

    # first index of x inside the kernel
    # TODO: optimize this, x is sorted
    i = findall(xi -> xnew - m <= xi <= xnew + m, x2)::Array{Int}
    ai = a(xnew .- x2[i], n, m)
    ynew = sum(ai .* y2[i])

    if isnan(ynew)
        @error "ynew is NaN" xnew ynew ai i x2[i] n m
        error("ynew is NaN")
    end

    return ynew
end

"""
interpolate for a vector of new values
"""
function interpolate(xnew::AbstractVector,
                     x::AbstractVector,
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
        ynew[i] = interpolate(xnew[i], x, y, n, m)
    end

    return ynew
end
