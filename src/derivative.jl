
# To calculate the derivative of a kernel interpolation function with respect to
# kernel \( k \), follow these steps:

# 1. **Define the kernel interpolation function**: Suppose \( f(x) = \sum_i w_i
# k(x, x_i) \), where \( w_i \) are weights and \( x_i \) are data points.

# 2. **Differentiate with respect to \( x \)**: Take the derivative of \( f(x) \)
# with respect to \( x \):

#    \[
#         \frac{d}{dx} f(x)
#         = \frac{d}{dx} \left( \sum_i w_i k(x, x_i) \right)
#         = \sum_i w_i \frac{d}{dx} k(x, x_i)
#    \]

# 3. **Differentiate the kernel function**: If \( k(x, x_i) \) is a known function
# (e.g., Gaussian), differentiate it with respect to \( x \):

#    \[
#         \frac{d}{dx} k(x, x_i)
#    \]

# Putting it all together, you’ll have:

# \[
#     \frac{d}{dx} f(x) = \sum_i w_i \frac{d}{dx} k(x, x_i)
# \]

# Compute the derivatives and sum them up as specified.


# Kernel interpolations does the same, see
# https://joshualampert.github.io/KernelInterpolation.jl/stable/tutorial_differentiating_interpolation/#Applying-differential-operators
