"""
    SRDFilter

This module provides functions for smoothing data with a modified sinc kernel
(Schmid, Rath, and , 2022).

Code translated from the GNU Octave code from the supported information to Julia
by Guido Kraemer and Miguel Mahecha.
"""
module SRDFilter

using LinearAlgebra
export smoothMS, smoothMS1

#### this file contains the method as described in Schmid et al. (2022)
include("filter.jl")

#### these are for the continuous interpolation version
include("mskernel.jl")
include("extrapolation.jl")
include("interpolation.jl")


end # module SRDFilter
