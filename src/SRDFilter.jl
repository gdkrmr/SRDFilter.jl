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

include("filter.jl")
include("continuous.jl")


end # module SRDFilter
