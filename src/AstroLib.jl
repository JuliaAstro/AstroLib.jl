# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

__precompile__()

module AstroLib
using Dates, LinearAlgebra, Printf
using StaticArrays
using Interpolations

# Note on function definitions in this package.  Most functions are defined as
# follows:
#
#     _funcname{T<:AbstractFloat}(x::T, y::T) = ...
#     funcname(x::Real, y::Real) = _funcname(promote(float(x), float(y))...)
#
# "_funcname" is the function with the implementation of the full algorithm.
# "_funcname" doesn't have optional argument nor keywords.  Instead, "funcname"
# is the function exposed to the user which ensures all input arguments are
# converted to the same floating point type, in order to improve both
# type-stability and precision of calculations.  "funcname" may have optional
# arguments or keywords.  In addition, the use of AbstractFloat type in
# "_funcname" enables compatibility with Measurements.jl package, which defines
# a new type Measurement{T<:AbstractFloat} <: AbstractFloat.

# Core functions
include("misc.jl")
include("utils.jl")
# New types
include("types.jl")
# Common constants
include("common.jl")

end # module
