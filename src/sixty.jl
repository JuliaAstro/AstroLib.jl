# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

# XXX: trailsign keyword not implemented, not sure it's worth the effort.
# Possible strategy: hold the sign of number in a variable (you can use
# copysign(1, number)), set result[1] equal to dd and multiply the first
# non-zero element (use findfirst for that) by the sign of number.
function sixty(number::T) where {T<:AbstractFloat}
    dd = trunc(abs(number))
    mm = trunc(abs(60 * number) - 60 * dd)
    ss = abs(3600 * number)
    return SVector{3,T}(trunc(number), mm, ss - 3600 * dd - 60 * mm)
end

"""
    sixty(number) -> [deg, min, sec]

### Purpose ###

Converts a decimal number to sexagesimal.

### Explanation ###

The reverse of `ten` function.

### Argument ###

* `number`: decimal number to be converted to sexagesimal.

### Output ###

An array of three `AbstractFloat`, that are the sexagesimal counterpart
(degrees, minutes, seconds) of `number`.

### Example ###

```jldoctest
julia> using AstroLib

julia> sixty(-0.615)
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 -0.0
 36.0
 54.0
```

### Notes ###

Code of this function is based on IDL Astronomy User's Library.
"""
sixty(n::Real) = sixty(float(n))
