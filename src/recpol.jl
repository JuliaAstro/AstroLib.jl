# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function _recpol(x::T, y::T, degrees::Bool) where {T<:AbstractFloat}
    if degrees
        return hypot(x, y), atand(y, x)
    else
        return hypot(x, y), atan(y, x)
    end
end

"""
    recpol(x, y[, degrees=false]) -> radius, angle

### Purpose ###

Convert 2D rectangular coordinates to polar coordinates.

### Explanation ###

This is the partial inverse function of `polrec`.

### Arguments ###

* `x`: the abscissa coordinate of the point.  It may be a scalar or an array.
* `y`: the ordinate coordinate of the point.  It may be a scalar or an array
  of the same lenth as `x`.
* `degrees` (optional boolean keyword): if `true`, the output `angle` is given
  in degrees, otherwise in radians.  It defaults to `false`.

Mandatory arguments may also be passed as the 2-tuple `(x, y)`, so that it is
possible to execute `polrec(recpol(x, y))`.

### Output ###

A 2-tuple `(radius, angle)` with the polar coordinates of the input.  The
coordinate `angle` coordinate lies in the range ``[-π, π]`` if
`degrees=false`, or ``[-180, 180]`` when `degrees=true`.

If `x` and `y` are arrays, `radius` and `angle` are arrays of the same length as
`radius` and `angle`.

### Example ###

Calculate polar coordinates ``(r, φ)`` of point with rectangular
coordinates ``(x, y) = (2.24, -1.87)``.

```jldoctest
julia> using AstroLib

julia> r, phi = recpol(2.24, -1.87)
(2.9179616172938263, -0.6956158538564537)
```

Angle ``φ`` is given in radians.

"""
recpol(x::Real, y::Real; degrees::Bool=false) =
    _recpol(promote(float(x), float(y))..., degrees)

recpol(xy::Tuple{Real, Real}; degrees::Bool=false) =
    recpol(xy..., degrees=degrees)

function recpol(x::AbstractArray{X}, y::AbstractArray{Y};
                degrees::Bool=false) where {X<:Real, Y<:Real}
    if length(x) != length(y)
        throw(DimensionMismatch("x and y arrays must have the same length"))
    end
    typex = float(X)
    r = similar(x, typex)
    a = similar(x, typex)
    for i in eachindex(x)
        r[i], a[i] = recpol(x[i], y[i], degrees=degrees)
    end
    return r, a
end
