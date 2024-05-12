# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function aitoff(l::T, b::T) where {T<:AbstractFloat}
    l = rem(l, T(360), RoundNearest)
    alpha2 = deg2rad(l/2)
    delta = deg2rad(b)
    r2 = sqrt(T(2))
    f = 2*r2/pi
    sin_alpha2, cos_alpha2 = sincos(alpha2)
    sin_delta, cos_delta = sincos(delta)
    denom = sqrt(1 + cos_delta * cos_alpha2)
    return rad2deg(cos_delta * sin_alpha2 * 2 * r2 / denom / f),
           rad2deg(sin_delta * r2 / denom / f)
end

"""
    aitoff(l, b) -> x, y

### Purpose ###

Convert longitude `l` and latitude `b` to `(x, y)` using an Aitoff projection.

### Explanation ###

This function can be used to create an all-sky map in Galactic coordinates with
an equal-area Aitoff projection.  Output map coordinates are zero longitude
centered.

### Arguments ###

* `l`: longitude, scalar or vector, in degrees.
* `b`: latitude, number of elements as `l`, in degrees.

Coordinates can be given also as a 2-tuple `(l, b)`.

### Output ###

2-tuple `(x, y)`.

* `x`: x coordinate, same number of elements as `l`.  `x` is normalized to be in
  \$[-180, 180]\$.
* `y`: y coordinate, same number of elements as `l`.  `y` is normalized to be in
  \$[-90, 90]\$.

### Example ###

Get \$(x ,y)\$ Aitoff coordinates of Sirius, whose Galactic coordinates are
\$(227.23, -8.890)\$.

```jldoctest
julia> using AstroLib

julia> x, y = aitoff(227.23, -8.890)
(-137.92196683723276, -11.772527357473054)
```

### Notes ###

See AIPS memo No. 46
(ftp://ftp.aoc.nrao.edu/pub/software/aips/TEXT/PUBL/AIPSMEMO46.PS), page 4, for
details of the algorithm.  This version of `aitoff` assumes the projection is
centered at b=0 degrees.

Code of this function is based on IDL Astronomy User's Library.
"""
aitoff(l::Real, b::Real) = aitoff(promote(float(l), float(b))...)

aitoff(lb::Tuple{Real, Real}) = aitoff(lb...)

function aitoff(l::AbstractArray{L}, b::AbstractArray{B}) where {L<:Real,B<:Real}
    if length(l) != length(b)
        throw(ArgumentError("l and b arrays must have the same length"))
    end
    typel = float(L)
    x = similar(l, typel)
    y = similar(b, typel)
    for i in eachindex(l)
        x[i], y[i] = aitoff(l[i], b[i])
    end
    return x, y
end
