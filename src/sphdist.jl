# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function _sphdist(long1::T, lat1::T, long2::T, lat2::T,
                  degrees::Bool) where {T<:AbstractFloat}
    # Convert both points to rectangular coordinates.
    rxy, z1 = polrec(1,   lat1,  degrees=degrees)
    x1, y1  = polrec(rxy, long1, degrees=degrees)
    rxy, z2 = polrec(1,   lat2,  degrees=degrees)
    x2, y2  = polrec(rxy, long2, degrees=degrees)
    # Compute vector dot product for both points.
    cs = x1*x2 + y1*y2 + z1*z2
    # Compute the vector cross product for both points.
    xc = y1*z2 - z1*y2
    yc = z1*x2 - x1*z2
    zc = x1*y2 - y1*x2
    sn = norm((xc, yc, zc))
    # Convert to polar coordinates.
    radius, angle = recpol(cs, sn)
    return degrees ? rad2deg(angle) : angle
end

"""
    sphdist(long1, lat1, long2, lat2[, degrees=true]) -> angular_distance

### Purpose ###

Angular distance between points on a sphere.

### Arguments ###

* `long1`:  longitude of point 1
* `lat1`: latitude of point 1
* `long2`: longitude of point 2
* `lat2`: latitude of point 2
* `degrees` (optional boolean keyword): if `true`, all angles, including the
  output distance, are assumed to be in degrees, otherwise they are all in
  radians.  It defaults to `false`.

### Output ###

Angular distance on a sphere between points 1 and 2, as an `AbstractFloat`.  It
is expressed in radians unless `degrees` keyword is set to `true`.

### Example ###

```jldoctest
julia> using AstroLib

julia> sphdist(120, -43, 175, +22, degrees=true)
82.33037963318735
```

### Notes ###

* `gcirc` function is similar to `sphdist`, but may be more suitable for
  astronomical applications.

Code of this function is based on IDL Astronomy User's Library.
"""
sphdist(long1::Real, lat1::Real,
        long2::Real, lat2::Real; degrees::Bool=false) =
            _sphdist(promote(float(long1), float(lat1),
                             float(long2), float(lat2))... , degrees)
