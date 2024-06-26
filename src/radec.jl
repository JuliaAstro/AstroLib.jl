# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function _radec(ra::T, dec::T, hours::Bool) where {T<:AbstractFloat}
    # Compute right ascension.
    if hours
        ra_hr, ra_min, ra_sec = sixty(mod(ra, 24))
    else
        ra_hr, ra_min, ra_sec = sixty(mod(ra, 360) / 15)
    end
    # Compute declination.
    dec_deg, dec_min, dec_sec = sixty(dec)
    return ra_hr, ra_min, ra_sec, dec_deg, dec_min, dec_sec
end

"""
    radec(ra::Real, dec::Real[, hours=false]) -> ra_hours, ra_minutes, ra_seconds, dec_degrees, dec_minutes, dec_seconds

### Purpose ###

Convert right ascension and declination from decimal to sexagesimal units.

### Explanation ###

The conversion is to sexagesimal hours for right ascension, and sexagesimal
degrees for declination.

### Arguments ###

* `ra`: decimal right ascension, scalar or array.  It is expressed in degrees,
  unless the optional keyword `hours` is set to `true`.
* `dec`: declination in decimal degrees, scalar or array, same number of
  elements as `ra`.
* `hours` (optional boolean keyword): if `false` (the default), `ra` is assumed
  to be given in degrees, otherwise `ra` is assumed to be expressed in hours.

### Output ###

A 6-tuple of `AbstractFloat`:
```
(ra_hours, ra_minutes, ra_seconds, dec_degrees, dec_minutes, dec_seconds)
```

If `ra` and `dec` are arrays, also each element of the output 6-tuple are arrays
of the same dimension.

### Example ###

Position of Sirius in the sky is (ra, dec) = (6.7525, -16.7161), with right
ascension expressed in hours.  Its sexagesimal representation is given by

```jldoctest
julia> using AstroLib

julia> radec(6.7525, -16.7161, hours=true)
(6.0, 45.0, 9.0, -16.0, 42.0, 57.9600000000064)
```
"""
radec(ra::Real, dec::Real; hours::Bool=false) =
    _radec(promote(float(ra), float(dec))..., hours)

function radec(ra::AbstractArray{R}, dec::AbstractArray{D};
               hours::Bool=false) where {R<:Real, D<:Real}
    if length(ra) != length(dec)
        throw(DimensionMismatch("ra and dec arrays should be of the same length"))
    end
    typera = float(R)
    ra_hr   = similar(ra, typera)
    ra_min  = similar(ra, typera)
    ra_sec  = similar(ra, typera)
    dec_deg = similar(ra, typera)
    dec_min = similar(ra, typera)
    dec_sec = similar(ra, typera)
    for i in eachindex(ra)
        ra_hr[i], ra_min[i], ra_sec[i], dec_deg[i], dec_min[i], dec_sec[i] =
            radec(ra[i], dec[i], hours=hours)
    end
    return ra_hr, ra_min, ra_sec, dec_deg, dec_min, dec_sec
end
