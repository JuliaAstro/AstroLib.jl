# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function geo2eci(lat::T, long::T, alt::T, jd::T) where {T<:AbstractFloat}
    Re    = planets["earth"].eqradius / 1000
    sin_lat, cos_lat = sincos(deg2rad(lat))
    long  = deg2rad(long)
    gst   = ct2lst(zero(T), jd)
    sid_angle = gst*pi/12 # Sidereal angle.
    theta = long + sid_angle # Azimuth
    altRe = alt + Re
    r     = altRe * cos_lat
    sin_theta, cos_theta = sincos(theta)
    return r * cos_theta, r * sin_theta, altRe * sin_lat
end

"""
    geo2eci(latitude, longitude, altitude, jd) -> x, y, z

### Purpose ###

Convert geographic spherical coordinates to Earth-centered inertial coordinates.

### Explanation ###

Converts from geographic spherical coordinates (latitude, longitude, altitude)
to ECI (Earth-Centered Inertial) (x, y, z) rectangular coordinates.  Julian days
is also needed.

Geographic coordinates assume the Earth is a perfect sphere, with radius equal
to its equatorial radius.  ECI coordinates are in km from Earth center at epoch
TOD (True of Date).

### Arguments ###

* `latitude`: geographic latitude, in degrees.
* `longitude`: geographic longitude, in degrees.
* `altitude`: geographic altitude, in kilometers.
* `jd`: Julian days.

The three coordinates can be passed as a 3-tuple `(latitude, longitude,
altitude)`.  In addition, `latitude`, `longitude`, `altitude`, and `jd` can be
given as arrays of the same length.

### Output ###

The 3-tuple of ECI (x, y, z) coordinates, in kilometers.  The TOD epoch is the
supplied `jd` time.

If geographical coordinates are given as arrays, a 3-tuple of arrays of the same
length is returned.

### Example ###

Obtain the ECI coordinates of the intersection of the equator and Greenwich's
meridian on 2015-06-30T14:03:12.857

```jldoctest
julia> using AstroLib

julia> geo2eci(0, 0, 0, jdcnv("2015-06-30T14:03:12.857"))
(-4024.8671780315185, 4947.835465127513, 0.0)
```

### Notes ###

`eci2geo` converts Earth-centered inertial coordinates to geographic spherical
coordinates.

Code of this function is based on IDL Astronomy User's Library.
"""
geo2eci(lat::Real, long::Real, alt::Real, jd::Real) =
    geo2eci(promote(float(lat), float(long), float(alt), float(jd))...)

geo2eci(lla::Tuple{Real, Real, Real}, jd::Real) =
    geo2eci(lla..., jd)

function geo2eci(lat::AbstractArray{LA}, long::AbstractArray{<:Real},
                 alt::AbstractArray{<:Real},
                 jd::AbstractArray{<:Real}) where {LA<:Real}
    if !(length(lat) == length(long) == length(alt) == length(jd))
        throw(DimensionMismatch(
            "lat, long, alt, and jd arrays must have the same length"))
    end
    typela = float(LA)
    x = similar(lat, typela)
    y = similar(lat, typela)
    z = similar(lat, typela)
    for i in eachindex(lat)
        x[i], y[i], z[i] = geo2eci(lat[i], long[i], alt[i], jd[i])
    end
    return x, y, z
end
