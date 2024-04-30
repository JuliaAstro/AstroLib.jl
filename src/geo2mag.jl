# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function _geo2mag(lat::T, long::T, pole_lat::T, pole_long::T) where {T<:AbstractFloat}
    r    = 1 # Distance from planet center.  Value unimportant -- just need a length for
             # conversion to rectangular coordinates
    sin_lat, cos_lat = sincos(deg2rad(lat))
    sin_long, cos_long = sincos(deg2rad(long))
    x    = r * cos_lat * cos_long
    y    = r * cos_lat * sin_long
    z    = r * sin_lat

    # Compute first rotation matrix: rotation around plane of the equator, from
    # the Greenwich meridian to the meridian containing the magnetic dipole
    # pole.
    sin_pole_long, cos_pole_long = sincos(pole_long)
    geolong2maglong = SMatrix{3,3}(cos_pole_long, -sin_pole_long, zero(T),
                                   sin_pole_long,  cos_pole_long, zero(T),
                                   zero(T),              zero(T),  one(T))

    # Second rotation: in the plane of the current meridian from geographic pole
    # to magnetic dipole pole.
    s, c = sincos(T(pi) / 2 - pole_lat)
    tomaglat = SMatrix{3,3}( c,       zero(T),       s,
                             zero(T),  one(T), zero(T),
                            -s,       zero(T),       c)
    out = tomaglat * geolong2maglong * SVector(x, y, z)

    maglat  = rad2deg(atan(out[3], hypot(out[1], out[2])))
    maglong = rad2deg(atan(out[2], out[1]))
    # I don't care about that one...just put it there for completeness' sake
    # magalt  = norm(out) - r
    return maglat, maglong
end

"""
    geo2mag(latitude, longitude[, year]) -> geomagnetic_latitude, geomagnetic_longitude

### Purpose ###

Convert from geographic to geomagnetic coordinates.

### Explanation ###

Converts from geographic (latitude, longitude) to geomagnetic (latitude,
longitude).  Altitude is not involved in this function.

### Arguments ###

* `latitude`: geographic latitude (North), in degrees.
* `longitude`: geographic longitude (East), in degrees.
* `year` (optional numerical argument): the year in which to perform conversion.
  If omitted, defaults to current year.

The coordinates can be passed as arrays of the same length.

### Output ###

The 2-tuple of magnetic (latitude, longitude) coordinates, in degrees.

If geographical coordinates are given as arrays, a 2-tuple of arrays of the same
length is returned.

### Example ###

Kyoto has geographic coordinates 35° 00' 42'' N, 135° 46' 06'' E, find its
geomagnetic coordinates in 2016:

```jldoctest
julia> using AstroLib

julia> geo2mag(ten(35,0,42), ten(135,46,6), 2016)
(36.86579228937769, -60.1840605366516)
```

### Notes ###

This function uses list of North Magnetic Pole positions provided by World
Magnetic Model (https://www.ngdc.noaa.gov/geomag/data/poles/NP.xy).

`mag2geo` converts geomagnetical coordinates to geographic coordinates.

Code of this function is based on IDL Astronomy User's Library.
"""
geo2mag(lat::Real, long::Real, year::Real=Dates.year(Dates.now())) =
    _geo2mag(promote(float(lat), float(long),
                     deg2rad(POLELATLONG[year][1]::AbstractFloat),
                     deg2rad(POLELATLONG[year][2]::AbstractFloat))...)

function geo2mag(lat::AbstractArray{LA}, long::AbstractArray{<:Real},
                 year::Real=Dates.year(Dates.now())) where {LA<:Real}
    @assert length(lat) == length(long)
    typela   = float(LA)
    maglat   = similar(lat, typela)
    maglong  = similar(lat, typela)
    polelat  = deg2rad(POLELATLONG[year][1])
    polelong = deg2rad(POLELATLONG[year][2])
    for i in eachindex(lat)
        maglat[i], maglong[i] =
            _geo2mag(promote(float(lat[i]), float(long[i]),
                             polelat, polelong)...)
    end
    return maglat, maglong
end
