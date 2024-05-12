# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function _mag2geo(lat::T, long::T, pole_lat::T, pole_long::T) where {T<:AbstractFloat}
    r    = 1 # Distance from planet center.  Value unimportant -- just need
             # a length for conversion to rectangular coordinates
    sin_lat, cos_lat = sincos(deg2rad(lat))
    sin_long, cos_long = sincos(deg2rad(long))

    # convert to rectangular coordinates
    #   x-axis: defined by the vector going from Earth's center towards
    #        the intersection of the equator and Greenwich's meridian.
    #   z-axis: axis of the geographic poles
    #   y-axis: defined by y=z^x
    x = r * cos_lat * cos_long
    y = r * cos_lat * sin_long
    z = r * sin_lat

    # First rotation: in the plane of the current meridian from magnetic pole to
    # geographic pole.
    s, c = sincos(T(pi) / 2 - pole_lat)
    togeolat = SMatrix{3,3}(c,       zero(T),      -s,
                            zero(T),  one(T), zero(T),
                            s,       zero(T),       c)

    # Second rotation matrix: rotation around plane of the equator, from the
    # meridian containing the magnetic poles to the Greenwich meridian.
    sin_pole_long, cos_pole_long = sincos(pole_long)
    maglong2geolong = SMatrix{3,3}( cos_pole_long, sin_pole_long, zero(T),
                                   -sin_pole_long, cos_pole_long, zero(T),
                                    zero(T),             zero(T),  one(T))
    out = maglong2geolong * togeolat * SVector(x, y, z)

    geolat  = rad2deg(atan(out[3], hypot(out[1], out[2])))
    geolong = rad2deg(atan(out[2], out[1]))
    # I don't care about that one...just put it there for completeness' sake
    # magalt  = norm(out) - r
    return geolat, geolong
end

"""
    mag2geo(latitude, longitude[, year]) -> geographic_latitude, geographic_longitude

### Purpose ###

Convert from geomagnetic to geographic coordinates.

### Explanation ###

Converts from geomagnetic (latitude, longitude) to geographic (latitude,
longitude).  Altitude is not involved in this function.

### Arguments ###

* `latitude`: geomagnetic latitude (North), in degrees.
* `longitude`: geomagnetic longitude (East), in degrees.
* `year` (optional numerical argument): the year in which to perform conversion.
  If omitted, defaults to current year.

The coordinates can be passed as arrays of the same length.

### Output ###

The 2-tuple of geographic (latitude, longitude) coordinates, in degrees.

If geomagnetic coordinates are given as arrays, a 2-tuple of arrays of the same
length is returned.

### Example ###

Find position of North Magnetic Pole in 2016

```jldoctest
julia> using AstroLib

julia> mag2geo(90, 0, 2016)
(86.395, -166.29000000000002)
```

### Notes ###

This function uses list of North Magnetic Pole positions provided by World
Magnetic Model (https://www.ngdc.noaa.gov/geomag/data/poles/NP.xy).

`geo2mag` converts geographic coordinates to geomagnetic coordinates.

Code of this function is based on IDL Astronomy User's Library.
"""
mag2geo(lat::Real, long::Real, year::Real=Dates.year(Dates.now())) =
    _mag2geo(promote(float(lat), float(long),
                     deg2rad(POLELATLONG[year][1]::Float64),
                     deg2rad(POLELATLONG[year][2]::Float64))...)

function mag2geo(lat::AbstractArray{LA}, long::AbstractArray{LO},
                 year::Real=Dates.year(Dates.now())) where {LA<:Real, LO<:Real}
    if length(lat) != length(long)
        throw(ArgumentError("lat and long arrays should be of the same length"))
    end
    typela   = float(LA)
    geolat   = similar(lat, typela)
    geolong  = similar(lat, typela)
    polelat  = deg2rad(POLELATLONG[year][1])
    polelong = deg2rad(POLELATLONG[year][2])
    for i in eachindex(lat)
        geolat[i], geolong[i] =
            _mag2geo(promote(float(lat[i]), float(long[i]),
                             polelat, polelong)...)
    end
    return geolat, geolong
end
