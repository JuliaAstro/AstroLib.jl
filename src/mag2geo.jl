# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function _mag2geo{T<:AbstractFloat}(lat::T, long::T, pole_lat::T, pole_long::T)
    r       = 1.0 # Distance from planet center.  Value unimportant -- just need
                  # a length for conversion to rectangular coordinates
    lat  = deg2rad(lat)
    long = deg2rad(long)
    alt  = lat*0.0 + r

    # convert to rectangular coordinates
    #   x-axis: defined by the vector going from Earth's center towards
    #        the intersection of the equator and Greenwich's meridian.
    #   z-axis: axis of the geographic poles
    #   y-axis: defined by y=z^x
    x = alt*cos(lat)*cos(long)
    y = alt*cos(lat)*sin(long)
    z = alt*sin(lat)

    # First rotation: in the plane of the current meridian from magnetic pole to
    # geographic pole.
    togeolat = Array{T}(3, 3)
    togeolat[:,1] = [cos(pi/2 - pole_lat), 0.0, -sin(pi/2 - pole_lat)]
    togeolat[:,2] = [                 0.0, 1.0,                   0.0]
    togeolat[:,3] = [sin(pi/2 - pole_lat), 0.0,  cos(pi/2 - pole_lat)]
    out = togeolat * [x, y, z]

    # Second rotation matrix: rotation around plane of the equator, from the
    # meridian containing the magnetic poles to the Greenwich meridian.
    maglong2geolong = Array{T}(3, 3)
    maglong2geolong[:,1] = [ cos(pole_long), sin(pole_long), 0.0]
    maglong2geolong[:,2] = [-sin(pole_long), cos(pole_long), 0.0]
    maglong2geolong[:,3] = [            0.0,            0.0, 1.0]
    out = maglong2geolong * out

    geolat  = rad2deg(atan2(out[3], hypot(out[1], out[2])))
    geolong = rad2deg(atan2(out[2], out[1]))
    # I don't care about that one...just put it there for completeness' sake
    # magalt  = vecnorm(out) - r
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

``` julia
mag2geo(90, 0, 2016)
# => (86.395,-166.29000000000002)
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

function mag2geo{LA<:Real, LO<:Real}(lat::AbstractArray{LA},
                                     long::AbstractArray{LO},
                                     year::Real=Dates.year(Dates.now()))
    @assert length(lat) == length(long)
    typela   = typeof(float(one(LA)))
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
