# This file is a part of AstroLib.jl. License is MIT "Expat".

#TODO: Use full JPL ephemeris for high precision
function _planet_coords(date::T, num::Integer) where {T<:AbstractFloat}
    rad, long, lat = helio(date, num, true)
    rade, longe, late = helio(date, 3, true)
    sin_lat, cos_lat = sincos(lat)
    sin_long, cos_long = sincos(long)
    sin_late, cos_late = sincos(late)
    sin_longe, cos_longe = sincos(longe)
    x = rad * cos_lat * cos_long - rade * cos_late * cos_longe
    y = rad * cos_lat * sin_long - rade * cos_late * sin_longe
    z = rad * sin_lat - rade * sin_late
    lamda = rad2deg(atan(y,x))
    beta = rad2deg(atan(z, hypot(x,y)))
    ra, dec = euler(lamda, beta, 4)
    return ra, dec
end

"""
    planet_coords(date, num)

### Purpose ###

Find right ascension and declination for the planets when provided a date as input.

### Explanation ###

This function uses the [`helio`](@ref) to get the heliocentric ecliptic coordinates of the
planets at the given date which it then converts these to geocentric ecliptic
coordinates ala "Astronomical Algorithms" by Jean Meeus (1991, p 209).
These are then converted to right ascension and declination using [`euler`](@ref).

The accuracy between the years 1800 and 2050 is better than 1 arcminute for the
terrestial planets, but reaches 10 arcminutes for Saturn. Before 1850 or after 2050
the accuracy can get much worse.

### Arguments ###

* `date`: Can be either a single date or an array of dates. Each element can be
  either a `DateTime` type or Julian Date. It can be a scalar or vector.
* `num`: integer denoting planet number, scalar or vector
  1 = Mercury, 2 = Venus, ... 9 = Pluto. If not in that change, then the
  program will throw an error.

### Output ###

* `ra`: right ascension of planet(J2000), in degrees
* `dec`: declination of the planet(J2000), in degrees

### Example ###

Find the RA, Dec of Venus on 1992 Dec 20

```jldoctest
julia> using AstroLib, Dates

julia> adstring(planet_coords(DateTime(1992,12,20),2))
" 21 05 02.8  -18 51 41"
```

### Notes ###

Code of this function is based on IDL Astronomy User's Library.
"""
planet_coords(date::Real, num::Integer) = _planet_coords(float(date), num)

function planet_coords(date::AbstractVector{R},
                       num::AbstractVector{<:Integer}) where {R<:Real}
    if length(date) != length(num)
      throw(ArgumentError("date and num arrays should be of the same length"))
    end
    typedate = float(R)
    ra_out  = similar(date, typedate)
    dec_out = similar(date, typedate)
    for i in eachindex(date)
        ra_out[i], dec_out[i] = planet_coords(date[i], num[i])
    end
    return ra_out, dec_out
end

planet_coords(date::DateTime, num::Integer) = planet_coords(jdcnv(date), num)
