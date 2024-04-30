# This file is a part of AstroLib.jl. License is MIT "Expat".

function _hor2eq(alt::T, az::T, jd::T, lat::T, lon::T, altitude::T,
                 pressure::T, temperature::T, ws::Bool, B1950::Bool,
                 precession::Bool, nutate::Bool, aberration::Bool,
                 refract::Bool) where {T<:AbstractFloat}

    if refract
        alt = co_refract(alt, altitude, pressure, temperature)
    end

    if ws
        az -= 180
    end
    _, _, eps, d_psi = co_nutate(jd, 45, 45)
    last = 15 * ct2lst(lon, jd) + d_psi * cos(eps) / 3600
    ha, dec = altaz2hadec(alt, az, lat)
    ra = mod(last - ha, 360)
    dra1, ddec1, eps = co_nutate(jd, ra, dec)

    if aberration
        dra2, ddec2 = co_aberration(jd, ra, dec, eps)
        ra -= dra2 / 3600
        dec -= ddec2 / 3600
    end

    if nutate
       ra -= dra1
       dec -= ddec1
    end
    j_now = (jd - J2000) / JULIANYEAR + 2000

    if precession
        if B1950
            ra, dec = precess(ra, dec, j_now, 1950, FK4=true)
        else
            ra, dec = precess(ra, dec, j_now, 2000)
        end
    end
    return ra, dec, ha
end

hor2eq(alt::Real, az::Real, jd::Real, lat::Real, lon::Real, altitude::Real=0;
       ws::Bool=false, B1950::Bool=false, precession::Bool=true, nutate::Bool=true,
       aberration::Bool=true, refract::Bool=true, pressure::Real=NaN,
       temperature::Real=NaN) =
           _hor2eq(promote(float(alt), float(az), float(jd), float(lat), float(lon),
                   float(altitude), float(temperature), float(pressure))..., ws, B1950,
                   precession, nutate, aberration, refract)

hor2eq(alt::Real, az::Real, jd::Real, obsname::AbstractString; kwargs...) =
    hor2eq(alt, az, jd, observatories[obsname].latitude, observatories[obsname].longitude,
           observatories[obsname].altitude; kwargs...)


"""
    hor2eq(alt, az, jd[, obsname;
           ws=false, B1950=false, precession=true, nutate=true,
           aberration=true, refract=true,
           lat=NaN, lon=NaN, altitude=0,
           pressure=NaN, temperature=NaN]) -> ra, dec, ha

    hor2eq(alt, az, jd, lat, lon[, altitude=0;
           ws=false, B1950=false, precession=true, nutate=true,
           aberration=true, refract=true,
           pressure=NaN, temperature=NaN]) -> ra, dec, ha


### Purpose

Converts local horizon coordinates (alt-az) to equatorial (ra-dec) coordinates.

### Explanation

This is a function to calculate equatorial (ra,dec) coordinates from
horizon (alt,az) coords. It is accurate to about 1 arcsecond or better.
It performs precession, nutation, aberration, and refraction corrections.

### Arguments

This function has two base methods.  With one you can specify the name of the observatory,
if present in `AstroLib.observatories`, with the other one you can provide the coordinates
of the observing site and, optionally, the altitude.

Common mandatory arguments:

* `alt`: altitude of horizon coords, in degrees
* `az`: azimuth angle measured East from North (unless ws is `true`), in degrees
* `jd`: julian date

Other positional arguments:

* `obsname`: set this to a valid observatory name in `AstroLib.observatories`.

or

* `lat`: north geodetic latitude of location, in degrees.
* `lon`: AST longitude of location, in degrees. You can specify west longitude with a
  negative sign.
* `altitude`: the altitude of the observing location, in meters.  It is `0` by default

Optional keyword arguments:

* `ws` (optional boolean keyword): set this to `true` to get the azimuth measured
  westward from south. This is `false` by default
* `B1950` (optional boolean keyword): set this to `true` if the ra and dec
  are specified in B1950 (FK4 coordinates) instead of J2000 (FK5). This is `false` by
  default
* `precession` (optional boolean keyword): set this to `false` for no precession,
  `true` by default
* `nutate` (optional boolean keyword): set this to `false` for no nutation,
  `true` by default
* `aberration` (optional boolean keyword): set this to `false` for no aberration
  correction, `true` by default
* `refract` (optional boolean keyword): set this to `false` for no refraction
  correction, `true` by default
* `pressure` (optional keyword): the pressure at the observing location, in millibars.
  Default value is `NaN`
* `temperature` (optional keyword): the temperature at the observing location, in Kelvins.
  Default value is `NaN`

### Output

* `ra`: right ascension of object, in degrees (FK5)
* `dec`: declination of the object, in degrees (FK5)
* `ha`: hour angle, in degrees

### Example

You are at Kitt Peak National Observatory, looking at a star at azimuth
angle 264d 55m 06s and elevation 37d 54m 41s (in the visible). Today is Dec 25, 2041
and the local time is 10 PM precisely. What is the right ascension and declination
(J2000) of the star you're looking at? The temperature here is about 0 Celsius,
and the pressure is 781 millibars. The Julian date for this time is 2466879.7083333

```jldoctest
julia> using AstroLib

julia> ra_o, dec_o = hor2eq(ten(37,54,41), ten(264,55,06), 2466879.7083333,
                            "kpno", pressure = 781, temperature = 273)
(3.3224480269254713, 15.19061543702944, 54.61174536229464)

julia> adstring(ra_o, dec_o)
" 00 13 17.4  +15 11 26"
```

### Notes

Code of this function is based on IDL Astronomy User's Library.
"""
function hor2eq end

# TODO: Make hor2eq type-stable, which it isn't currently because of keyword arguments
# Note that the inner function `_hor2eq` is type stable
