# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function sunpos(jd::Real; radians::Bool = false)
    # Number of Julian centuries since 1899-12-31T12:00:00
    t = (jd - 2415020) / JULIANCENTURY
    # Sun's mean longitude
    l = (279.696678 + mod(36000.768925*t, 360)) * 3600
    # Allow for ellipticity of the orbit (equation of centre) using the Earth's
    # mean anomaly ME
    me = 358.475844 + mod(35999.049750*t, 360)
    ellcor = (6910.1 - 17.2 * t) * sind(me) + 72.3 * sind(2 * me)
    l += ellcor
    # Allow for the Venus perturbations using the mean anomaly of Venus MV
    mv = 212.603219 + mod(58517.803875*t, 360)
    vencorr = 4.8*cosd(299.1017 + mv - me) +
        5.5*cosd(148.3133 +  2 * mv - 2 * me) +
        2.5*cosd(315.9433 +  2 * mv - 3 * me) +
        1.6*cosd(345.2533 +  3 * mv - 4 * me) +
        1.0*cosd(318.15   +  3 * mv - 5 * me)
    l += vencorr
    # Allow for the Mars perturbations using the mean anomaly of Mars MM.
    mm = 319.529425  +  mod(19139.858500*t, 360)
    marscorr = 2 * cosd(343.8883 - 2 * mm + 2 * me) +
        1.8 * cosd(200.4017 - 2 * mm + me)
    l += marscorr
    # Allow for the Jupiter perturbations using the mean anomaly of Jupiter MJ
    mj = 225.328328 + mod(3034.6920239 * t, 360)
    jupcorr = 7.2*cosd(179.5317 - mj + me) +
        2.6*cosd(263.2167 - mj) +
        2.7*cosd( 87.1450 - 2 * mj + 2 * me) +
        1.6*cosd(109.4933 - 2 * mj + me)
    l += jupcorr
    # Allow for the Moons perturbations using the mean elongation of the Moon
    # from the Sun D
    d = 350.7376814 + mod(445267.11422 * t, 360)
    mooncorr  = 6.5*sind(d)
    l += mooncorr
    # Allow for long period terms
    longterm = 6.4*sind(231.19 + 20.20*t)
    l += longterm
    l  = mod(l + 2592000, 1296000)
    longmed = sec2rad(l)
    # Allow for Aberration
    l -= 20.5
    # Allow for Nutation using the longitude of the Moons mean node OMEGA
    sin_omega, cos_omega = sincos(deg2rad(259.183275 - mod(1934.142008 * t, 360)))
    l -= 17.2 * sin_omega
    # True Obliquity
    oblt = 23.452294 - 0.0130125 * t + 9.2 * cos_omega / 3600
    # Right Ascension and Declination
    l /= 3600
    oblt = deg2rad(oblt)
    sin_oblt, cos_oblt = sincos(oblt)
    sin_l, cos_l = sincos(deg2rad(l))
    ra = mod2pi(atan(sin_l * cos_oblt, cos_l))
    dec = asin(sin_l * sin_oblt)
    if radians
        return ra, dec, longmed, oblt
    else
        return rad2deg(ra), rad2deg(dec), rad2deg(longmed), rad2deg(oblt)
    end
end

"""
    sunpos(jd[, radians=true]) -> ra, dec, elong, obliquity

### Purpose ###

Compute the right ascension and declination of the Sun at a given date.

### Arguments ###

* `jd`: the Julian date of when you want to calculate Sun position.  It can be
  either a scalar or a vector.  Use `jdcnv` to get the Julian date for a given
  date and time.
* `radians` (optional boolean keyword): if set to `true`, all output quantities
  are given in radians.  The default is `false`, so all quantities are given in
  degrees.

### Output ###

The 4-tuple `(ra, dec, elong, obliquity)`:

* `ra`: the right ascension of the Sun at that date
* `dec`: the declination of the Sun at that date
* `elong`: ecliptic longitude of the Sun at that date
* `obliquity`: the obliquity of the ecliptic

All quantities are given in degrees, unless `radians` keyword is set to `true`
(see "Arguments" section).  If `jd` is an array, arrays of the same given as
`jd` are returned.

### Method ###

Uses a truncated version of Newcomb's Sun.  Adapted from the IDL routine SUN_POS
by CD Pike, which was adapted from a FORTRAN routine by B. Emerson (RGO).

### Example ###

(1) Find the apparent right ascension and declination of the Sun on May 1, 1982

```jldoctest
julia> using AstroLib

julia> adstring(sunpos(jdcnv(1982, 5, 1))[1:2], precision=2)
" 02 31 32.614  +14 54 34.92"
```

The Astronomical Almanac gives `02 31 32.58 +14 54 34.9` so the error for this
case is < 0.5".

(2) Plot the apparent right ascension, in hours, and declination of the Sun, in
degrees, for every day in 2016.  Use
[PyPlot.jl](https://github.com/JuliaPlots/Plots.jl/) for plotting.

```julia
using PyPlot
using Dates

days = DateTime(2016):Day(1):DateTime(2016, 12, 31);
ra, declin = sunpos(jdcnv.(days));
plot(days, ra/15); plot(days, declin)
```

### Notes ###

Patrick Wallace (Rutherford Appleton Laboratory, UK) has tested the accuracy of
a C adaptation of the present algorithm and found the following results.  From
1900-2100 `sunpos` gave 7.3 arcsec maximum error, 2.6 arcsec RMS.  Over the
shorter interval 1950-2050 the figures were 6.4 arcsec max, 2.2 arcsec RMS.

The returned `ra` and `dec` are in the given date's equinox.

Code of this function is based on IDL Astronomy User's Library.
"""

function sunpos(jd::AbstractVector{J}; radians::Bool=false) where {J<:Real}
    typej = float(J)
    ra = similar(jd, typej)
    dec = similar(jd, typej)
    longmed = similar(jd, typej)
    oblt = similar(jd, typej)
    for i in eachindex(jd)
        ra[i], dec[i], longmed[i], oblt[i] = sunpos(jd[i], radians=radians)
    end
    return ra, dec, longmed, oblt
end
