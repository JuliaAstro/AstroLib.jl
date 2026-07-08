# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function _xyz(jd::T, equinox::T) where {T <: AbstractFloat}
    t = (jd - 15020) / JULIANCENTURY # Reduced Julian century from 1900

    # NOTE: longitude arguments below are given in *equinox* of date.  Precess
    # these to equinox 1950 to give everything an even footing.  Compute
    # argument of precession from equinox of date back to 1950
    pp = (1.396041 + 0.000308 * (t + 0.5)) * (t - 0.499998)
    # Compute mean solar longitude, precessed back to 1950
    el = @evalpoly t (279.696678 - pp) 36000.76892 0.000303

    # Compute Mean longitude of the Moon
    c = @evalpoly t (270.434164 - pp) (480960.0 + 307.883142) -0.001133

    # Compute longitude of Moon's ascending node
    n = @evalpoly t (259.183275 - pp) -(1800.0 + 134.142008) 0.002078

    # Compute mean solar anomaly
    g = @evalpoly t 358.475833 35999.04975 -0.00015

    # Compute the mean jupiter anomaly
    j = @evalpoly t 225.444651 2880.0 154.906654

    # Compute mean anomaly of Venus
    v = @evalpoly t 212.603219 (58320.0 + 197.803875) 0.001286

    # Compute mean anomaly of Mars
    m = @evalpoly t 319.529425 (19080.0 + 59.8585) 0.000181

    # Convert degrees to radians for trig functions
    el = deg2rad(el)
    g = deg2rad(g)
    j = deg2rad(j)
    c = deg2rad(c)
    v = deg2rad(v)
    n = deg2rad(n)
    m = deg2rad(m)

    # Calculate x, y, z using trigonometric series
    sin1, cos1 = sincos(el)
    sin2, cos2 = sincos(g - el)
    sin3, cos3 = sincos(g + el)
    sin4, cos4 = sincos(g + g + el)
    sin5, cos5 = sincos(g + g - el)
    sin6, cos6 = sincos(g - el - j)
    sin7, cos7 = sincos(2 * g + el - 2 * v)
    sin8, cos8 = sincos(c)
    sin9, cos9 = sincos(c - 2 * el)
    sin10, cos10 = sincos(4 * g + el - 8 * m + 3 * j)
    sin11, cos11 = sincos(4 * g - el - 8 * m + 3 * j)
    sin12, cos12 = sincos(g + el - v)
    sin13, cos13 = sincos(2 * g - el - 2 * j)

    x = 0.99986 * cos1 -
        0.025127 * cos2 +
        0.008374 * cos3 +
        0.000105 * cos4 +
        0.000063 * t * cos2 +
        0.000035 * cos5 -
        0.000026 * sin6 -
        0.000021 * t * cos3 +
        0.000018 * sin7 +
        0.000017 * cos8 -
        0.000014 * cos9 +
        0.000012 * cos10 -
        0.000012 * cos11 -
        0.000012 * cos12 +
        0.000011 * cos7 +
        0.000011 * cos13

    y = 0.917308 * sin1 +
        0.023053 * sin2 +
        0.007683 * sin3 +
        0.000097 * sin4 -
        0.000057 * t * sin2 -
        0.000032 * sin5 -
        0.000024 * cos6 -
        0.000019 * t * sin3 -
        0.000017 * cos7 +
        0.000016 * sin8 +
        0.000013 * sin9 +
        0.000011 * sin10 +
        0.000011 * sin11 -
        0.000011 * sin12 +
        0.00001 * sin7 -
        0.00001 * sin13

    z = 0.397825 * sin1 +
        0.009998 * sin2 +
        0.003332 * sin3 +
        0.000042 * sin4 -
        0.000025 * t * sin2 -
        0.000014 * sin5 -
        0.00001 * cos6

    # Precess to new equator?  Avoid useless calculations.
    if isfinite(equinox) && equinox != 1950
        x, y, z = precess_xyz(x, y, z, 1950.0, equinox)
    end

    # Velocities
    xvel = -0.0172 * sin(el) -
        0.000288 * sin(g + el) -
        0.000005 * sin(2.0 * g + el) -
        0.000004 * sin(c) +
        0.000003 * sin(c - 2.0 * el) +
        0.000001 * t * sin(g + el) -
        0.000001 * sin(2.0 * g - el)

    yvel = 0.01578 * cos(el) +
        0.000264 * cos(g + el) +
        0.000005 * cos(2.0 * g + el) +
        0.000004 * cos(c) +
        0.000003 * cos(c - 2.0 * el) -
        0.000001 * t * cos(g + el)

    zvel = 0.006843 * cos(el) +
        0.000115 * cos(g + el) +
        0.000002 * cos(2.0 * g + el) +
        0.000002 * cos(c) +
        0.000001 * cos(c - 2.0 * el)

    # Precess to new equator?  Avoid useless calculations.
    if isfinite(equinox) && equinox != 1950
        xvel, yvel, zvel = precess_xyz(xvel, yvel, zvel, 1950, equinox)
    end

    return x, y, z, xvel, yvel, zvel
end

"""
    xyz(jd[, equinox]) -> x, y, z, v_x, v_y, v_z

### Purpose ###

Calculate geocentric ``x``, ``y``, and ``z`` and velocity coordinates of the
Sun.

### Explanation ###

Calculates geocentric ``x``, ``y``, and ``z`` vectors and velocity coordinates
(``dx``, ``dy`` and ``dz``) of the Sun.  (The positive ``x`` axis is directed
towards the equinox, the ``y``-axis, towards the point on the equator at right
ascension 6h, and the ``z`` axis toward the north pole of the equator).  Typical
position accuracy is ``<10^{-4}`` AU (15000 km).

### Arguments ###

* `jd`: number of Reduced Julian Days for the wanted date.  It can be either a
  scalar or a vector.
* `equinox` (optional numeric argument): equinox of output. Default is 1950.

You can use `juldate` to get the number of Reduced Julian Days for the selected
dates.

### Output ###

The 6-tuple ``(x, y, z, v_x, v_y, v_z)``, where

* ``x, y, z``: scalars or vectors giving heliocentric rectangular
  coordinates (in AU) for each date supplied.  Note that ``\\sqrt{x^2 + y^2 + z^2}``
  gives the Earth-Sun distance for the given date.
* ``v_x, v_y, v_z``: velocity vectors corresponding to ``x, y``, and
  ``z``.

### Example ###

What were the rectangular coordinates and velocities of the Sun on
1999-01-22T00:00:00 (= JD 2451200.5) in J2000 coords?  Note: Astronomical
Almanac (AA) is in TDT, so add 64 seconds to UT to convert.

```jldoctest
julia> using AstroLib, Dates

julia> jd = juldate(DateTime(1999, 1, 22))
51200.5

julia> xyz(jd + 64/86400, 2000)
(0.514568709240398, -0.7696326261820209, -0.33376880143023935, 0.014947267514079971, 0.008314838205477328, 0.003606857607575486)
```

Compare to Astronomical Almanac (1999 page C20)

```plain
            x  (AU)        y  (AU)     z (AU)
xyz:      0.51456871   -0.76963263  -0.33376880
AA:       0.51453130   -0.7697110   -0.3337152
abs(err): 0.00003739    0.00007839   0.00005360
abs(err)
    (km):   5609          11759         8040
```

NOTE: Velocities in AA are for Earth/Moon barycenter
      (a very minor offset) see AA 1999 page E3
```plain
           x vel (AU/day) y vel (AU/day)   z vel (AU/day)
xyz:      -0.014947268   -0.0083148382    -0.0036068576
AA:       -0.01494574    -0.00831185      -0.00360365
abs(err):  0.000001583    0.0000029886     0.0000032076
abs(err)
 (km/sec): 0.00265        0.00519          0.00557
```

### Notes ###

Code of this function is based on IDL Astronomy User's Library.
"""
xyz(jd::Real, equinox::Real = NaN) = _xyz(promote(float(jd), float(equinox))...)

# Can't use @vectorize_1arg because of the optional keyword.
function xyz(jd::AbstractArray{J}, equinox::Real = NaN) where {J <: Real}
    typej = float(J)
    x = similar(jd, typej)
    y = similar(jd, typej)
    z = similar(jd, typej)
    xvel = similar(jd, typej)
    yvel = similar(jd, typej)
    zvel = similar(jd, typej)
    for i in eachindex(jd)
        x[i], y[i], z[i], xvel[i], yvel[i], zvel[i] = xyz(jd[i], equinox)
    end
    return x, y, z, xvel, yvel, zvel
end
