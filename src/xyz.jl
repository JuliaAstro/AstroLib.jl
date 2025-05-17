# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function _xyz{T<:AbstractFloat}(jd::T, equinox::T)
    t = (float(jd) - 15020.0)*inv(JULIANYEAR*100) # Reduced Julian century from 1900

    # NOTE: longitude arguments below are given in *equinox* of date.  Precess
    # these to equinox 1950 to give everything an even footing.  Compute
    # argument of precession from equinox of date back to 1950
    pp = (1.396041 + 0.000308*(t + 0.5))*(t-0.499998)
    # Compute mean solar longitude, precessed back to 1950
    el = 279.696678 + 36000.76892*t + 0.000303*t*t - pp

    # Compute Mean longitude of the Moon
    c = 270.434164 + 480960.0*t + 307.883142*t - 0.001133*t*t - pp

    # Compute longitude of Moon's ascending node
    n = 259.183275 - 1800.0*t - 134.142008*t + 0.002078*t*t - pp

    # Compute mean solar anomaly
    g = 358.475833 + 35999.04975*t - 0.00015*t*t

    # Compute the mean jupiter anomaly
    j = 225.444651 + 2880.0*t + 154.906654*t*t

    # Compute mean anomaly of Venus
    v = 212.603219 + 58320.0*t + 197.803875*t + 0.001286*t*t

    # Compute mean anomaly of Mars
    m = 319.529425 + 19080.0*t + 59.8585*t + 0.000181*t*t

    # Convert degrees to radians for trig functions
    el = deg2rad(el)
    g  = deg2rad(g)
    j  = deg2rad(j)
    c  = deg2rad(c)
    v  = deg2rad(v)
    n  = deg2rad(n)
    m  = deg2rad(m)

    # Calculate x, y, z using trigonometric series
    x = 0.999860*cos(el) -
        0.025127*cos(g - el) +
        0.008374*cos(g + el) +
        0.000105*cos(g + g + el) +
        0.000063*t*cos(g - el) +
        0.000035*cos(g + g - el) -
        0.000026*sin(g - el - j) -
        0.000021*t*cos(g + el) +
        0.000018*sin(2.0*g + el - 2.0*v) +
        0.000017*cos(c) -
        0.000014*cos(c - 2.0*el) +
        0.000012*cos(4.0*g + el - 8.0*m + 3.0*j) -
        0.000012*cos(4.0*g - el - 8.0*m + 3.0*j) -
        0.000012*cos(g + el - v) +
        0.000011*cos(2.0*g + el - 2.0*v) +
        0.000011*cos(2.0*g - el - 2.0*j)

    y = 0.917308*sin(el) +
        0.023053*sin(g - el) +
        0.007683*sin(g + el) +
        0.000097*sin(g + g + el) -
        0.000057*t*sin(g - el) -
        0.000032*sin(g + g - el) -
        0.000024*cos(g - el - j) -
        0.000019*t*sin(g + el) -
        0.000017*cos(2.0*g + el - 2.0*v) +
        0.000016*sin(c) +
        0.000013*sin(c - 2.0*el ) +
        0.000011*sin(4.0*g + el - 8.0*m + 3.0*j) +
        0.000011*sin(4.0*g - el - 8.0*m + 3.0*j) -
        0.000011*sin(g + el - v) +
        0.000010*sin(2.0*g + el - 2.0*v ) -
        0.000010*sin(2.0*g - el - 2.0*j )

    z = 0.397825*sin(el)     +
        0.009998*sin(g-el)   +
        0.003332*sin(g+el)   +
        0.000042*sin(g+g+el) -
        0.000025*t*sin(g-el) -
        0.000014*sin(g+g-el) -
        0.000010*cos(g-el-j)

    # Precess to new equator?  Avoid useless calculations.
    if isfinite(equinox) && equinox != 1950
        x, y, z = precess_xyz(x, y, z, 1950.0, equinox)
    end

    # Velocities
    xvel = -0.017200 * sin(el)   -
        0.000288*sin(g + el)     -
        0.000005*sin(2.0*g + el) -
        0.000004*sin(c)          +
        0.000003*sin(c - 2.0*el) +
        0.000001*t*sin(g + el)   -
        0.000001*sin(2.0*g - el)

    yvel = 0.015780*cos(el)       +
        0.000264*cos(g + el)      +
        0.000005*cos(2.0*g + el)  +
        0.000004*cos(c)           +
        0.000003*cos(c - 2.0*el)  -
        0.000001*t*cos(g + el)

    zvel = 0.006843*cos(el)      +
        0.000115*cos(g  + el)    +
        0.000002*cos(2.0*g + el) +
        0.000002*cos(c)          +
        0.000001*cos(c - 2.0*el)

    # Precess to new equator?  Avoid useless calculations.
    if isfinite(equinox) && equinox != 1950
        xvel, yvel, zvel = precess_xyz(xvel, yvel, zvel, 1950, equinox)
    end

    return x, y, z, xvel, yvel, zvel
end

"""
    xyz(jd[, equinox]) -> x, y, z, v_x, v_y, v_z

### Purpose ###

Calculate geocentric \$x\$, \$y\$, and \$z\$ and velocity coordinates of the
Sun.

### Explanation ###

Calculates geocentric \$x\$, \$y\$, and \$z\$ vectors and velocity coordinates
(\$dx\$, \$dy\$ and \$dz\$) of the Sun.  (The positive \$x\$ axis is directed
towards the equinox, the \$y\$-axis, towards the point on the equator at right
ascension 6h, and the \$z\$ axis toward the north pole of the equator).  Typical
position accuracy is \$<10^{-4}\$ AU (15000 km).

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

``` julia
jd = juldate(DateTime(1999, 1, 22))
# => 51200.5
xyz(jd + 64./86400., 2000)
# => (0.5145687092402946,-0.7696326261820777,-0.33376880143026394,0.014947267514081075,0.008314838205475709,0.003606857607574784)
```

Compare to Astronomical Almanac (1999 page C20)

                x  (AU)        y  (AU)     z (AU)
    xyz:      0.51456871   -0.76963263  -0.33376880
    AA:       0.51453130   -0.7697110   -0.3337152
    abs(err): 0.00003739    0.00007839   0.00005360
    abs(err)
        (km):   5609          11759         8040

NOTE: Velocities in AA are for Earth/Moon barycenter
      (a very minor offset) see AA 1999 page E3

               x vel (AU/day) y vel (AU/day)   z vel (AU/day)
    xyz:      -0.014947268   -0.0083148382    -0.0036068576
    AA:       -0.01494574    -0.00831185      -0.00360365
    abs(err):  0.000001583    0.0000029886     0.0000032076
    abs(err)
     (km/sec): 0.00265        0.00519          0.00557


### Notes ###

Code of this function is based on IDL Astronomy User's Library.
"""
xyz(jd::Real, equinox::Real=NaN) = _xyz(promote(float(jd), float(equinox))...)

# Can't use @vectorize_1arg because of the optional keyword.
function xyz{J<:Real}(jd::AbstractArray{J}, equinox::Real=NaN)
    typej = typeof(float(one(J)))
    x     = similar(jd, typej)
    y     = similar(jd, typej)
    z     = similar(jd, typej)
    xvel  = similar(jd, typej)
    yvel  = similar(jd, typej)
    zvel  = similar(jd, typej)
    for i in eachindex(jd)
        x[i], y[i], z[i], xvel[i], yvel[i], zvel[i] = xyz(jd[i], equinox)
    end
    return x, y, z, xvel, yvel, zvel
end
