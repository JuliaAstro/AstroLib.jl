# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

const d_lng = [0, -2, 0, 0, 0, 0, -2, 0, 0, -2, -2, -2, 0, 2, 0, 2, 0, 0, -2, 0,
               2, 0, 0, -2, 0, -2, 0, 0, 2,  -2, 0, -2, 0, 0, 2, 2, 0, -2, 0, 2,
               2, -2, -2, 2, 2, 0, -2, -2, 0, -2, -2, 0, -1, -2, 1, 0, 0, -1, 0,
               0,  2, 0, 2]
const M_lng = [0, 0, 0, 0, 1, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 2, 0, 2, 1, 0, -1, 0, 0, 0, 1, 1, -1, 0, 0, 0,
               0, 0, 0, -1, -1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, -1, 1, -1, -1, 0,
               -1]
const Mprime_lng = [0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, -1, 0, 1, -1, -1, 1, 2,
                    -2, 0, 2, 2, 1, 0, 0, -1, 0, -1,  0, 0, 1, 0, 2, -1, 1, 0,
                    1, 0, 0, 1, 2, 1, -2, 0, 1, 0, 0, 2, 2, 0, 1, 1, 0, 0, 1,
                    -2, 1, 1, 1, -1, 3, 0]
const F_lng = [0, 2, 2, 0, 0, 0, 2, 2, 2, 2, 0, 2, 2, 0, 0, 2, 0, 2, 0, 2, 2, 2,
               0, 2, 2, 2, 2, 0, 0, 2, 0, 0,  0, -2, 2, 2, 2, 0, 2, 2, 0, 2, 2,
               0, 0, 0, 2, 0, 2, 0, 2, -2, 0, 0, 0, 2, 2, 0, 0, 2, 2, 2, 2]
const ω_lng = [1, 2, 2, 2, 0, 0, 2, 1, 2, 2, 0, 1, 2, 0, 1, 2, 1, 1, 0, 1, 2, 2,
               0, 2, 0, 0, 1, 0, 1, 2, 1, 1, 1, 0, 1, 2, 2, 0, 2, 1, 0, 2, 1, 1,
               1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 2, 2, 2, 2]
const sin_lng = [-171996, -13187, -2274, 2062, 1426, 712, -517, -386, -301, 217,
                 -158, 129, 123, 63, 63, -59, -58, -51, 48, 46, -38, -31, 29,
                 29, 26, -22, 21, 17, 16, -16, -15, -13, -12, 11, -10, -8, 7,
                 -7, -7, -7, 6, 6, 6, -6, -6, 5, -5, -5, -5, 4, 4, 4, -4, -4,
                 -4, 3, -3, -3, -3, -3, -3, -3, -3]
const sdelt = [-174.2, -1.6, -0.2, 0.2, -3.4, 0.1, 1.2, -0.4, 0, -0.5, 0, 0.1,
               0.0, 0.0, 0.1, 0.0, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, -0.1, 0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
const cos_lng = [92025, 5736, 977, -895, 54, -7, 224, 200, 129, -95, 0, -70,
                 -53, 0, -33, 26, 32, 27, 0, -24, 16, 13, 0, -12, 0, 0, -10, 0,
                 -8, 7, 9, 7, 6, 0, 5, 3, -3, 0, 3, 3, 0, -3, -3, 3, 3, 0, 3, 3,
                 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
const cdelt = [8.9, -3.1, -0.5, 0.5, -0.1, 0.0, -0.6, 0.0, -0.1, 0.3, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

function _nutate{T<:AbstractFloat}(jd::T)
    # Number of Julian centuries since 2000-01-01T12:00:00
    t = (jd - J2000)*inv(JULIANYEAR*100)
    # Mean elongation of the Moon
    d = deg2rad(cirrange(@evalpoly(t, 297.85036, 445267.111480, -0.0019142, inv(189474))))
    # Sun's mean anomaly
    M = deg2rad(cirrange(@evalpoly(t, 357.52772, 35999.050340, -0.0001603, -inv(3e5))))
    # Moon's mean anomaly
    Mprime = deg2rad(cirrange(@evalpoly(t, 134.96298, 477198.867398, 0.0086972, inv(5.625e4))))
    # Moon's argument of latitude
    F = deg2rad(cirrange(@evalpoly(t, 93.27191, 483202.017538, -0.0036825, -inv(3.27270e5))))
    # Longitude of the ascending node of the Moon's mean orbit on the ecliptic,
    # measured from the mean equinox of the date
    ω = deg2rad(cirrange(@evalpoly(t, 125.04452, -1934.136261, 0.0020708, inv(4.5e5))))
    arg = d_lng*d + M_lng*M + Mprime_lng*Mprime + F_lng*F + ω_lng*ω
    long  = 0.0001*sum((sdelt*t + sin_lng).*sin.(arg))
    obliq = 0.0001*sum((cdelt*t + cos_lng).*cos.(arg))
    return long, obliq
end

"""
    nutate(jd) -> long, obliq

### Purpose ###

Return the nutation in longitude and obliquity for a given Julian date.

### Arguments ###

* `jd`: Julian ephemeris date, it can be either a scalar or a vector

### Output ###

The 2-tuple `(long, obliq)`, where

* `long`: the nutation in longitude
* `obl`: the nutation in latitude

If `jd` is an array, `long` and `obl` are arrays of the same length.

### Method ###

Uses the formula in Chapter 22 of ``Astronomical Algorithms'' by Jean Meeus
(1998, 2nd ed.) which is based on the 1980 IAU Theory of Nutation and includes
all terms larger than 0.0003".

### Example ###

(1) Find the nutation in longitude and obliquity 1987 on Apr 10 at Oh.  This is
example 22.a from Meeus

``` julia
jd = jdcnv(1987, 4, 10);
nutate(jd)
# => (-3.787931077110755,9.442520698644401)
```

(2) Plot the daily nutation in longitude and obliquity during the 21st century.
Use [PyPlot.jl](https://github.com/stevengj/PyPlot.jl) for plotting.

``` julia
using PyPlot
years = DateTime(2000):DateTime(2100);
long, obl = nutate(jdcnv(years));
plot(years, long); plot(years, obl)
```

You can see both the dominant large scale period of nutation, of 18.6 years, and
smaller oscillations with shorter periods.

### Notes ###

Code of this function is based on IDL Astronomy User's Library.
"""
nutate(jd::Real) = _nutate(float(jd))

function nutate{J<:Real}(jd::AbstractArray{J})
    typej = float(J)
    long = similar(jd, typej)
    obliq = similar(jd, typej)
    for i in eachindex(jd)
        long[i], obliq[i] = nutate(jd[i])
    end
    return long, obliq
end
