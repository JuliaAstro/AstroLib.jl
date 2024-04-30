# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

const d_lng = SVector(0, -2, 0, 0, 0, 0, -2, 0, 0, -2, -2, -2, 0, 2, 0, 2, 0, 0, -2, 0,
                      2, 0, 0, -2, 0, -2, 0, 0, 2,  -2, 0, -2, 0, 0, 2, 2, 0, -2, 0, 2,
                      2, -2, -2, 2, 2, 0, -2, -2, 0, -2, -2, 0, -1, -2, 1, 0, 0, -1, 0,
                      0,  2, 0, 2)
const M_lng = SVector(0, 0, 0, 0, 1, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 2, 0, 2, 1, 0, -1, 0, 0, 0, 1, 1, -1, 0, 0, 0,
                      0, 0, 0, -1, -1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, -1, 1, -1, -1, 0,
                      -1)
const Mprime_lng = SVector(0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, -1, 0, 1, -1, -1, 1, 2,
                           -2, 0, 2, 2, 1, 0, 0, -1, 0, -1,  0, 0, 1, 0, 2, -1, 1, 0,
                           1, 0, 0, 1, 2, 1, -2, 0, 1, 0, 0, 2, 2, 0, 1, 1, 0, 0, 1,
                           -2, 1, 1, 1, -1, 3, 0)
const F_lng = SVector(0, 2, 2, 0, 0, 0, 2, 2, 2, 2, 0, 2, 2, 0, 0, 2, 0, 2, 0, 2, 2, 2,
                      0, 2, 2, 2, 2, 0, 0, 2, 0, 0,  0, -2, 2, 2, 2, 0, 2, 2, 0, 2, 2,
                      0, 0, 0, 2, 0, 2, 0, 2, -2, 0, 0, 0, 2, 2, 0, 0, 2, 2, 2, 2)
const ω_lng = SVector(1, 2, 2, 2, 0, 0, 2, 1, 2, 2, 0, 1, 2, 0, 1, 2, 1, 1, 0, 1, 2, 2,
                      0, 2, 0, 0, 1, 0, 1, 2, 1, 1, 1, 0, 1, 2, 2, 0, 2, 1, 0, 2, 1, 1,
                      1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 2, 2, 2, 2)
const sin_lng = SVector(-171996, -13187, -2274, 2062, 1426, 712, -517, -386, -301, 217,
                        -158, 129, 123, 63, 63, -59, -58, -51, 48, 46, -38, -31, 29,
                        29, 26, -22, 21, 17, 16, -16, -15, -13, -12, 11, -10, -8, 7,
                        -7, -7, -7, 6, 6, 6, -6, -6, 5, -5, -5, -5, 4, 4, 4, -4, -4,
                        -4, 3, -3, -3, -3, -3, -3, -3, -3)
const sdelt = SVector(-174.2, -1.6, -0.2, 0.2, -3.4, 0.1, 1.2, -0.4, 0, -0.5, 0, 0.1,
                      0.0, 0.0, 0.1, 0.0, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, 0.0, -0.1, 0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
const cos_lng = SVector(92025, 5736, 977, -895, 54, -7, 224, 200, 129, -95, 0, -70,
                        -53, 0, -33, 26, 32, 27, 0, -24, 16, 13, 0, -12, 0, 0, -10, 0,
                        -8, 7, 9, 7, 6, 0, 5, 3, -3, 0, 3, 3, 0, -3, -3, 3, 3, 0, 3, 3,
                        3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
const cdelt = SVector(8.9, -3.1, -0.5, 0.5, -0.1, 0.0, -0.6, 0.0, -0.1, 0.3, 0.0, 0.0,
                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

function nutate(jd::AbstractFloat)
    # Number of Julian centuries since 2000-01-01T12:00:00
    t = (jd - J2000) / JULIANCENTURY
    # Mean elongation of the Moon
    d = deg2rad(mod(@evalpoly(t, 297.85036, 445267.111480, -0.0019142, inv(189474)), 360))
    # Sun's mean anomaly
    M = deg2rad(mod(@evalpoly(t, 357.52772, 35999.050340, -0.0001603, -inv(3e5)), 360))
    # Moon's mean anomaly
    Mprime = deg2rad(mod(@evalpoly(t, 134.96298, 477198.867398, 0.0086972, inv(5.625e4)), 360))
    # Moon's argument of latitude
    F = deg2rad(mod(@evalpoly(t, 93.27191, 483202.017538, -0.0036825, -inv(3.27270e5)), 360))
    # Longitude of the ascending node of the Moon's mean orbit on the ecliptic,
    # measured from the mean equinox of the date
    ω = deg2rad(mod(@evalpoly(t, 125.04452, -1934.136261, 0.0020708, inv(4.5e5)), 360))
    arg = d_lng * d + M_lng * M + Mprime_lng * Mprime + F_lng * F + ω_lng * ω
    sc_arg = sincos.(arg)
    s = getindex.(sc_arg, 1)
    c = getindex.(sc_arg, 2)
    long  = dot((sdelt .* t .+ sin_lng), s) / 10000
    obliq = dot((cdelt .* t .+ cos_lng), c) / 10000
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

Uses the formula in Chapter 22 of "Astronomical Algorithms" by Jean Meeus
(1998, 2nd ed.) which is based on the 1980 IAU Theory of Nutation and includes
all terms larger than 0.0003".

### Example ###

- Find the nutation in longitude and obliquity 1987 on Apr 10 at 0h.  This is
  example 22.a from Meeus

  ```jldoctest
  julia> using AstroLib

  julia> jd = jdcnv(1987, 4, 10);

  julia> nutate(jd)
  (-3.787931077110494, 9.44252069864449)
  ```

- Plot the daily nutation in longitude and obliquity during the 21st century.
  Use [PyPlot.jl](https://github.com/JuliaPlots/Plots.jl/) for plotting.

  ```julia
  using PyPlot
  years = DateTime(2000):DateTime(2100);
  long, obl = nutate(jdcnv.(years));
  plot(years, long); plot(years, obl)
  ```

You can see both the dominant large scale period of nutation, of 18.6 years, and
smaller oscillations with shorter periods.

### Notes ###

Code of this function is based on IDL Astronomy User's Library.
"""
nutate(jd::Real) = nutate(float(jd))

function nutate(jd::AbstractArray{J}) where {J<:Real}
    typej = float(J)
    long = similar(jd, typej)
    obliq = similar(jd, typej)
    for i in eachindex(jd)
        long[i], obliq[i] = nutate(jd[i])
    end
    return long, obliq
end
