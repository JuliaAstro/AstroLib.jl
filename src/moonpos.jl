# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

const moon_d_lng = [0, 2, 2, 0, 0, 0, 2, 2, 2, 2, 0, 1, 0, 2, 0, 0, 4, 0, 4, 2,
                    2, 1, 1, 2, 2, 4, 2, 0, 2, 2, 1, 2, 0, 0, 2, 2, 2, 4, 0, 3,
                    2, 4, 0, 2, 2, 2, 4, 0, 4, 1, 2, 0, 1, 3, 4, 2, 0, 1, 2, 2]
const moon_M_lng = [0, 0, 0, 0, 1, 0, 0, -1, 0, -1, 1, 0, 1, 0, 0, 0, 0, 0, 0,
                    1, 1, 0, 1, -1, 0, 0, 0, 1, 0, -1, 0, -2, 1, 2, -2, 0, 0,
                    -1, 0, 0, 1, -1, 2, 2, 1, -1, 0, 0, -1, 0, 1, 0, 1, 0, 0,
                    -1, 2, 1, 0, 0]
const moon_Mprime_lng = [1, -1, 0, 2, 0, 0, -2, -1, 1, 0, -1, 0, 1, 0, 1, 1, -1,
                         3, -2, -1, 0, -1, 0, 1, 2, 0, -3, -2, -1, -2, 1, 0, 2,
                         0, -1, 1, 0, -1, 2, -1, 1, -2, -1, -1, -2, 0, 1, 4, 0,
                         -2, 0, 2, 1, -2, -3, 2, 1, -1, 3, -1]
const moon_F_lng = [0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, -2, 2, -2, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, -2, 2, 0,
                    2, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, -2, -2, 0, 0, 0, 0, 0,
                    0, 0, -2]
const moon_sin_lng = [6288774, 1274027, 658314, 213618, -185116, -114332, 58793,
                      57066, 53322, 45758, -40923, -34720, -30383, 15327,
                      -12528, 10980, 10675, 10034, 8548, -7888, -6766,
                      -5163, 4987, 4036, 3994, 3861, 3665, -2689, -2602, 2390,
                      -2348, 2236, -2120, -2069, 2048, -1773, -1595, 1215,
                      -1110, -892, -810, 759, -713, -700, 691, 596, 549, 537,
                      520, -487, -399, -381, 351, -340, 330, 327, -323, 299, 294,
                      0]*1.0
const moon_cos_lng = [-20905355, -3699111, -2955968, -569925, 48888, -3149,
                      246158, -152138, -170733, -204586, -129620, 108743,
                      104755, 10321, 0, 79661, -34782, -23210, -21636, 24208,
                      30824, -8379, -16675, -12831, -10445, -11650, 14403,
                      -7003, 0, 10056, 6322, -9884, 5751, 0, -4950, 4130, 0,
                      -3958, 0, 3258, 2616, -1897, -2117, 2354, 0, 0, -1423,
                      -1117, -1571, -1739, 0, -4421, 0, 0, 0, 0, 1165, 0, 0,
                      8752]*1.0
const moon_d_lat = [0, 0, 0, 2, 2, 2, 2, 0, 2, 0, 2, 2, 2, 2, 2, 2, 2, 0, 4, 0,
                    0, 0, 1, 0, 0, 0, 1, 0, 4, 4, 0, 4, 2, 2, 2, 2, 0, 2, 2, 2,
                    2, 4, 2, 2, 0, 2, 1, 1, 0, 2, 1, 2, 0, 4, 4, 1, 4, 1, 4, 2]
const moon_M_lat = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 1, -1, -1, -1, 1, 0,
                    1, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0,
                    1, 1, 0, -1, -2, 0, 1, 1, 1, 1, 1, 0, -1, 1, 0, -1, 0, 0, 0,
                    -1, -2]
const moon_Mprime_lat = [0, 1, 1, 0, -1, -1, 0, 2, 1, 2, 0, -2, 1, 0, -1, 0, -1,
                         -1, -1, 0, 0, -1, 0, 1, 1, 0, 0, 3, 0, -1, 1,  -2, 0,
                         2, 1, -2, 3, 2, -3, -1, 0, 0, 1, 0, 1, 1, 0, 0, -2, -1,
                         1, -2, 2, -2, -1, 1, 1, -1, 0, 0]
const moon_F_lat =[ 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, -1, -1, 1, -1, 1, 1, -1,
                    -1, -1, 1, 3, 1, 1, 1, -1, -1, -1, 1, -1, 1, -3, 1, -3, -1,
                    -1, 1, -1, 1, -1, 1, 1, 1, 1, -1, 3, -1, -1, 1, -1, -1, 1,
                    -1, 1, -1, -1, -1, -1, -1, -1, 1]
const moon_sin_lat = [5128122, 280602, 277693, 173237, 55413, 46271, 32573,
                      17198, 9266, 8822, 8216, 4324, 4200, -3359, 2463, 2211,
                      2065, -1870, 1828, -1794, -1749, -1565, -1491, -1475,
                      -1410, -1344, -1335, 1107, 1021, 833, 777, 671, 607, 596,
                      491, -451, 439, 422, 421, -366, -351, 331, 315, 302, -283,
                      -229, 223, 223, -220, -220, -185, 181, -177, 176, 166,
                      -164, 132, -119, 115, 107]*1.0

function _moonpos{T<:AbstractFloat}(jd::T, radians::Bool)
    # Number of Julian centuries since 2000-01-01T12:00:00
    t = (jd - J2000)*inv(JULIANYEAR*100)
    # Mean longitude of the moon referred to mean equinox of the date
    Lprimed = cirrange(@evalpoly(t, 218.3164477, 481267.88123421,
                                 -0.0015786, inv(538841), -inv(6.5194e7)))
    Lprime = deg2rad(Lprimed)
    # Mean elongation of the Moon
    d = deg2rad(cirrange(@evalpoly(t, 297.8501921, 445267.1114034, -0.0018819,
                                   inv(545868), -inv(1.13065e8))))
    # Sun's mean anomaly
    M = deg2rad(cirrange(@evalpoly(t, 357.5291092, 35999.0502909, -0.0001536,
                                   inv(2.449e7))))
    # Moon's mean anomaly
    Mprime = deg2rad(cirrange(@evalpoly(t, 134.9633964, 477198.8675055,
                                        0.0087414, inv(6.9699e4),
                                        -inv(1.4712e7))))
    # Moon's argument of latitude
    F = deg2rad(cirrange(@evalpoly(t, 93.2720950, 483202.0175233, -0.0036539,
                                   -inv(3.526e7), inv(8.6331e8))))
    # Eccentricity of Earth's orbit around the Sun
    E = 1.0 - 0.002516*t - 7.4e-6*t*t
    E2 = E*E
    # Additional arguments
    A1 = deg2rad(119.75 + 131.849*t)
    A2 = deg2rad(53.09 + 479264.290*t)
    A3 = deg2rad(313.45 + 481266.484*t)
    suml_add = 3958*sin(A1) + 1962*sin(Lprime - F) + 318*sin(A2)
    sumb_add = -2235*sin(Lprime) + 382*sin(A3) + 175*sin(A1 - F) +
        175*sin(A1 + F) + 127*sin(Lprime - Mprime) -
        115*sin(Lprime + Mprime)
    # Sum the period terms
    sinlng = copy(moon_sin_lng)
    coslng = copy(moon_cos_lng)
    sinlat = copy(moon_sin_lat)
    for i in eachindex(moon_M_lng)
        if abs(moon_M_lng[i]) == 1
            sinlng[i] = E*sinlng[i]
            coslng[i] = E*coslng[i]
        elseif abs(moon_M_lng[i]) == 2
            sinlng[i] = E2*sinlng[i]
            coslng[i] = E2*coslng[i]
        end
        if abs(moon_M_lat[i]) == 1
            sinlat[i] = E*sinlat[i]
        elseif abs(moon_M_lat[i]) == 2
            sinlat[i] = E2*sinlat[i]
        end
    end
    arg = moon_d_lng*d + moon_M_lng*M + moon_Mprime_lng*Mprime + moon_F_lng*F
    geolong = Lprimed + (sum(sinlng.*sin.(arg)) + suml_add)/1e6
    dis = 385000.56 + sum(coslng.*cos.(arg))/1e3
    arg = moon_d_lat*d + moon_M_lat*M + moon_Mprime_lat*Mprime + moon_F_lat*F
    geolat = (sum(sinlat.*sin.(arg)) + sumb_add)/1e6
    nlong, elong = nutate(jd)
    geolong = cirrange(geolong + nlong/3.6e3)
    λ = deg2rad(geolong)
    β = deg2rad(geolat)
    # Find mean obliquity and convert λ, β to right ascension and declination.
    ɛ = ten(23, 26) + @evalpoly(t/1e2, 21.448, -4680.93, -1.55, 1999.25, -51.38,
                                -249.67, -39.05, 7.12, 27.87, 5.79, 2.45)/3.6e3
    ɛ = deg2rad(ɛ + elong/3.6e3)
    ra = cirrange(atan2(sin(λ)*cos(ɛ) - tan(β)*sin(ɛ), cos(λ)), 2.*pi)
    dec = asin(sin(β)*cos(ɛ) + cos(β)*sin(ɛ)*sin(λ))
    if radians
        return ra, dec, dis, λ, β
    else
        return rad2deg(ra), rad2deg(dec), dis, geolong, geolat
    end
end

"""
    moonpos(jd[, radians=true]) -> ra, dec, dis, geolong, geolat

### Purpose ###

Compute the right ascension and declination of the Moon at specified Julian
date.

### Arguments ###

* `jd`: the Julian ephemeris date.  It can be either a scalar or an array
* `radians` (optional boolean keyword): if set to `true`, then all output
  angular quantities are given in radians rather than degrees.  The default is
  `false`

### Output ###

The 5-tuple `(ra, dec, dis, geolong, geolat)`:

* `ra`: apparent right ascension of the Moon in degrees, referred to the true
  equator of the specified date(s)
* `dec`: the declination of the Moon in degrees
* `dis`: the distance between the centre of the Earth and the centre of the Moon
  in kilometers
* `geolong`: apparent longitude of the moon in degrees, referred to the ecliptic
  of the specified date(s)
* `geolat`: apparent longitude of the moon in degrees, referred to the ecliptic
  of the specified date(s)

If `jd` is an array, then all output quantities are arrays of the same length as
`jd`.

### Method ###

Derived from the Chapront ELP2000/82 Lunar Theory (Chapront-Touze' and Chapront,
1983, 124, 50), as described by Jean Meeus in Chapter 47 of ``Astronomical
Algorithms'' (Willmann-Bell, Richmond), 2nd edition, 1998.  Meeus quotes an
approximate accuracy of 10" in longitude and 4" in latitude, but he does not
give the time range for this accuracy.

Comparison of the IDL procedure with the example in ``Astronomical Algorithms''
reveals a very small discrepancy (~1 km) in the distance computation, but no
difference in the position calculation.

### Example ###

(1) Find the position of the moon on April 12, 1992

``` julia
jd = jdcnv(1992, 4, 12);
adstring(moonpos(jd)[1:2],precision=1)
# => " 08 58 45.23  +13 46 06.1"
```

This is within 1" from the position given in the Astronomical Almanac.

(2) Plot the Earth-moon distance during 2016 with sampling of 6 hours.  Use
[PyPlot.jl](https://github.com/stevengj/PyPlot.jl) for plotting

``` julia
using PyPlot
points = DateTime(2016):Dates.Hour(6):DateTime(2017);
plot(points, moonpos(jdcnv(points))[3])
```

### Notes ###

Code of this function is based on IDL Astronomy User's Library.
"""
moonpos(jd::Real; radians::Bool=false) = _moonpos(float(jd), radians)

function moonpos{J<:Real}(jd::AbstractArray{J}; radians::Bool=false)
    typej = float(J)
    ra = similar(jd, typej)
    dec = similar(jd, typej)
    dis = similar(jd, typej)
    long = similar(jd, typej)
    lat = similar(jd, typej)
    for i in eachindex(jd)
        ra[i], dec[i], dis[i], long[i], lat[i] = moonpos(jd[i], radians=radians)
    end
    return ra, dec, dis, long, lat
end
