# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function _precess(ra::T, dec::T, equinox1::T, equinox2::T,
                  FK4::Bool, radians::Bool) where {T<:AbstractFloat}
    if !radians
        # Convert to radians.
        ra  = deg2rad(ra)
        dec = deg2rad(dec)
    end
    sin_dec, cos_dec = sincos(dec)
    sin_ra, cos_ra = sincos(ra)
    x = SVector(cos_dec * cos_ra, cos_dec * sin_ra, sin_dec)
    r = premat(equinox1, equinox2, FK4=FK4)
    x2 = r*x
    ra_out  = mod2pi(atan(x2[2], x2[1]))
    dec_out = asin(x2[3])
    if radians
        return ra_out, dec_out
    else
        return rad2deg(ra_out), rad2deg(dec_out)
    end
end

"""
    precess(ra, dec, equinox1, equinox2[, FK4=true, radians=true]) -> prec_ra, prec_dec

### Purpose ###

Precess coordinates from `equinox1` to `equinox2`.

### Explanation ###

The default `(ra, dec)` system is FK5 based on epoch J2000.0 but FK4 based on
B1950.0 is available via the `FK4` boolean keyword.

### Arguments ###

* `ra`: input right ascension, scalar or vector, in degrees, unless the
  `radians` keyword is set to `true`
* `dec`: input declination, scalar or vector, in degrees, unless the `radians`
  keyword is set to `true`
* `equinox1`: original equinox of coordinates, numeric scalar.
* `equinox2`: equinox of precessed coordinates.
* `FK4` (optional boolean keyword): if this keyword is set to `true`, the FK4
  (B1950.0) system precession angles are used to compute the precession matrix.
  When it is `false`, the default, use FK5 (J2000.0) precession angles.
* `radians` (optional boolean keyword): if this keyword is set to `true`, then
  the input and output right ascension and declination vectors are in radians
  rather than degrees.

### Output ###

The 2-tuple `(ra, dec)` of coordinates modified by precession.

### Example ###

The Pole Star has J2000.0 coordinates (2h, 31m, 46.3s, 89d 15' 50.6"); compute
its coordinates at J1985.0

```jldoctest
julia> using AstroLib

julia> ra, dec = ten(2,31,46.3)*15, ten(89,15,50.6)
(37.94291666666666, 89.26405555555556)

julia> adstring(precess(ra, dec, 2000, 1985), precision=1)
" 02 16 22.73  +89 11 47.3"
```

Precess the B1950 coordinates of Eps Ind (RA = 21h 59m,33.053s, DEC = (-56d,
59', 33.053") to equinox B1975.

```jldoctest
julia> using AstroLib

julia> ra, dec = ten(21, 59, 33.053) * 15, ten(-56, 59, 33.053)
(329.88772083333333, -56.992514722222225)

julia> adstring(precess(ra, dec, 1950, 1975, FK4=true), precision=1)
" 22 01 15.46  -56 52 18.7"
```

### Method ###

Algorithm from "Computational Spherical Astronomy" by Taff (1983),
p. 24. (FK4). FK5 constants from "Explanatory Supplement To The Astronomical
Almanac" 1992, page 104 Table 3.211.1
(https://archive.org/details/131123ExplanatorySupplementAstronomicalAlmanac).

### Notes ###

Accuracy of precession decreases for declination values near 90 degrees.
`precess` should not be used more than 2.5 centuries from 2000 on the FK5 system
(1950.0 on the FK4 system).  If you need better accuracy, use `bprecess` or
`jprecess` as needed.

Code of this function is based on IDL Astronomy User's Library.
"""
precess(ra::Real, dec::Real, equinox1::Real, equinox2::Real;
        FK4::Bool=false, radians::Bool=false) =
            _precess(promote(float(ra), float(dec),
                             float(equinox1), float(equinox2))...,
                     FK4, radians)

precess(radec::Tuple{Real, Real}, equinox1::Real, equinox2::Real;
        FK4::Bool=false,
        radians::Bool=false) =
            precess(radec[1], radec[2], equinox1, equinox2,
                    FK4=FK4, radians=radians)

function precess(ra::AbstractArray{R}, dec::AbstractArray{D}, equinox1::Real,
                 equinox2::Real; FK4::Bool=false,
                 radians::Bool=false) where {R<:Real, D<:Real}
    if length(ra) != length(dec)
        throw(ArgumentError("ra and dec arrays should be of the same length"))
    end
    typera = float(R)
    ra_out  = similar(ra,  typera)
    dec_out = similar(dec, typera)
    for i in eachindex(ra)
        ra_out[i], dec_out[i] = precess(ra[i], dec[i], equinox1, equinox2,
                                        FK4=FK4, radians=radians)
    end
    return ra_out, dec_out
end
