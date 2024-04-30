# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function _premat(equinox1::T, equinox2::T, FK4::Bool) where {T<:AbstractFloat}
    t = (equinox2 - equinox1) / 1000
    if FK4
        st = (equinox1 - 1900) / 1000
        #  Compute 3 rotation angles
        a = sec2rad(t*(23042.53 + st*(139.75 + 0.06*st)
                       + t*(30.23 - 0.27*st + 18.0*t)))
        b = sec2rad(t*t*(79.27 + 0.66*st + 0.32*t)) + a
        c = sec2rad(t*(20046.85 - st*(85.33 + 0.37*st)
                       + t*(-42.67 - 0.37*st -41.8*t)))
    else
        st = (equinox1 / 1000) - 2
        # Compute 3 rotation angles
        a = sec2rad(t *(23062.181 + st*(139.656 + 0.0139*st)
                        + t*(30.188 - 0.344*st + 17.998*t)))
        b = sec2rad(t*t*(79.280 + 0.410*st + 0.205*t)) + a
        c = sec2rad(t*(20043.109 - st*(85.33 + 0.217*st)
                       + t*(-42.665 - 0.217*st - 41.833*t)))
    end
    sa, ca = sincos(a)
    sb, cb = sincos(b)
    sc, cc = sincos(c)
    return @SMatrix [ca * cb * cc - sa * sb   -ca * sb - sa * cb * cc   -cb*sc;
                     sa * cb + ca * sb * cc    ca * cb - sa * sb * cc   -sb*sc;
                     ca * sc                  -sa * sc                   cc]
end

"""
    premat(equinox1, equinox2[, FK4=true]) -> precession_matrix

### Purpose ###

Return the precession matrix needed to go from `equinox1` to `equinox2`.

### Explanation ###

This matrix is used by `precess` and `baryvel` to precess astronomical
coordinates.

### Arguments ###

* `equinox1`: original equinox of coordinates.
* `equinox2`: equinox of precessed coordinates.
* `FK4` (optional boolean keyword): if this keyword is set to `true`, the FK4
  (B1950.0) system precession angles are used to compute the precession matrix.
  When it is `false`, the default, use FK5 (J2000.0) precession angles.

### Output ###

A 3×3 matrix, used to precess equatorial rectangular coordinates.

### Example ###

Return the precession matrix from 1950.0 to 1975.0 in the FK4 system

```jldoctest
julia> using AstroLib

julia> premat(1950, 1975, FK4=true)
3×3 StaticArraysCore.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 0.999981    -0.00558775  -0.00242909
 0.00558775   0.999984    -6.78691e-6
 0.00242909  -6.78633e-6   0.999997
```

### Method ###

FK4 constants from "Computational Spherical Astronomy" by Taff (1983),
p. 24. (FK4). FK5 constants from "Explanatory Supplement To The Astronomical
Almanac" 1992, page 104 Table 3.211.1
(https://archive.org/details/131123ExplanatorySupplementAstronomicalAlmanac).

### Notes ###

Code of this function is based on IDL Astronomy User's Library.
"""
premat(eq1::Real, eq2::Real; FK4::Bool=false) =
    _premat(promote(float(eq1), float(eq2))..., FK4)
