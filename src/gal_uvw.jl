# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

const A_G = SVector(0.0548755604, +0.4941094279, -0.8676661490,
                    0.8734370902, -0.4448296300, -0.1980763734,
                    0.4838350155,  0.7469822445, +0.4559837762)


function _gal_uvw(ra::T, dec::T, pmra::T, pmdec::T, vrad::T,
                  plx::T, lsr::Bool) where {T<:AbstractFloat}
    sindec, cosdec = sincos(deg2rad(dec))
    sinra, cosra = sincos(deg2rad(ra))
    k = 4.740470463533348 # = 149597870.7/(86400*365.25) = 1 AU/year in km/s
    vec2 = k*pmra/plx
    vec3 = k*pmdec/plx
    u = ( A_G[1]*cosra*cosdec + A_G[4]*sinra*cosdec + A_G[7]*sindec)*vrad +
        (-A_G[1]*sinra        + A_G[4]*cosra                       )*vec2 +
        (-A_G[1]*cosra*sindec - A_G[4]*sinra*sindec + A_G[7]*cosdec)*vec3
    v = ( A_G[2]*cosra*cosdec + A_G[5]*sinra*cosdec + A_G[8]*sindec)*vrad +
        (-A_G[2]*sinra        + A_G[5]*cosra                       )*vec2 +
        (-A_G[2]*cosra*sindec - A_G[5]*sinra*sindec + A_G[8]*cosdec)*vec3
    w = ( A_G[3]*cosra*cosdec + A_G[6]*sinra*cosdec + A_G[9]*sindec)*vrad +
        (-A_G[3]*sinra        + A_G[6]*cosra                       )*vec2 +
        (-A_G[3]*cosra*sindec - A_G[6]*sinra*sindec + A_G[9]*cosdec)*vec3
    if lsr
        u += -8.5
        v += 13.38
        w +=  6.49
    end
    return u, v, w
end

"""
    gal_uvw(ra, dec, pmra, pmdec, vrad, plx[, lsr=true]) -> u, v, w

### Purpose ###

Calculate the Galactic space velocity \$(u, v, w)\$ of a star.

### Explanation ###

Calculates the Galactic space velocity \$(u, v, w)\$ of a star given its (1)
coordinates, (2) proper motion, (3) parallax, and (4) radial velocity.

### Arguments ###

User must supply a position, proper motion, radial velocity and parallax.
Either scalars or arrays all of the same length can be supplied.

(1) Position:

* `ra`: right ascension, in degrees
* `dec`: declination, in degrees

(2) Proper Motion

* `pmra`: proper motion in right ascension in arc units (typically
  milli-arcseconds/yr).  If given \$\\mu_\\alpha\$ -- proper motion in seconds of
  time/year -- then this is equal to \$15 \\mu_\\alpha \\cos(\\text{dec})\$.
* `pmdec`: proper motion in declination (typically mas/yr).

(3) Radial Velocity

* `vrad`: velocity in km/s

(4) Parallax

* `plx`: parallax with same distance units as proper motion measurements
  typically milliarcseconds (mas)

If you know the distance in parsecs, then set `plx` to \$1000/\\text{distance}\$,
if proper motion measurements are given in milli-arcseconds/yr.

There is an additional optional keyword:

* `lsr` (optional boolean keyword): if this keyword is set to `true`, then the
  output velocities will be corrected for the solar motion
  \$(u, v, w)_\\odot = (-8.5, 13.38, 6.49)\$ (Coşkunoǧlu et al. 2011 MNRAS, 412,
  1237;
  DOI:[10.1111/j.1365-2966.2010.17983.x](http://dx.doi.org/10.1111/j.1365-2966.2010.17983.x))
  to the local standard of rest (LSR).  Note that the value of the solar motion
  through the LSR remains poorly determined.

### Output ###

The 3-tuple \$(u, v, w)\$

* ``u``: velocity (km/s) positive toward the Galactic *anti*center
* ``v``: velocity (km/s) positive in the direction of Galactic rotation
* ``w``: velocity (km/s) positive toward the North Galactic Pole

### Method ###

Follows the general outline of Johnson & Soderblom (1987, AJ, 93, 864;
DOI:[10.1086/114370](http://dx.doi.org/10.1086/114370)) except that \$u\$ is
positive outward toward the Galactic *anti*center, and the J2000 transformation
matrix to Galactic coordinates is taken from the introduction to the Hipparcos
catalog.

### Example ###

Compute the U,V,W coordinates for the halo star HD 6755.  Use values from
Hipparcos catalog, and correct to the LSR.

```jldoctest
julia> using AstroLib

julia> ra=ten(1,9,42.3)*15.; dec = ten(61,32,49.5);

julia> pmra = 627.89;  pmdec = 77.84; # mas/yr

julia> vrad = -321.4; dis = 129; # distance in parsecs

julia> u, v, w = gal_uvw(ra, dec, pmra, pmdec, vrad, 1e3/dis, lsr=true)
(118.2110474553902, -466.4828898385057, 88.16573278565097)
```

### Notes ###

This function does not take distance as input.  See "Arguments" section above
for how to provide it using parallax argument.

Code of this function is based on IDL Astronomy User's Library.
"""
gal_uvw(ra::Real, dec::Real, pmra::Real, pmdec::Real,
        vrad::Real, plx::Real; lsr::Bool=false) =
            _gal_uvw(promote(float(ra), float(dec), float(pmra), float(pmdec),
                             float(vrad), float(plx))..., lsr)

function gal_uvw(ra::AbstractArray{R}, dec::AbstractArray{<:Real},
                 pmra::AbstractArray{<:Real}, pmdec::AbstractArray{<:Real},
                 vrad::AbstractArray{<:Real}, plx::AbstractArray{<:Real};
                 lsr::Bool=false) where {R<:Real}
    if !(length(ra) == length(dec) == length(pmra) ==
         length(pmdec) == length(vrad) == length(plx))
        throw(ArgumentError(
            "ra, dec, pmra, pmdec, vrad, and plx arrays must all have the same length"))
    end
    typer = float(R)
    u = similar(ra, typer)
    v = similar(ra, typer)
    w = similar(ra, typer)
    for i in eachindex(ra)
        u[i], v[i], w[i] = gal_uvw(ra[i], dec[i], pmra[i], pmdec[i],
                                   vrad[i], plx[i], lsr=lsr)
    end
    return u, v, w
end
