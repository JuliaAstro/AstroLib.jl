# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function _posang{T<:AbstractFloat}(units::Integer, ra1::T, dec1::T, ra2::T, dec2::T)
    # Convert all quantities to radians.
    if units == 0
        # All radians
        ra1_rad  = ra1
        ra2_rad  = ra2
        dec1_rad = dec1
        dec2_rad = dec2
    elseif units == 1
        # Right ascensions are in hours, declinations in degrees.
        ra1_rad  = ra1*pi/12.0
        ra2_rad  = ra2*pi/12.0
        dec1_rad = deg2rad(dec1)
        dec2_rad = deg2rad(dec2)
    elseif units == 2
        # Right ascensions and declinations are in degrees.
        ra1_rad  = deg2rad(ra1)
        ra2_rad  = deg2rad(ra2)
        dec1_rad = deg2rad(dec1)
        dec2_rad = deg2rad(dec2)
    else
        # In any other case throw an error.
        error("units must be 0 (radians), 1 (hours, degrees) or 2 (degrees)")
    end
    radif = ra2_rad - ra1_rad
    angle = atan2(sin(radif), cos(dec1_rad)*tan(dec2_rad) -
                  sin(dec1_rad)*cos(radif))
    if units == 0
        return angle
    else
        return rad2deg(angle)
    end
end

"""
    posang(units, ra1, dec1, ra2, dec2) -> angular_distance

### Purpose ###

Compute rigorous position angle of point 2 relative to point 1.

### Explanation ###

Computes the rigorous position angle of point 2 (with given right ascension and
declination) using point 1 (with given right ascension and declination) as the
center.

### Arguments ###

* `units`: integer, can be either 0, or 1, or 2.  Describes units of inputs and
 output:
    * 0: everything (input right ascensions and declinations, and output
      distance) is radians
    * 1: right ascensions are in decimal hours, declinations in decimal degrees,
      output distance in degrees
    * 2: right ascensions and declinations are in degrees, output distance in
      degrees
* `ra1`:  right ascension or longitude of point 1
* `dec1`: declination or latitude of point 1
* `ra2`: right ascension or longitude of point 2
* `dec2`: declination or latitude of point 2

Both `ra1` and `dec1`, and `ra2` and `dec2` can be given as 2-tuples `(ra1,
dec1)` and `(ra2, dec2)`.

### Output ###

Angle of the great circle containing `[ra2, dec2]` from the meridian containing
`[ra1, dec1]`, in the sense north through east rotating about `[ra1, dec1]`.
See `units` argument above for units.

### Method ###

The "four-parts formula" from spherical trigonometry (p. 12 of Smart's Spherical
Astronomy or p. 12 of Green' Spherical Astronomy).


### Example ###

Mizar has coordinates (ra, dec) = (13h 23m 55.5s, +54° 55' 31'').  Its
companion, Alcor, has coordinates (ra, dec) = (13h 25m 13.5s, +54° 59' 17'').
Find the position angle of Alcor with respect to Mizar.

``` julia
posang(1, ten(13, 25, 13.5), ten(54, 59, 17), ten(13, 23, 55.5), ten(54, 55, 31))
# => -108.46011246802047
```

### Notes ###

* If `ra1`, `dec1` are scalars, and `ra2`, `dec2` are vectors, then the output
 is a vector giving the distance of each element of `ra2`, `dec2` to `ra1`,
 `dec1`.  Similarly, if `ra1`,`de1` are vectors, and `ra2`,` dec2` are scalars,
 then the output is a vector giving the distance of each element of `ra1`,
 `dec1` to `ra2`, `dec2`.  If both `ra1`, `dec1` and `ra2`, `dec2` are vectors
 then the output is a vector giving the distance of each element of `ra1`,
 `dec1` to the corresponding element of `ra2`, `dec2`.
* The function `sphdist` provides an alternate method of computing a spherical
 distance.
* Note that `posang` is not commutative: the position angle between A and B is
  \$\\theta\$, then the position angle between B and A is \$180 + \\theta\$.

Code of this function is based on IDL Astronomy User's Library.
"""
posang(units::Integer, ra1::Real, dec1::Real, ra2::Real, dec2::Real) =
    _posang(units, promote(float(ra1), float(dec1), float(ra2), float(dec2))...)

### Tuples input
posang(units::Integer, radec1::Tuple{Real, Real}, ra2::Real, dec2::Real) =
    posang(units, radec1..., ra2, dec2)

posang(units::Integer, ra1::Real, dec1::Real, radec2::Tuple{Real, Real}) =
    posang(units, ra1, dec1, radec2...)

posang(units::Integer, radec1::Tuple{Real, Real}, radec2::Tuple{Real, Real}) =
    posang(units, radec1..., radec2...)
