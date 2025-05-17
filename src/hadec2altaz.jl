# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function _hadec2altaz{T<:AbstractFloat}(ha::T, dec::T, lat::T, ws::Bool)
    sh = sind(ha)
    ch = cosd(ha)
    sd = sind(dec)
    cd = cosd(dec)
    sl = sind(lat)
    cl = cosd(lat)

    x = -ch*cd*sl + sd*cl
    y = -sh*cd
    z = ch*cd*cl + sd*sl
    r = hypot(x, y)

    # Now get altitude, azimuth
    az  = cirrange(rad2deg(atan2(y, x)))
    alt = rad2deg(atan2(z, r))
    # Convert azimuth to West from South, if desired
    if ws
        az = cirrange(az + 180.0)
    end
    return alt, az
end

"""
    hadec2altaz(ha, dec, lat[, ws=true]) -> alt, az

### Purpose ###

Convert Hour Angle and Declination to Horizon (Alt-Az) coordinates.

### Explanation ###

Can deal with the NCP singularity.  Intended mainly to be used by program
`eq2hor`.

### Arguments ###

Input coordinates may be either a scalar or an array, of the same dimension.

* `ha`: the local apparent hour angle, in degrees.  The hour angle is the time
  that right ascension of 0 hours crosses the local meridian.  It is
  unambiguously defined.
* `dec`: the local apparent declination, in degrees.
* `lat`: the local geodetic latitude, in degrees, scalar or array.
* `ws` (optional boolean keyword): if true, the output azimuth is measured West
  from South.  The default is to measure azimuth East from North.

`ha` and `dec` can be given as a 2-tuple `(ha, dec)`.

### Output ###

2-tuple `(alt, az)`

* `alt`: local apparent altitude, in degrees.
* `az`: the local apparent azimuth, in degrees.

The output coordinates are always floating points and have the same type (scalar
or array) as the input coordinates.

### Example ###

Arcturus is observed at an apparent hour angle of 336.6829 and a declination of
19.1825 while at the latitude of +43° 4' 42''.  What are the local altitude and
azimuth of this object?

``` julia
alt, az = hadec2altaz(336.6829, 19.1825, ten(43, 4, 42))
# => (59.08617155005683,133.3080693440254)
```

### Notes ###

`altaz2hadec` converts Horizon (Alt-Az) coordinates to Hour Angle and
Declination.

Code of this function is based on IDL Astronomy User's Library.
"""
hadec2altaz(ha::Real, dec::Real, lat::Real; ws::Bool=false) =
    _hadec2altaz(promote(float(ha), float(dec), float(lat))..., ws)

hadec2altaz(hadec::Tuple{Real, Real}, lat::Real; ws::Bool=false) =
    hadec2altaz(hadec..., lat, ws=ws)

function hadec2altaz{R1<:Real, R2<:Real, R3<:Real}(ha::AbstractArray{R1},
                                                   dec::AbstractArray{R2},
                                                   lat::AbstractArray{R3};
                                                   ws::Bool=false)
    @assert length(ha) == length(dec) == length(lat)
    typeha = typeof(float(one(R1)))
    alt = similar(ha, typeha)
    az  = similar(ha, typeha)
    for i in eachindex(ha)
        alt[i], az[i] = hadec2altaz(ha[i], dec[i], lat[i], ws=ws)
    end
    return alt, az
end
