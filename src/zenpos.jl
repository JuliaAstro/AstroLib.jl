# This file is a part of AstroLib.jl. License is MIT "Expat".

"""
    zenpos(jd, latitude, longitude) -> zenith_right_ascension, declination
    zenpos(date, latitude, longitude, tz) -> zenith_right_ascension, declination

### Purpose ###

Return the zenith right ascension and declination in radians for a given Julian date
or a local civil time and timezone.

### Explanation ###

The local sidereal time is computed with the help of [`ct2lst`](@ref), which is the right
ascension of the zenith. This and the observatories latitude (corresponding to the declination)
are converted to radians and returned as the zenith direction.

### Arguments ###

The function can be called in two different ways. The arguments common to
both methods are `latitude` and `longitude`:

* `latitude` : latitude of the desired location.
* `longitude` : longitude of the desired location.

The zenith direction can be computed either by providing the Julian date:

* `jd` : the Julian date of the date and time for which the zenith position is
  desired.

or the time zone and the date:

* `tz`: the time zone (in hours) of the desired location (e.g. 4 = EDT, 5 = EST)
* `date`: the local civil time with type `DateTime`.

### Output ###

A 2-tuple `(ra, dec)`:

* `ra` : the right ascension (in radians) of the zenith.
* `dec` : the declination (in radians) of the zenith.

### Example ###

```jldoctest
julia> using AstroLib, Dates

julia> zenpos(DateTime(2017, 04, 25, 18, 59), 43.16, -24.32, 4)
(0.946790432684706, 0.7532841051607526)

julia> zenpos(jdcnv(2016, 05, 05, 13, 41), ten(35,0,42), ten(135,46,6))
(3.5757821152779536, 0.6110688599440813)
```

### Notes ###

Code of this function is based on IDL Astronomy User's Library.
"""
function zenpos end

_zenpos(latitude::T, longitude::T, rest...) where {T <: AbstractFloat} =
    (ct2lst(longitude, rest...) / 12 * pi, deg2rad(latitude))

zenpos(jd::Real, latitude::Real, longitude::Real) =
    _zenpos(promote(float(latitude), float(longitude), float(jd))...)

zenpos(date::DateTime, latitude::Real, longitude::Real, tz::Real) =
    _zenpos(promote(float(latitude), float(longitude), float(tz))..., date)
