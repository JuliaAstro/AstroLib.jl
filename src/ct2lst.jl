# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function _ct2lst{T<:AbstractFloat}(long::T, jd::T)
    t0 = jd - J2000
    t  = t0*inv(JULIANYEAR*100)
    # Compute GST in seconds.
    θ = ct2lst_c[1] + (ct2lst_c[2]*t0) + t*t*(ct2lst_c[3] - t*inv(ct2lst_c[4]))
    return cirrange((θ + long)/15, 24)
end

"""
    ct2lst(longitude, jd) -> local_sidereal_time
    ct2lst(longitude, tz, date) -> local_sidereal_time

### Purpose ###

Convert from Local Civil Time to Local Mean Sidereal Time.

### Arguments ###

The function can be called in two different ways.  The only argument common to
both methods is `longitude`:

* `longitude`: the longitude in degrees (east of Greenwich) of the place for
  which the local sidereal time is desired, scalar.  The Greenwich mean sidereal
  time (GMST) can be found by setting longitude = `0`.

The civil date to be converted to mean sidereal time can be specified either by
providing the Julian days:

* `jd`: this is number of Julian days for the date to be converted.  It can be a
  scalar or an array.

or the time zone and the date:

* `tz`: the time zone of the site in hours, positive East of the
  Greenwich meridian (ahead of GMT).  Use this parameter to easily
  account for Daylight Savings time (e.g. -4=EDT, -5 = EST/CDT), scalar.
  * `date`: this is the local civil time with type `DateTime`.  It can
  be a scalar or an array.

### Output ###

The local sidereal time for the date/time specified in hours.  This is a scalar
or an array of the same length as `jd` or `date`.

### Method ###

The Julian days of the day and time is question is used to determine the number
of days to have passed since 2000-01-01.  This is used in conjunction with the
GST of that date to extrapolate to the current GST; this is then used to get the
LST.  See Astronomical Algorithms by Jean Meeus, p. 84 (Eq. 11-4) for the
constants used.

### Example ###

Find the Greenwich mean sidereal time (GMST) on 2008-07-30 at 15:53 in
Baltimore, Maryland (longitude=-76.72 degrees).  The timezone is EDT or tz=-4

``` julia
lst = ct2lst(-76.72, -4, DateTime(2008, 7, 30, 15, 53))
# => 11.356505172312609
sixty(lst)
# => 3-element Array{Float64,1}:
#     11.0    # Hours
#     21.0    # Minutes
#     23.4186 # Seconds
```

Find the Greenwich mean sidereal time (GMST) on 2015-11-24 at 13:21 in
Heidelberg, Germany (longitude=08° 43' E).  The timezone is CET or tz=1.
Provide `ct2lst` only with the longitude of the place and the number of
Julian days.

``` julia
# Convert longitude to decimals.
longitude=ten(8, 43);
# Get number of Julian days.  Remember to subtract the time zone in
# order to convert local time to UTC.
jd = jdcnv(DateTime(2015, 11, 24, 13, 21) - Dates.Hour(1));
# Calculate Greenwich Mean Sidereal Time.
lst = ct2lst(longitude, jd)
# => 17.140685171005316
sixty(lst)
# => 3-element Array{Float64,1}:
#     17.0    # Hours
#     8.0    # Minutes
#     26.4666 # Seconds
```

### Notes ###

Code of this function is based on IDL Astronomy User's Library.
"""
ct2lst(long::Real, jd::Real) = _ct2lst(promote(float(long), float(jd))...)

function ct2lst{T<:AbstractFloat}(long::T, tz::T,
                                  date::DateTime)
    # In order to handle time zones, package "TimeZones.jl" is much better, bur
    # here we need only to add the time zone to UTC time.  All time zones I know
    # are either integer, or ±30 or ±45 minutes, so it should be safe enough to
    # convert hours to minutes and subtract minutes from "time".
    date = date - Dates.Minute(round(Int64, tz*60))
    return ct2lst(long, jdcnv(date))
end

ct2lst(long::Real, tz::Real, date::DateTime) =
    ct2lst(promote(float(long), float(tz))..., date)
