# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function _helio_jd{T<:AbstractFloat}(date::T, ra::T, dec::T, B1950::Bool, diff::Bool)
    # Because `xyz' uses default B1950 coordinates, we'll convert everything to
    # B1950.
    if ! B1950
        ra, dec = bprecess(ra, dec)
    end
    delta_t = (date - 33282.42345905)*inv(JULIANYEAR*100)
    epsilon_sec = @evalpoly(delta_t, 44.836, -46.8495, -0.00429, 0.00181)
    epsilon = deg2rad(23.433333 + epsilon_sec/3600.0)
    ra = deg2rad(ra)
    dec = deg2rad(dec)
    x, y, z = xyz(date)
    time = -499.00522*(cos(dec)*cos(ra)*x +
                       (tan(epsilon)*sin(dec) + cos(dec)*sin(ra))*y)
    if diff
        return time
    else
        return date + time/86400.0
    end
end

"""
    helio_jd(date, ra, dec[, B1950=true, diff=false]) -> jd_helio
    helio_jd(date, ra, dec[, B1950=true, diff=true]) -> time_diff

### Purpose ###

Convert geocentric (reduced) Julian date to heliocentric Julian date.

### Explanation ###

This procedure corrects for the extra light travel time between the Earth and the
Sun.

An online calculator for this quantity is available at
http://www.physics.sfasu.edu/astro/javascript/hjd.html

Users requiring more precise calculations and documentation should look at the
IDL code available at http://astroutils.astronomy.ohio-state.edu/time/

### Arguments ###

* `date`: reduced Julian date (= JD - 2400000), it can be either a scalar or
  vector.  You can use `juldate()` to calculate the reduced Julian date.
* `ra` and `dec`: scalars giving right ascension and declination in degrees.
  Default equinox is J2000.
* `B1950` (optional boolean keyword): if set to `true`, then input coordinates
  are assumed to be in equinox B1950 coordinates.  Default is `false`.
* `diff` (optional boolean keyword): if set to `true`, the function returns the
  time difference (heliocentric JD - geocentric JD) in seconds.  Default is
  `false`.

### Output ###

The return value depends on the value of `diff` optional keywords:

* if `diff` is `false` (default), then the heliocentric reduced Julian date is
  returned.
* if `diff` is `true`, then the time difference in seconds between the
  geocentric and heliocentric Julian date is returned.

### Example ###

What is the heliocentric Julian date of an observation of V402 Cygni (J2000: RA
= 20 9 7.8, Dec = 37 09 07) taken on June 15, 2016 at 11:40 UT?

``` julia
jd = juldate(2016, 6, 15, 11, 40);
helio_jd(jd, ten(20, 9, 7.8)*15, ten(37, 9, 7))
# => 57554.98808289718
```

### Notes ###

Wayne Warren (Raytheon ITSS) has compared the results of this algorithm with the
FORTRAN subroutines in the STARLINK SLALIB library (see
http://star-www.rl.ac.uk/).

```
                                                 Time Diff (sec)
     Date               RA(2000)   Dec(2000)  STARLINK      IDL

1999-10-29T00:00:00.0  21 08 25.  -67 22 00.  -59.0        -59.0
1999-10-29T00:00:00.0  02 56 33.4 +00 26 55.  474.1        474.1
1940-12-11T06:55:00.0  07 34 41.9 -00 30 42.  366.3        370.2
1992-02-29T03:15:56.2  12 56 27.4 +42 10 17.  350.8        350.9
2000-03-01T10:26:31.8  14 28 36.7 -20 42 11.  243.7        243.7
2100-02-26T09:18:24.2  08 26 51.7 +85 47 28.  104.0        108.8
```

Code of this function is based on IDL Astronomy User's Library.
"""
helio_jd(date::Real, ra::Real, dec::Real;
         B1950::Bool=false, diff::Bool=false) =
             _helio_jd(promote(float(date), float(ra), float(dec))...,
                       B1950, diff)
