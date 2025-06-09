# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function _mphase{T<:AbstractFloat}(jd::T)
    ram, decm, dism = moonpos(jd, radians=true)
    ras, decs = sunpos(jd, radians=true)
    # phi: geocentric elongation of the Moon from the Sun
    # inc: selenocentric (Moon centered) elongation of the Earth from the Sun
    phi = acos(sin(decs)*sin(decm) + cos(decs)*cos(decm)*cos(ras - ram))
    # "dism" is in kilometers, AU in meters
    inc = atan2(AU*sin(phi), dism*1e3 - AU*cos(phi))
    return (1.0 + cos(inc))/2.0
end

"""
    mphase(jd) -> k

### Purpose ###

Return the illuminated fraction of the Moon at given Julian date(s).

### Arguments ###

* `jd`: the Julian ephemeris date.  It can be either a scalar or an array.

### Output ###

The illuminated fraction \$k\$ of Moon's disk, with \$0 \\leq k \\leq 1\$. \$k
= 0\$ indicates a new moon, while \$k = 1\$ stands for a full moon.  If `jd` is
given as an array, an array of the same number of elements as `jd` is returned.

### Method ###

Algorithm from Chapter 46 of "Astronomical Algorithms" by Jean Meeus
(Willmann-Bell, Richmond) 1991.  `sunpos` and `moonpos` are used to get
positions of the Sun and the Moon, and the Moon distance.  The selenocentric
elongation of the Earth from the Sun (phase angle) is then computed, and used to
determine the illuminated fraction.

### Example ###

Plot the illuminated fraction of the Moon for every day in January 2018 with a
hourly sampling.  Use [PyPlot.jl](https://github.com/stevengj/PyPlot.jl) for
plotting

``` julia
using PyPlot
points = DateTime(2018,01,01):Dates.Hour(1):DateTime(2018,01,31,23,59,59);
plot(points, mphase(jdcnv(points)))
```

Note that in this calendar month there are two full moons, this event is called
[blue moon](https://en.wikipedia.org/wiki/Blue_moon).

### Notes ###

Code of this function is based on IDL Astronomy User's Library.
"""
mphase(jd::Real) = _mphase(float(jd))
