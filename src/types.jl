# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

import Base.show

##### Observatory
"""
Type holding information about an observing site.  Its fields are:

* `name`: the name of the site
* `latitude`: North-ward latitude of the site in degrees
* `longitude`: East-ward longitude of the site in degrees
* `altitude`: altitude of the site in meters
* `tz`: the number of hours of offset from UTC
"""
struct Observatory
    name::String
    latitude::Float64
    longitude::Float64
    altitude::Float64
    tz::Float64 # There are non-integer time zones

    # Define constructor that automatically converts longitude and latitude with
    # "ten", for convenience.
    Observatory(name, lat, long, alt, tz) = new(
        String(name),
        Float64(ten(lat)), Float64(ten(long)),
        Float64(float(alt)),
        Float64(ten(tz))
    )
end

# New type representation
function show(io::IO, obs::Observatory)
    println(io, "Observatory: ", obs.name)
    println(io, "latitude:    ", obs.latitude, "°N")
    println(io, "longitude:   ", obs.longitude, "°E")
    println(io, "altitude:    ", obs.altitude, " m")
    tzdec, tzint = modf(obs.tz)
    tzstr = tzdec == 0 ? "" : @sprintf(":%d", abs(tzdec * 60))
    print(io, "time zone:   ", @sprintf("UTC%+d", tzint), tzstr)
    return
end

##### Planet
"""
Type holding information about a planet.  Its fields are:

Designation:

* `name`: the name

Physical characteristics:

* `radius`: mean radius in meters
* `eqradius`: equatorial radius in meters
* `polradius`: polar radius in meters
* `mass`: mass in kilogram

Orbital characteristics (epoch J2000):

* `ecc`: eccentricity of the orbit
* `axis`: semi-major axis of the orbit in meters
* `period`: sidereal orbital period in seconds

Position characteristics (epoch J2000):

* `inc`: inclination in degrees
* `asc_long`: longitude of the ascending node in degrees
* `per_long`: longitude of perihelion in degrees
* `mean_long`: mean longitude in degrees
"""
Base.@kwdef struct Planet
    name::String
    radius::Float64
    eqradius::Float64
    polradius::Float64
    mass::Float64
    ecc::Float64
    axis::Float64
    period::Float64
    inc::Float64
    asc_long::Float64
    per_long::Float64
    mean_long::Float64
end

# New type representation
function show(io::IO, pl::Planet)
    println(io, "Planet:                      ", uppercasefirst(pl.name))
    println(io, "mean radius:                 ", pl.radius, " m")
    println(io, "equatorial radius:           ", pl.eqradius, " m")
    println(io, "polar radius:                ", pl.polradius, " m")
    println(io, "mass:                        ", pl.mass, " kg")
    println(io, "eccentricity:                ", pl.ecc)
    println(io, "semi-major axis:             ", pl.axis, " m")
    println(io, "period:                      ", pl.period, " s")
    println(io, "inclination:                 ", pl.inc, " °")
    println(io, "longitude of ascending node: ", pl.asc_long, " °")
    println(io, "longitude of perihelion:     ", pl.per_long, " °")
    print(io,   "mean longitude:              ", pl.mean_long, " °")
    return
end

export Observatory, Planet
