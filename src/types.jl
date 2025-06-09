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
immutable Observatory
    name::AbstractString
    latitude::Real
    longitude::Real
    altitude::Real
    tz::Real # There are non-integer time zones
    # Define constructor that automatically converts longitude and latitude with
    # "ten", for convenience.
    Observatory(name, lat, long, alt, tz) =
        new(name, ten(lat), ten(long), float(alt), ten(tz))
end

# New type representation
function show(io::IO, obs::Observatory)
    println(io, "Observatory: ", obs.name)
    println(io, "latitude:    ", obs.latitude, "°N")
    println(io, "longitude:   ", obs.longitude, "°E")
    println(io, "altitude:    ", obs.altitude, " m")
    tzdec, tzint = modf(obs.tz)
    print(io,   "time zone:   ", @sprintf("UTC%+d", tzint),
          tzdec == 0 ? "" : @sprintf(":%d", abs(tzdec*60)))
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
"""
immutable Planet
    name::AbstractString
    radius::Real
    eqradius::Real
    polradius::Real
    mass::Real
    ecc::Real
    axis::Real
    period::Real
    Planet(name, radius, eqradius, polradius, mass, ecc, axis, period) =
        new(name, float(radius), float(eqradius),
            float(polradius), float(mass), float(ecc),
            float(axis), float(period))
end

# New type representation
function show(io::IO, pl::Planet)
    println(io, "Planet:            ", ucfirst(pl.name))
    println(io, "mean radius:       ", pl.radius, " m")
    println(io, "equatorial radius: ", pl.eqradius, " m")
    println(io, "polar radius:      ", pl.polradius, " m")
    println(io, "mass:              ", pl.mass, " kg")
    println(io, "eccentricity:      ", pl.ecc)
    println(io, "semi-major axis:   ", pl.axis, " m")
    print(io,   "period:            ", pl.axis, " s")
end

export Observatory, Planet
