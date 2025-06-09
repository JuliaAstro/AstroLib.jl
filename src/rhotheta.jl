# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

function _rhotheta{T<:AbstractFloat}(period::T, periastron::T, eccentricity::T,
                                     semimajor_axis::T, inclination::T,
                                     omega::T, omega2::T, epoch::T)
    rho = theta = -one(period)
    # See chapter 55.
    n = 360.0/period
    M = deg2rad(n*(epoch - periastron))
    E = kepler_solver(M, eccentricity)
    r  = semimajor_axis*(1.0 - eccentricity*cos(E))
    nu = trueanom(E, eccentricity)
    # Convert variables in radians.
    omega2      = deg2rad(omega2)
    inclination = deg2rad(inclination)
    omega       = deg2rad(omega)
    theta = omega + atan2(sin(nu + omega2)*cos(inclination), cos(nu + omega2))
    rho   = r*cos(nu + omega2)/cos(theta - omega)
    # Convert theta to degrees and for it to be in [0, 360) range.
    theta = cirrange(rad2deg(theta))
    return rho, theta
end

"""
    rhotheta(period, periastron, eccentricity, semimajor_axis, inclination, omega, omega2, epoch) -> rho, theta

### Purpose ###

Calculate the separation and position angle of a binary star.

### Explanation ###

This function will return the separation \$\\rho\$ and position angle
\$\\theta\$ of a visual binary star derived from its orbital elements.  The
algorithms described in the following book will be used: Meeus J., 1992,
Astronomische Algorithmen, Barth.  Compared to the examples given at page 400
and no discrepancy found.

### Arguments ###

* `period`: period [year]
* `periastro`: time of periastron passage [year]
* `eccentricity`: eccentricity of the orbit
* `semimajor_axis`: semi-major axis [arc second]
* `inclination`: inclination angle [degree]
* `omega`: node [degree]
* `omega2`: longitude of periastron [degree]
* `epoch`: epoch of observation [year]

All input parameters have to be scalars.

### Output ###

The 2-tuple \$(\\rho, \\theta)\$, where

* ``\\rho`` is separation [arc second], and
* ``\\theta`` is position angle (degree).

### Example ###

Find the position of Eta Coronae Borealis at the epoch 2016

``` julia
ρ, θ = rhotheta(41.623, 1934.008, 0.2763, 0.907, 59.025, 23.717, 219.907, 2016)
# => (0.6351167848228113,214.42513388052114)
```

### Notes ###

Code of this function is based on IDL Astronomy User's Library.
"""
rhotheta(period::Real, periastron::Real, eccentricity::Real,
         semimajor_axis::Real, inclination::Real,
         omega::Real, omega2::Real, epoch::Real) =
             _rhotheta(promote(float(period), float(periastron),
                               float(eccentricity), float(semimajor_axis),
                               float(inclination), float(omega),
                               float(omega2), float(epoch))...)
